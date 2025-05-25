use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::collections::HashMap;
use std::path::Path;
use clap::{App, Arg};
use log::{info, warn};
use serde::{Serialize, Deserialize};

#[derive(Debug, Serialize, Deserialize, Clone)]
struct SparseEntry {
    gene_idx: usize,
    cell_idx: usize,
    count: f64,
}

#[derive(Debug, Serialize, Deserialize)]
struct CellCoordinate {
    cell_id: usize,
    pc1: f64,
    pc2: f64,
    pc3: f64,
}

#[derive(Debug)]
struct AnalysisConfig {
    input_file: String,
    output_file: String,
    num_components: usize,
    // Note: These fields are kept for future implementation
    #[allow(dead_code)]
    min_genes_per_cell: usize,
    #[allow(dead_code)]
    min_cells_per_gene: usize,
}

fn main() -> Result<(), Box<dyn Error>> {
    // Initialize logger
    env_logger::init();
    info!("🧬 Single-Cell RNA-seq Analyzer v1.0");

    // Parse command line arguments
    let matches = App::new("Single-Cell RNA-seq Analyzer")
        .version("1.0")
        .author("Bioinformatics Team")
        .about("Performs dimensionality reduction on single-cell RNA-seq data")
        .arg(Arg::with_name("input")
            .short("i")
            .long("input")
            .takes_value(true)
            .help("Input sparse matrix file (TSV format)")
            .default_value("data/raw/sparse_counts.tsv"))
        .arg(Arg::with_name("output")
            .short("o")
            .long("output")
            .takes_value(true)
            .help("Output coordinates file")
            .default_value("data/processed/cell_coords.tsv"))
        .arg(Arg::with_name("components")
            .short("c")
            .long("components")
            .takes_value(true)
            .help("Number of principal components")
            .default_value("3"))
        .arg(Arg::with_name("min-genes-per-cell")
            .long("min-genes-per-cell")
            .takes_value(true)
            .help("Minimum genes per cell for filtering")
            .default_value("200"))
        .arg(Arg::with_name("min-cells-per-gene")
            .long("min-cells-per-gene")
            .takes_value(true)
            .help("Minimum cells per gene for filtering")
            .default_value("3"))
        .get_matches();

    let config = AnalysisConfig {
        input_file: matches.value_of("input").unwrap().to_string(),
        output_file: matches.value_of("output").unwrap().to_string(),
        num_components: matches.value_of("components").unwrap().parse()?,
        min_genes_per_cell: matches.value_of("min-genes-per-cell").unwrap().parse()?,
        min_cells_per_gene: matches.value_of("min-cells-per-gene").unwrap().parse()?,
    };

    info!("📂 Input file: {}", config.input_file);
    info!("📁 Output file: {}", config.output_file);

    // Load sparse matrix data
    info!("📊 Loading sparse matrix data...");
    let entries = load_sparse_data(&config.input_file)?;
    info!("✅ Loaded {} entries", entries.len());

    // Calculate matrix dimensions
    let max_gene = entries.iter().map(|e| e.gene_idx).max().unwrap_or(0);
    let max_cell = entries.iter().map(|e| e.cell_idx).max().unwrap_or(0);
    info!("🔢 Matrix dimensions: {} genes × {} cells", max_gene + 1, max_cell + 1);

    // Calculate quality metrics
    info!("📈 Calculating quality metrics...");
    let qc_metrics = calculate_quality_metrics(&entries, max_gene + 1, max_cell + 1);
    print_quality_metrics(&qc_metrics);

    // Apply basic filtering (currently returns all entries)
    info!("🔍 Applying quality filters...");
    let filtered_entries = apply_basic_filters(&entries);
    info!("✅ After filtering: {} entries", filtered_entries.len());

    // Perform mock dimensionality reduction
    info!("🎯 Performing dimensionality reduction...");
    let coordinates = perform_mock_pca(&filtered_entries, max_cell + 1, config.num_components)?;
    info!("✅ Generated coordinates for {} cells", coordinates.len());

    // Save results
    info!("💾 Saving results...");
    save_coordinates(&coordinates, &config.output_file)?;
    info!("✅ Results saved to: {}", config.output_file);

    // Print summary
    print_analysis_summary(&coordinates);
    
    info!("🎉 Analysis completed successfully!");
    Ok(())
}

fn load_sparse_data(file_path: &str) -> Result<Vec<SparseEntry>, Box<dyn Error>> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);
    let mut entries = Vec::new();
    let mut line_count = 0;

    for line in reader.lines() {
        line_count += 1;
        let line = line?;
        
        // Skip header line
        if line_count == 1 && line.contains("gene_idx") {
            continue;
        }

        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() >= 3 {
            match (parts[0].parse::<usize>(), parts[1].parse::<usize>(), parts[2].parse::<f64>()) {
                (Ok(gene_idx), Ok(cell_idx), Ok(count)) => {
                    if count > 0.0 {
                        entries.push(SparseEntry { gene_idx, cell_idx, count });
                    }
                }
                _ => warn!("Skipping invalid line {}: {}", line_count, line),
            }
        }
    }

    if entries.is_empty() {
        return Err("No valid entries found in input file".into());
    }

    Ok(entries)
}

fn calculate_quality_metrics(entries: &[SparseEntry], n_genes: usize, n_cells: usize) -> HashMap<String, f64> {
    let mut metrics = HashMap::new();
    
    // Calculate basic statistics
    let total_counts: f64 = entries.iter().map(|e| e.count).sum();
    let total_entries = entries.len();
    
    // Calculate per-cell statistics
    let mut cell_counts: HashMap<usize, f64> = HashMap::new();
    let mut cell_genes: HashMap<usize, usize> = HashMap::new();
    
    for entry in entries {
        *cell_counts.entry(entry.cell_idx).or_insert(0.0) += entry.count;
        *cell_genes.entry(entry.cell_idx).or_insert(0) += 1;
    }
    
    let mean_counts_per_cell = if !cell_counts.is_empty() {
        cell_counts.values().sum::<f64>() / cell_counts.len() as f64
    } else { 0.0 };
    
    let mean_genes_per_cell = if !cell_genes.is_empty() {
        cell_genes.values().sum::<usize>() as f64 / cell_genes.len() as f64
    } else { 0.0 };
    
    // Calculate sparsity
    let matrix_size = n_genes * n_cells;
    let sparsity = 1.0 - (total_entries as f64 / matrix_size as f64);
    
    metrics.insert("total_counts".to_string(), total_counts);
    metrics.insert("total_entries".to_string(), total_entries as f64);
    metrics.insert("unique_cells".to_string(), cell_counts.len() as f64);
    metrics.insert("unique_genes".to_string(), 
        entries.iter().map(|e| e.gene_idx).collect::<std::collections::HashSet<_>>().len() as f64);
    metrics.insert("mean_counts_per_cell".to_string(), mean_counts_per_cell);
    metrics.insert("mean_genes_per_cell".to_string(), mean_genes_per_cell);
    metrics.insert("sparsity".to_string(), sparsity);
    
    metrics
}

fn print_quality_metrics(metrics: &HashMap<String, f64>) {
    info!("=== QUALITY CONTROL METRICS ===");
    info!("Total entries: {:.0}", metrics.get("total_entries").unwrap_or(&0.0));
    info!("Unique cells: {:.0}", metrics.get("unique_cells").unwrap_or(&0.0));
    info!("Unique genes: {:.0}", metrics.get("unique_genes").unwrap_or(&0.0));
    info!("Total counts: {:.0}", metrics.get("total_counts").unwrap_or(&0.0));
    info!("Mean counts per cell: {:.2}", metrics.get("mean_counts_per_cell").unwrap_or(&0.0));
    info!("Mean genes per cell: {:.1}", metrics.get("mean_genes_per_cell").unwrap_or(&0.0));
    info!("Matrix sparsity: {:.3}", metrics.get("sparsity").unwrap_or(&0.0));
}

fn apply_basic_filters(entries: &[SparseEntry]) -> Vec<SparseEntry> {
    // For now, return all entries (implement proper filtering later)
    // In a full implementation, you would:
    // 1. Calculate genes per cell and cells per gene
    // 2. Filter based on min_genes_per_cell and min_cells_per_gene
    // 3. Remove low-quality cells and rarely expressed genes
    
    entries.to_vec()
}

fn perform_mock_pca(entries: &[SparseEntry], n_cells: usize, num_components: usize) -> Result<Vec<CellCoordinate>, Box<dyn Error>> {
    info!("🎲 Performing mock PCA (replace with real PCA implementation)");
    
    // Calculate per-cell statistics for mock PCA coordinates
    let mut cell_stats: HashMap<usize, (f64, usize)> = HashMap::new();
    
    for entry in entries {
        let stats = cell_stats.entry(entry.cell_idx).or_insert((0.0, 0));
        stats.0 += entry.count;  // Total counts
        stats.1 += 1;            // Number of genes
    }
    
    let mut coordinates = Vec::new();
    
    for cell_id in 0..n_cells {
        let (total_counts, n_genes) = cell_stats.get(&cell_id).unwrap_or(&(0.0, 0));
        
        // Mock PCA coordinates based on cell statistics
        // In real implementation, this would be proper SVD-based PCA
        let pc1 = if *total_counts > 0.0 { 
            (total_counts / 1000.0).ln_1p() - 5.0 
        } else { -5.0 };
        
        let pc2 = (*n_genes as f64 / 100.0) - 5.0;
        
        let pc3 = if num_components > 2 {
            (cell_id as f64 * 0.01) % 10.0 - 5.0
        } else { 0.0 };
        
        coordinates.push(CellCoordinate {
            cell_id,
            pc1,
            pc2,
            pc3,
        });
    }
    
    Ok(coordinates)
}

fn save_coordinates(coordinates: &[CellCoordinate], output_file: &str) -> Result<(), Box<dyn Error>> {
    // Create output directory if it doesn't exist
    if let Some(parent) = Path::new(output_file).parent() {
        std::fs::create_dir_all(parent)?;
    }
    
    let file = File::create(output_file)?;
    let mut writer = BufWriter::new(file);
    
    // Write header
    writeln!(writer, "cell_id\tpc1\tpc2\tpc3")?;
    
    // Write coordinates
    for coord in coordinates {
        writeln!(writer, "{}\t{:.6}\t{:.6}\t{:.6}", 
                coord.cell_id, coord.pc1, coord.pc2, coord.pc3)?;
    }
    
    writer.flush()?;
    Ok(())
}

fn print_analysis_summary(coordinates: &[CellCoordinate]) {
    if coordinates.is_empty() {
        return;
    }
    
    let pc1_values: Vec<f64> = coordinates.iter().map(|c| c.pc1).collect();
    let pc2_values: Vec<f64> = coordinates.iter().map(|c| c.pc2).collect();
    let pc3_values: Vec<f64> = coordinates.iter().map(|c| c.pc3).collect();
    
    let pc1_mean = pc1_values.iter().sum::<f64>() / pc1_values.len() as f64;
    let pc2_mean = pc2_values.iter().sum::<f64>() / pc2_values.len() as f64;
    let pc3_mean = pc3_values.iter().sum::<f64>() / pc3_values.len() as f64;
    
    info!("=== ANALYSIS SUMMARY ===");
    info!("Cells analyzed: {}", coordinates.len());
    info!("PC1 mean: {:.4}, range: [{:.4}, {:.4}]", 
          pc1_mean,
          pc1_values.iter().fold(f64::INFINITY, |a, &b| a.min(b)),
          pc1_values.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b)));
    info!("PC2 mean: {:.4}, range: [{:.4}, {:.4}]", 
          pc2_mean,
          pc2_values.iter().fold(f64::INFINITY, |a, &b| a.min(b)),
          pc2_values.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b)));
    info!("PC3 mean: {:.4}, range: [{:.4}, {:.4}]", 
          pc3_mean,
          pc3_values.iter().fold(f64::INFINITY, |a, &b| a.min(b)),
          pc3_values.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b)));
}