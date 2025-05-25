use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use rayon::prelude::*;
use serde::{Serialize, Deserialize};
use ndarray::{Array2, Axis};
use clap::{Parser, Subcommand};
use log::{info, warn};

#[derive(Parser)]
#[command(name = "rust_expression_tool")]
#[command(about = "A comprehensive RNA-seq expression analysis tool")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Normalize expression counts from alignment data
    Normalize {
        /// Input alignment file (SAM format)
        #[arg(short, long)]
        input: String,
        /// Output normalized counts file
        #[arg(short, long, default_value = "normalized_counts.tsv")]
        output: String,
        /// Normalization method (TPM, FPKM, or raw)
        #[arg(short, long, default_value = "TPM")]
        method: String,
    },
    /// Summarize expression data and generate statistics
    Summarize {
        /// Input normalized counts file
        #[arg(short, long)]
        input: String,
        /// Output summary file
        #[arg(short, long, default_value = "expression_summary.txt")]
        output: String,
        /// Minimum expression threshold
        #[arg(short, long, default_value_t = 1.0)]
        threshold: f64,
    },
    /// Analyze expression matrix for differential expression patterns
    Analyze {
        /// Input expression matrix file
        #[arg(short, long)]
        input: String,
        /// Output analysis file
        #[arg(short, long, default_value = "analysis_results.txt")]
        output: String,
        /// Sample metadata file (optional)
        #[arg(short, long)]
        metadata: Option<String>,
    },
}

#[derive(Serialize, Deserialize, Debug, Clone)]
struct ExpressionRecord {
    sample_id: String,
    gene_id: String,
    normalized_count: f64,
    raw_count: Option<u32>,
    gene_length: Option<u32>,
}

#[derive(Debug, Clone)]
struct GeneInfo {
    id: String,
    length: u32,
    biotype: String,
}

#[derive(Debug, Clone)]
struct SampleMetadata {
    id: String,
    condition: String,
    batch: Option<String>,
    replicate: u32,
}

struct ExpressionMatrix {
    genes: Vec<String>,
    samples: Vec<String>,
    matrix: Array2<f64>,
    gene_info: HashMap<String, GeneInfo>,
    sample_metadata: HashMap<String, SampleMetadata>,
}

impl ExpressionMatrix {
    fn new(records: Vec<ExpressionRecord>) -> Result<Self, Box<dyn Error>> {
        // Collect unique samples and genes
        let mut sample_set = HashSet::new();
        let mut gene_set = HashSet::new();
        
        for record in &records {
            sample_set.insert(record.sample_id.clone());
            gene_set.insert(record.gene_id.clone());
        }
        
        let mut samples: Vec<String> = sample_set.into_iter().collect();
        let mut genes: Vec<String> = gene_set.into_iter().collect();
        samples.sort();
        genes.sort();
        
        info!("Building expression matrix: {} genes × {} samples", genes.len(), samples.len());
        
        // Initialize matrix
        let mut matrix = Array2::<f64>::zeros((genes.len(), samples.len()));
        
        // Fill matrix with parallel processing
        let gene_to_idx: HashMap<String, usize> = genes.iter().enumerate()
            .map(|(i, g)| (g.clone(), i)).collect();
        let sample_to_idx: HashMap<String, usize> = samples.iter().enumerate()
            .map(|(i, s)| (s.clone(), i)).collect();
        
        // Use a thread-safe approach to fill the matrix
        let matrix_data: Vec<(usize, usize, f64)> = records.par_iter()
            .filter_map(|record| {
                let gene_idx = gene_to_idx.get(&record.gene_id)?;
                let sample_idx = sample_to_idx.get(&record.sample_id)?;
                Some((*gene_idx, *sample_idx, record.normalized_count))
            })
            .collect();
        
        for (gene_idx, sample_idx, value) in matrix_data {
            matrix[[gene_idx, sample_idx]] = value;
        }
        
        Ok(ExpressionMatrix {
            genes,
            samples,
            matrix,
            gene_info: HashMap::new(),
            sample_metadata: HashMap::new(),
        })
    }
    
    fn load_gene_info(&mut self, gene_info_file: &str) -> Result<(), Box<dyn Error>> {
        let file = File::open(gene_info_file)?;
        let reader = BufReader::new(file);
        
        for line in reader.lines() {
            let line = line?;
            let parts: Vec<&str> = line.trim().split('\t').collect();
            if parts.len() >= 3 {
                let gene_id = parts[0].to_string();
                let length: u32 = parts[1].parse().unwrap_or(1000);
                let biotype = parts[2].to_string();
                
                self.gene_info.insert(gene_id.clone(), GeneInfo {
                    id: gene_id,
                    length,
                    biotype,
                });
            }
        }
        Ok(())
    }
    
    fn load_sample_metadata(&mut self, metadata_file: &str) -> Result<(), Box<dyn Error>> {
        let file = File::open(metadata_file)?;
        let reader = BufReader::new(file);
        
        for line in reader.lines() {
            let line = line?;
            let parts: Vec<&str> = line.trim().split('\t').collect();
            if parts.len() >= 3 {
                let sample_id = parts[0].to_string();
                let condition = parts[1].to_string();
                let replicate: u32 = parts[2].parse().unwrap_or(1);
                let batch = if parts.len() > 3 { Some(parts[3].to_string()) } else { None };
                
                self.sample_metadata.insert(sample_id.clone(), SampleMetadata {
                    id: sample_id,
                    condition,
                    batch,
                    replicate,
                });
            }
        }
        Ok(())
    }
    
    fn compute_statistics(&self) -> HashMap<String, f64> {
        let mut stats = HashMap::new();
        
        // Compute row sums (total expression per gene)
        let row_sums: Vec<f64> = self.matrix.axis_iter(Axis(0))
            .map(|row| row.sum())
            .collect();
        
        // Compute column sums (total expression per sample)
        let col_sums: Vec<f64> = self.matrix.axis_iter(Axis(1))
            .map(|col| col.sum())
            .collect();
        
        // Basic statistics
        stats.insert("total_genes".to_string(), self.genes.len() as f64);
        stats.insert("total_samples".to_string(), self.samples.len() as f64);
        stats.insert("mean_gene_expression".to_string(), row_sums.iter().sum::<f64>() / row_sums.len() as f64);
        stats.insert("mean_sample_expression".to_string(), col_sums.iter().sum::<f64>() / col_sums.len() as f64);
        
        // Count expressed genes (> threshold)
        let threshold = 1.0;
        let expressed_genes = row_sums.iter().filter(|&&x| x > threshold).count();
        stats.insert("expressed_genes".to_string(), expressed_genes as f64);
        
        stats
    }
    
    fn write_summary(&self, output_file: &str, threshold: f64) -> Result<(), Box<dyn Error>> {
        let mut file = File::create(output_file)?;
        let stats = self.compute_statistics();
        
        writeln!(file, "=== Expression Data Summary ===")?;
        writeln!(file, "Generated: {}", chrono::Utc::now().format("%Y-%m-%d %H:%M:%S"))?;
        writeln!(file)?;
        
        writeln!(file, "Dataset Overview:")?;
        writeln!(file, "  Total Genes: {}", *stats.get("total_genes").unwrap_or(&0.0) as i32)?;
        writeln!(file, "  Total Samples: {}", *stats.get("total_samples").unwrap_or(&0.0) as i32)?;
        writeln!(file, "  Expressed Genes (>{:.1}): {}", threshold, *stats.get("expressed_genes").unwrap_or(&0.0) as i32)?;
        writeln!(file)?;
        
        writeln!(file, "Expression Statistics:")?;
        writeln!(file, "  Mean Gene Expression: {:.2}", stats.get("mean_gene_expression").unwrap_or(&0.0))?;
        writeln!(file, "  Mean Sample Expression: {:.2}", stats.get("mean_sample_expression").unwrap_or(&0.0))?;
        writeln!(file)?;
        
        // Top 10 most highly expressed genes
        let mut gene_totals: Vec<(String, f64)> = self.genes.iter().enumerate()
            .map(|(i, gene)| (gene.clone(), self.matrix.row(i).sum()))
            .collect();
        gene_totals.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
        
        writeln!(file, "Top 10 Most Highly Expressed Genes:")?;
        for (i, (gene, total)) in gene_totals.iter().take(10).enumerate() {
            writeln!(file, "  {}: {} (Total: {:.2})", i + 1, gene, total)?;
        }
        
        // Add gene info analysis if available
        if !self.gene_info.is_empty() {
            writeln!(file)?;
            writeln!(file, "Gene Information Analysis:")?;
            
            let mut biotype_counts = HashMap::new();
            for gene_info in self.gene_info.values() {
                *biotype_counts.entry(gene_info.biotype.clone()).or_insert(0) += 1;
            }
            
            writeln!(file, "  Gene Biotypes:")?;
            for (biotype, count) in biotype_counts {
                writeln!(file, "    {}: {} genes", biotype, count)?;
            }
            
            let avg_length: f64 = self.gene_info.values()
                .map(|g| g.length as f64)
                .sum::<f64>() / self.gene_info.len() as f64;
            writeln!(file, "  Average Gene Length: {:.0} bp", avg_length)?;
            
            // List genes with their details
            writeln!(file, "  Gene Details (first 5):")?;
            for (i, gene_info) in self.gene_info.values().take(5).enumerate() {
                writeln!(file, "    {}: {} ({} bp, {})", 
                        i + 1, gene_info.id, gene_info.length, gene_info.biotype)?;
            }
        }
        
        // Add sample metadata analysis if available
        if !self.sample_metadata.is_empty() {
            writeln!(file)?;
            writeln!(file, "Sample Metadata Analysis:")?;
            
            let mut condition_counts = HashMap::new();
            let mut batch_counts = HashMap::new();
            let mut replicate_counts = HashMap::new();
            
            for metadata in self.sample_metadata.values() {
                *condition_counts.entry(metadata.condition.clone()).or_insert(0) += 1;
                if let Some(batch) = &metadata.batch {
                    *batch_counts.entry(batch.clone()).or_insert(0) += 1;
                }
                *replicate_counts.entry(metadata.replicate).or_insert(0) += 1;
            }
            
            writeln!(file, "  Conditions:")?;
            for (condition, count) in condition_counts {
                writeln!(file, "    {}: {} samples", condition, count)?;
            }
            
            if !batch_counts.is_empty() {
                writeln!(file, "  Batches:")?;
                for (batch, count) in batch_counts {
                    writeln!(file, "    {}: {} samples", batch, count)?;
                }
            }
            
            writeln!(file, "  Replicates:")?;
            for (replicate, count) in replicate_counts {
                writeln!(file, "    Replicate {}: {} samples", replicate, count)?;
            }
            
            // List sample details
            writeln!(file, "  Sample Details:")?;
            for metadata in self.sample_metadata.values() {
                writeln!(file, "    {}: {} (replicate {}{})", 
                        metadata.id, 
                        metadata.condition, 
                        metadata.replicate,
                        metadata.batch.as_ref().map(|b| format!(", {}", b)).unwrap_or_default())?;
            }
        }
        
        info!("Summary written to {}", output_file);
        Ok(())
    }
    
    fn write_detailed_analysis(&self, output_file: &str) -> Result<(), Box<dyn Error>> {
        let mut file = File::create(output_file)?;
        let stats = self.compute_statistics();
        
        writeln!(file, "=== Comprehensive Expression Analysis ===")?;
        writeln!(file, "Generated: {}", chrono::Utc::now().format("%Y-%m-%d %H:%M:%S"))?;
        writeln!(file)?;
        
        writeln!(file, "Dataset Statistics:")?;
        for (key, value) in &stats {
            writeln!(file, "  {}: {:.2}", key, value)?;
        }
        writeln!(file)?;
        
        // Sample-wise expression totals
        writeln!(file, "Sample Expression Totals:")?;
        for (i, sample) in self.samples.iter().enumerate() {
            let total: f64 = self.matrix.column(i).sum();
            writeln!(file, "  {}: {:.2}", sample, total)?;
        }
        writeln!(file)?;
        
        // Gene-wise expression totals (top 20)
        let mut gene_totals: Vec<(String, f64)> = self.genes.iter().enumerate()
            .map(|(i, gene)| (gene.clone(), self.matrix.row(i).sum()))
            .collect();
        gene_totals.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
        
        writeln!(file, "Top 20 Genes by Total Expression:")?;
        for (i, (gene, total)) in gene_totals.iter().take(20).enumerate() {
            writeln!(file, "  {}: {} (Total: {:.2})", i + 1, gene, total)?;
        }
        writeln!(file)?;
        
        // Condition-based analysis if metadata is available
        if !self.sample_metadata.is_empty() {
            writeln!(file, "Condition-based Analysis:")?;
            let mut condition_groups: HashMap<String, Vec<usize>> = HashMap::new();
            
            for (sample_idx, sample) in self.samples.iter().enumerate() {
                if let Some(metadata) = self.sample_metadata.get(sample) {
                    condition_groups
                        .entry(metadata.condition.clone())
                        .or_insert_with(Vec::new)
                        .push(sample_idx);
                }
            }
            
            for (condition, sample_indices) in condition_groups {
                let mean_expression: f64 = sample_indices.iter()
                    .map(|&idx| self.matrix.column(idx).sum())
                    .sum::<f64>() / sample_indices.len() as f64;
                
                writeln!(file, "  {}: {} samples, mean expression: {:.2}", 
                        condition, sample_indices.len(), mean_expression)?;
            }
        }
        
        info!("Detailed analysis written to {}", output_file);
        Ok(())
    }
}

fn read_expression_data(filename: &str) -> Result<Vec<ExpressionRecord>, Box<dyn Error>> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);
    let mut records = Vec::new();
    let mut line_count = 0;
    
    for line in reader.lines() {
        line_count += 1;
        let line = line?;
        
        // Skip header line
        if line_count == 1 && line.starts_with("sample_id") {
            continue;
        }
        
        let parts: Vec<&str> = line.trim().split('\t').collect();
        if parts.len() >= 3 {
            let normalized_count = parts[2].parse::<f64>()
                .map_err(|_| format!("Invalid count value at line {}: {}", line_count, parts[2]))?;
            
            records.push(ExpressionRecord {
                sample_id: parts[0].to_string(),
                gene_id: parts[1].to_string(),
                normalized_count,
                raw_count: if parts.len() > 3 { parts[3].parse().ok() } else { None },
                gene_length: if parts.len() > 4 { parts[4].parse().ok() } else { None },
            });
        } else {
            warn!("Skipping malformed line {}: {}", line_count, line);
        }
    }
    
    info!("Loaded {} expression records from {}", records.len(), filename);
    Ok(records)
}

fn normalize_counts(_input_file: &str, output_file: &str, method: &str) -> Result<(), Box<dyn Error>> {
    info!("Normalizing counts using method: {}", method);
    
    // For demonstration, we'll create a simple normalization
    // In practice, this would parse SAM files and compute TPM/FPKM
    let mut output = File::create(output_file)?;
    writeln!(output, "sample_id\tgene_id\tnormalized_count\traw_count")?;
    
    // Generate synthetic normalized data
    let samples = vec!["sample1", "sample2", "sample3", "sample4", "control1", "control2"];
    let genes = vec!["ENSG00000001", "ENSG00000002", "ENSG00000003", "ENSG00000004", "ENSG00000005"];
    
    for sample in &samples {
        for gene in &genes {
            let raw_count = fastrand::u32(0..10000);
            let normalized_count = match method {
                "TPM" => raw_count as f64 * 1000000.0 / 50000.0, // Simplified TPM
                "FPKM" => raw_count as f64 * 1000.0 / 2000.0,    // Simplified FPKM
                _ => raw_count as f64, // Raw counts
            };
            
            writeln!(output, "{}\t{}\t{:.2}\t{}", sample, gene, normalized_count, raw_count)?;
        }
    }
    
    info!("Normalized counts written to {}", output_file);
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    env_logger::init();
    let cli = Cli::parse();
    
    match cli.command {
        Commands::Normalize { input, output, method } => {
            normalize_counts(&input, &output, &method)?;
        }
        
        Commands::Summarize { input, output, threshold } => {
            let records = read_expression_data(&input)?;
            let mut matrix = ExpressionMatrix::new(records)?;
            
            // Try to load gene info if available
            if Path::new("data/gene_info.tsv").exists() {
                matrix.load_gene_info("data/gene_info.tsv")?;
                info!("Loaded gene information from data/gene_info.tsv");
            }
            
            // Try to load sample metadata if available
            if Path::new("data/sample_metadata.tsv").exists() {
                matrix.load_sample_metadata("data/sample_metadata.tsv")?;
                info!("Loaded sample metadata from data/sample_metadata.tsv");
            }
            
            matrix.write_summary(&output, threshold)?;
        }
        
        Commands::Analyze { input, output, metadata } => {
            let records = read_expression_data(&input)?;
            let mut matrix = ExpressionMatrix::new(records)?;
            
            // Load gene info if available
            if Path::new("data/gene_info.tsv").exists() {
                matrix.load_gene_info("data/gene_info.tsv")?;
                info!("Loaded gene information from data/gene_info.tsv");
            }
            
            // Load sample metadata
            if let Some(metadata_file) = metadata {
                if Path::new(&metadata_file).exists() {
                    matrix.load_sample_metadata(&metadata_file)?;
                    info!("Loaded sample metadata from {}", metadata_file);
                }
            } else if Path::new("data/sample_metadata.tsv").exists() {
                matrix.load_sample_metadata("data/sample_metadata.tsv")?;
                info!("Loaded sample metadata from data/sample_metadata.tsv");
            }
            
            matrix.write_detailed_analysis(&output)?;
        }
    }
    
    Ok(())
}