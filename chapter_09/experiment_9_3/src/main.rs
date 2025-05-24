use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use clap::{Arg, Command};
use serde::{Serialize, Deserialize};
use ndarray::{Array2, s};
use log::{info, warn, error};
use env_logger;

#[derive(Serialize, Deserialize, Debug, Clone)]
struct GeneCount {
    gene_id: String,
    sample_id: String,
    count: f64,
}

#[derive(Debug)]
struct NormalizationStats {
    total_genes: usize,
    total_samples: usize,
    zero_counts: usize,
    size_factors: Vec<f64>,
    geometric_means_computed: usize,
}

fn main() -> Result<(), Box<dyn Error>> {
    env_logger::init();
    
    let matches = Command::new("RNA-seq Normalizer")
        .version("1.0")
        .author("Bioinformatics Pipeline")
        .about("Performs DESeq2-style normalization on RNA-seq count data")
        .arg(Arg::new("input")
            .short('i')
            .long("input")
            .value_name("FILE")
            .help("Input TSV file with gene counts")
            .default_value("raw_counts.tsv"))
        .arg(Arg::new("output")
            .short('o')
            .long("output")
            .value_name("FILE")
            .help("Output file for normalized counts")
            .default_value("normalized_counts.tsv"))
        .arg(Arg::new("stats")
            .short('s')
            .long("stats")
            .value_name("FILE")
            .help("Output file for normalization statistics")
            .default_value("normalization_stats.txt"))
        .arg(Arg::new("min_count")
            .long("min-count")
            .value_name("NUMBER")
            .help("Minimum count threshold for filtering")
            .default_value("1.0"))
        .arg(Arg::new("pseudocount")
            .long("pseudocount")
            .value_name("NUMBER")
            .help("Pseudocount to add for geometric mean calculation")
            .default_value("1.0"))
        .get_matches();

    let input_file = matches.get_one::<String>("input").unwrap();
    let output_file = matches.get_one::<String>("output").unwrap();
    let stats_file = matches.get_one::<String>("stats").unwrap();
    let min_count: f64 = matches.get_one::<String>("min_count").unwrap().parse()?;
    let pseudocount: f64 = matches.get_one::<String>("pseudocount").unwrap().parse()?;

    info!("Starting RNA-seq normalization pipeline");
    info!("Input file: {}", input_file);
    info!("Output file: {}", output_file);
    info!("Min count threshold: {}", min_count);
    info!("Pseudocount: {}", pseudocount);

    // Read and validate input data
    let data = read_count_data(input_file)?;
    info!("Read {} count entries from input file", data.len());

    if data.is_empty() {
        error!("No valid data found in input file");
        return Err("Empty dataset".into());
    }

    // Create count matrix
    let (matrix, genes, samples) = create_count_matrix(&data)?;
    info!("Created count matrix: {} genes x {} samples", genes.len(), samples.len());

    // Quality control checks
    perform_quality_checks(&matrix, &genes, &samples, min_count)?;

    // Calculate normalization factors
    let size_factors = calculate_size_factors(&matrix, &genes, &samples, pseudocount)?;
    info!("Calculated size factors for {} samples", size_factors.len());

    // Normalize the matrix
    let normalized_matrix = normalize_counts(&matrix, &size_factors)?;
    info!("Normalization completed");

    // Write results
    write_normalized_counts(&normalized_matrix, &genes, &samples, output_file)?;
    info!("Normalized counts written to {}", output_file);

    // Write statistics
    let stats = NormalizationStats {
        total_genes: genes.len(),
        total_samples: samples.len(),
        zero_counts: count_zeros(&matrix),
        size_factors: size_factors.clone(),
        geometric_means_computed: genes.len(),
    };
    write_statistics(&stats, &genes, &samples, stats_file)?;
    info!("Statistics written to {}", stats_file);

    info!("Pipeline completed successfully");
    Ok(())
}

fn read_count_data(filename: &str) -> Result<Vec<GeneCount>, Box<dyn Error>> {
    if !Path::new(filename).exists() {
        return Err(format!("Input file '{}' does not exist", filename).into());
    }

    let file = File::open(filename)?;
    let reader = BufReader::new(file);
    let mut data = Vec::new();
    let mut line_number = 0;
    let mut header_skipped = false;

    for line in reader.lines() {
        line_number += 1;
        let l = line?;
        
        // Skip empty lines
        if l.trim().is_empty() {
            continue;
        }
        
        // Skip header line if it contains non-numeric data in third column
        if !header_skipped {
            let parts: Vec<&str> = l.split('\t').collect();
            if parts.len() >= 3 && parts[2].parse::<f64>().is_err() {
                header_skipped = true;
                info!("Skipped header line: {}", l);
                continue;
            }
            header_skipped = true;
        }

        let parts: Vec<&str> = l.split('\t').collect();
        if parts.len() >= 3 {
            match parts[2].parse::<f64>() {
                Ok(count) => {
                    if count >= 0.0 && count.is_finite() {
                        data.push(GeneCount {
                            gene_id: parts[0].trim().to_string(),
                            sample_id: parts[1].trim().to_string(),
                            count,
                        });
                    } else {
                        warn!("Invalid count value at line {}: {}", line_number, count);
                    }
                }
                Err(_) => {
                    warn!("Could not parse count at line {}: '{}'", line_number, parts[2]);
                }
            }
        } else {
            warn!("Insufficient columns at line {}: expected 3, found {}", line_number, parts.len());
        }
    }

    Ok(data)
}

fn create_count_matrix(data: &[GeneCount]) -> Result<(Array2<f64>, Vec<String>, Vec<String>), Box<dyn Error>> {
    // Collect unique genes and samples using HashSet for efficiency
    let mut gene_set = HashSet::new();
    let mut sample_set = HashSet::new();
    
    for entry in data {
        gene_set.insert(entry.gene_id.clone());
        sample_set.insert(entry.sample_id.clone());
    }

    let mut genes: Vec<String> = gene_set.into_iter().collect();
    let mut samples: Vec<String> = sample_set.into_iter().collect();
    
    genes.sort();
    samples.sort();

    // Create lookup maps for efficient indexing
    let gene_to_idx: HashMap<String, usize> = genes.iter().enumerate()
        .map(|(i, g)| (g.clone(), i)).collect();
    let sample_to_idx: HashMap<String, usize> = samples.iter().enumerate()
        .map(|(i, s)| (s.clone(), i)).collect();

    // Initialize matrix
    let mut matrix = Array2::<f64>::zeros((genes.len(), samples.len()));
    
    // Fill matrix
    for entry in data {
        if let (Some(&g_idx), Some(&s_idx)) = (
            gene_to_idx.get(&entry.gene_id),
            sample_to_idx.get(&entry.sample_id)
        ) {
            matrix[[g_idx, s_idx]] = entry.count;
        }
    }

    Ok((matrix, genes, samples))
}

fn perform_quality_checks(
    matrix: &Array2<f64>,
    genes: &[String],
    samples: &[String],
    min_count: f64,
) -> Result<(), Box<dyn Error>> {
    info!("Performing quality control checks");
    
    // Check for empty samples
    for (s_idx, sample) in samples.iter().enumerate() {
        let col_sum: f64 = matrix.column(s_idx).sum();
        if col_sum < min_count {
            warn!("Sample '{}' has very low total counts: {}", sample, col_sum);
        }
    }
    
    // Check for genes with all zero counts
    let mut zero_genes = 0;
    for (g_idx, _gene) in genes.iter().enumerate() {
        let row_sum: f64 = matrix.row(g_idx).sum();
        if row_sum == 0.0 {
            zero_genes += 1;
        }
    }
    
    if zero_genes > 0 {
        warn!("{} genes have zero counts across all samples", zero_genes);
    }
    
    info!("Quality control checks completed");
    Ok(())
}

fn calculate_size_factors(
    matrix: &Array2<f64>,
    genes: &[String],
    samples: &[String],
    pseudocount: f64,
) -> Result<Vec<f64>, Box<dyn Error>> {
    info!("Calculating geometric means");
    
    // Calculate geometric means for each gene
    let mut geom_means = Vec::with_capacity(genes.len());
    let mut valid_genes = 0;
    
    for g_idx in 0..genes.len() {
        let row = matrix.slice(s![g_idx, ..]);
        let mut log_sum = 0.0;
        let mut valid_samples = 0;
        
        for &count in row.iter() {
            if count > 0.0 {
                log_sum += (count + pseudocount).ln();
                valid_samples += 1;
            }
        }
        
        let gm = if valid_samples > 0 {
            (log_sum / valid_samples as f64).exp() - pseudocount
        } else {
            0.0
        };
        
        geom_means.push(gm.max(0.0));
        if gm > 0.0 {
            valid_genes += 1;
        }
    }
    
    info!("Calculated geometric means for {}/{} genes", valid_genes, genes.len());
    
    // Calculate size factors
    info!("Calculating size factors");
    let mut size_factors = Vec::with_capacity(samples.len());
    
    for s_idx in 0..samples.len() {
        let mut ratios = Vec::new();
        
        for g_idx in 0..genes.len() {
            let count_val = matrix[[g_idx, s_idx]];
            if geom_means[g_idx] > 0.0 && count_val > 0.0 {
                ratios.push(count_val / geom_means[g_idx]);
            }
        }
        
        if ratios.is_empty() {
            warn!("No valid ratios for sample '{}', using size factor of 1.0", samples[s_idx]);
            size_factors.push(1.0);
        } else {
            ratios.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            let median = if ratios.len() % 2 == 0 {
                (ratios[ratios.len() / 2 - 1] + ratios[ratios.len() / 2]) / 2.0
            } else {
                ratios[ratios.len() / 2]
            };
            size_factors.push(median.max(0.001)); // Avoid division by very small numbers
        }
    }
    
    Ok(size_factors)
}

fn normalize_counts(
    matrix: &Array2<f64>,
    size_factors: &[f64],
) -> Result<Array2<f64>, Box<dyn Error>> {
    let mut normalized = matrix.clone();
    
    for s_idx in 0..size_factors.len() {
        if size_factors[s_idx] > 0.0 {
            for g_idx in 0..matrix.nrows() {
                normalized[[g_idx, s_idx]] /= size_factors[s_idx];
            }
        }
    }
    
    Ok(normalized)
}

fn write_normalized_counts(
    matrix: &Array2<f64>,
    genes: &[String],
    samples: &[String],
    filename: &str,
) -> Result<(), Box<dyn Error>> {
    let mut file = File::create(filename)?;
    
    // Write header
    writeln!(file, "gene_id\tsample_id\tnormalized_count")?;
    
    // Write data
    for g_idx in 0..genes.len() {
        for s_idx in 0..samples.len() {
            writeln!(
                file,
                "{}\t{}\t{:.6}",
                genes[g_idx],
                samples[s_idx],
                matrix[[g_idx, s_idx]]
            )?;
        }
    }
    
    Ok(())
}

fn write_statistics(
    stats: &NormalizationStats,
    _genes: &[String],
    samples: &[String],
    filename: &str,
) -> Result<(), Box<dyn Error>> {
    let mut file = File::create(filename)?;
    
    writeln!(file, "RNA-seq Normalization Statistics")?;
    writeln!(file, "================================")?;
    writeln!(file, "Total genes: {}", stats.total_genes)?;
    writeln!(file, "Total samples: {}", stats.total_samples)?;
    writeln!(file, "Zero counts: {}", stats.zero_counts)?;
    writeln!(file, "Geometric means computed: {}", stats.geometric_means_computed)?;
    writeln!(file)?;
    
    writeln!(file, "Size Factors:")?;
    for (i, &sf) in stats.size_factors.iter().enumerate() {
        if i < samples.len() {
            writeln!(file, "{}\t{:.6}", samples[i], sf)?;
        }
    }
    
    Ok(())
}

fn count_zeros(matrix: &Array2<f64>) -> usize {
    matrix.iter().filter(|&&x| x == 0.0).count()
}