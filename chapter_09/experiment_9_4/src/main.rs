use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use clap::{Arg, Command};
use serde::{Serialize, Deserialize};
use statrs::distribution::{StudentsT, ContinuousCDF};
use log::{info, warn, error};
use env_logger;

#[derive(Serialize, Deserialize, Debug, Clone)]
struct GeneCount {
    gene_id: String,
    sample_id: String,
    group: String,
    count: f64,
}

#[derive(Debug, Clone)]
struct DifferentialResult {
    gene_id: String,
    control_mean: f64,
    treatment_mean: f64,
    log2_fold_change: f64,
    p_value: f64,
    adjusted_p_value: f64,
    significant: bool,
}

#[derive(Debug)]
struct AnalysisStats {
    total_genes: usize,
    total_samples: usize,
    control_samples: usize,
    treatment_samples: usize,
    significant_genes: usize,
    upregulated_genes: usize,
    downregulated_genes: usize,
}

fn main() -> Result<(), Box<dyn Error>> {
    env_logger::init();
    
    let matches = Command::new("Differential Expression Analyzer")
        .version("1.0")
        .author("Bioinformatics Pipeline")
        .about("Performs differential expression analysis on normalized RNA-seq count data")
        .arg(Arg::new("input")
            .short('i')
            .long("input")
            .value_name("FILE")
            .help("Input TSV file with normalized counts")
            .default_value("normalized_counts.tsv"))
        .arg(Arg::new("metadata")
            .short('m')
            .long("metadata")
            .value_name("FILE")
            .help("Sample metadata file with group assignments")
            .default_value("sample_metadata.tsv"))
        .arg(Arg::new("output")
            .short('o')
            .long("output")
            .value_name("FILE")
            .help("Output file for differential expression results")
            .default_value("differential_expression.tsv"))
        .arg(Arg::new("stats")
            .short('s')
            .long("stats")
            .value_name("FILE")
            .help("Output file for analysis statistics")
            .default_value("de_analysis_stats.txt"))
        .arg(Arg::new("alpha")
            .long("alpha")
            .value_name("NUMBER")
            .help("Significance threshold for adjusted p-values")
            .default_value("0.05"))
        .arg(Arg::new("min_count")
            .long("min-count")
            .value_name("NUMBER")
            .help("Minimum mean count threshold for analysis")
            .default_value("10.0"))
        .arg(Arg::new("control_group")
            .long("control")
            .value_name("STRING")
            .help("Name of control group")
            .default_value("Control"))
        .arg(Arg::new("treatment_group")
            .long("treatment")
            .value_name("STRING")
            .help("Name of treatment group")
            .default_value("Treatment"))
        .get_matches();

    let input_file = matches.get_one::<String>("input").unwrap();
    let metadata_file = matches.get_one::<String>("metadata").unwrap();
    let output_file = matches.get_one::<String>("output").unwrap();
    let stats_file = matches.get_one::<String>("stats").unwrap();
    let alpha: f64 = matches.get_one::<String>("alpha").unwrap().parse()?;
    let min_count: f64 = matches.get_one::<String>("min_count").unwrap().parse()?;
    let control_group = matches.get_one::<String>("control_group").unwrap();
    let treatment_group = matches.get_one::<String>("treatment_group").unwrap();

    info!("Starting differential expression analysis");
    info!("Input file: {}", input_file);
    info!("Metadata file: {}", metadata_file);
    info!("Control group: {}", control_group);
    info!("Treatment group: {}", treatment_group);
    info!("Significance threshold: {}", alpha);
    info!("Minimum count threshold: {}", min_count);

    // Read sample metadata
    let sample_groups = read_sample_metadata(metadata_file)?;
    info!("Read metadata for {} samples", sample_groups.len());

    // Read normalized count data
    let data = read_count_data(input_file, &sample_groups)?;
    info!("Read {} count entries from input file", data.len());

    if data.is_empty() {
        error!("No valid data found in input file");
        return Err("Empty dataset".into());
    }

    // Organize data by gene and group
    let (gene_data, genes) = organize_data_by_gene(&data, control_group, treatment_group)?;
    info!("Organized data for {} genes", genes.len());

    // Filter genes by minimum count threshold
    let filtered_genes = filter_genes_by_count(&gene_data, &genes, min_count)?;
    info!("Filtered to {} genes meeting minimum count threshold", filtered_genes.len());

    // Perform differential expression analysis
    let de_results = perform_differential_analysis(&gene_data, &filtered_genes)?;
    info!("Completed differential expression analysis for {} genes", de_results.len());

    // Apply multiple testing correction
    let corrected_results = apply_benjamini_hochberg_correction(de_results, alpha)?;
    info!("Applied Benjamini-Hochberg correction");

    // Generate analysis statistics
    let stats = generate_analysis_stats(&corrected_results, &sample_groups, control_group, treatment_group)?;
    info!("Generated analysis statistics");

    // Write results
    write_results(&corrected_results, output_file)?;
    info!("Results written to {}", output_file);

    // Write statistics
    write_analysis_stats(&stats, &corrected_results, stats_file)?;
    info!("Statistics written to {}", stats_file);

    info!("Differential expression analysis completed successfully");
    Ok(())
}

fn read_sample_metadata(filename: &str) -> Result<HashMap<String, String>, Box<dyn Error>> {
    if !Path::new(filename).exists() {
        return Err(format!("Metadata file '{}' does not exist", filename).into());
    }

    let file = File::open(filename)?;
    let reader = BufReader::new(file);
    let mut sample_groups = HashMap::new();
    let mut header_skipped = false;

    for (line_number, line) in reader.lines().enumerate() {
        let l = line?;
        if l.trim().is_empty() {
            continue;
        }

        // Skip header line
        if !header_skipped {
            header_skipped = true;
            continue;
        }

        let parts: Vec<&str> = l.split('\t').collect();
        if parts.len() >= 2 {
            let sample_id = parts[0].trim().to_string();
            let group = parts[1].trim().to_string();
            sample_groups.insert(sample_id, group);
        } else {
            warn!("Invalid metadata format at line {}: expected 2+ columns, found {}", line_number + 1, parts.len());
        }
    }

    Ok(sample_groups)
}

fn read_count_data(filename: &str, sample_groups: &HashMap<String, String>) -> Result<Vec<GeneCount>, Box<dyn Error>> {
    if !Path::new(filename).exists() {
        return Err(format!("Input file '{}' does not exist", filename).into());
    }

    let file = File::open(filename)?;
    let reader = BufReader::new(file);
    let mut data = Vec::new();
    let mut header_skipped = false;

    for (line_number, line) in reader.lines().enumerate() {
        let l = line?;
        if l.trim().is_empty() {
            continue;
        }

        // Skip header line
        if !header_skipped {
            header_skipped = true;
            continue;
        }

        let parts: Vec<&str> = l.split('\t').collect();
        if parts.len() >= 3 {
            let gene_id = parts[0].trim().to_string();
            let sample_id = parts[1].trim().to_string();
            
            if let Some(group) = sample_groups.get(&sample_id) {
                match parts[2].parse::<f64>() {
                    Ok(count) => {
                        if count >= 0.0 && count.is_finite() {
                            data.push(GeneCount {
                                gene_id,
                                sample_id,
                                group: group.clone(),
                                count,
                            });
                        } else {
                            warn!("Invalid count value at line {}: {}", line_number + 1, count);
                        }
                    }
                    Err(_) => {
                        warn!("Could not parse count at line {}: '{}'", line_number + 1, parts[2]);
                    }
                }
            } else {
                warn!("Sample '{}' not found in metadata", sample_id);
            }
        } else {
            warn!("Insufficient columns at line {}: expected 3, found {}", line_number + 1, parts.len());
        }
    }

    Ok(data)
}

fn organize_data_by_gene(
    data: &[GeneCount], 
    control_group: &str, 
    treatment_group: &str
) -> Result<(HashMap<String, (Vec<f64>, Vec<f64>)>, Vec<String>), Box<dyn Error>> {
    let mut gene_data: HashMap<String, (Vec<f64>, Vec<f64>)> = HashMap::new();
    let mut genes = HashSet::new();

    for entry in data {
        genes.insert(entry.gene_id.clone());
        
        let (control_counts, treatment_counts) = gene_data.entry(entry.gene_id.clone())
            .or_insert_with(|| (Vec::new(), Vec::new()));

        if entry.group == control_group {
            control_counts.push(entry.count);
        } else if entry.group == treatment_group {
            treatment_counts.push(entry.count);
        }
    }

    let mut genes_vec: Vec<String> = genes.into_iter().collect();
    genes_vec.sort();

    Ok((gene_data, genes_vec))
}

fn filter_genes_by_count(
    gene_data: &HashMap<String, (Vec<f64>, Vec<f64>)>,
    genes: &[String],
    min_count: f64,
) -> Result<Vec<String>, Box<dyn Error>> {
    let mut filtered_genes = Vec::new();

    for gene in genes {
        if let Some((control_counts, treatment_counts)) = gene_data.get(gene) {
            let control_mean = if control_counts.is_empty() { 0.0 } else { 
                control_counts.iter().sum::<f64>() / control_counts.len() as f64 
            };
            let treatment_mean = if treatment_counts.is_empty() { 0.0 } else { 
                treatment_counts.iter().sum::<f64>() / treatment_counts.len() as f64 
            };

            if control_mean >= min_count || treatment_mean >= min_count {
                filtered_genes.push(gene.clone());
            }
        }
    }

    Ok(filtered_genes)
}

fn perform_differential_analysis(
    gene_data: &HashMap<String, (Vec<f64>, Vec<f64>)>,
    genes: &[String],
) -> Result<Vec<DifferentialResult>, Box<dyn Error>> {
    let mut results = Vec::new();

    for gene in genes {
        if let Some((control_counts, treatment_counts)) = gene_data.get(gene) {
            if control_counts.len() < 2 || treatment_counts.len() < 2 {
                warn!("Gene {} has insufficient replicates for statistical testing", gene);
                continue;
            }

            let control_mean = control_counts.iter().sum::<f64>() / control_counts.len() as f64;
            let treatment_mean = treatment_counts.iter().sum::<f64>() / treatment_counts.len() as f64;

            // Calculate log2 fold change (add pseudocount to avoid log(0))
            let log2_fold_change = if control_mean > 0.0 {
                ((treatment_mean + 1.0) / (control_mean + 1.0)).log2()
            } else {
                0.0
            };

            // Perform Welch's t-test
            let p_value = if control_counts.len() >= 2 && treatment_counts.len() >= 2 {
                welch_t_test(control_counts, treatment_counts)?
            } else {
                1.0 // No significance if insufficient samples
            };

            results.push(DifferentialResult {
                gene_id: gene.clone(),
                control_mean,
                treatment_mean,
                log2_fold_change,
                p_value,
                adjusted_p_value: p_value, // Will be corrected later
                significant: false, // Will be determined after correction
            });
        }
    }

    Ok(results)
}

fn welch_t_test(group1: &[f64], group2: &[f64]) -> Result<f64, Box<dyn Error>> {
    if group1.len() < 2 || group2.len() < 2 {
        return Ok(1.0);
    }

    // Calculate means
    let mean1 = group1.iter().sum::<f64>() / group1.len() as f64;
    let mean2 = group2.iter().sum::<f64>() / group2.len() as f64;
    
    // Calculate variances
    let var1 = if group1.len() > 1 {
        let sum_sq_diff: f64 = group1.iter().map(|x| (x - mean1).powi(2)).sum();
        sum_sq_diff / (group1.len() - 1) as f64
    } else {
        0.0
    };
    
    let var2 = if group2.len() > 1 {
        let sum_sq_diff: f64 = group2.iter().map(|x| (x - mean2).powi(2)).sum();
        sum_sq_diff / (group2.len() - 1) as f64
    } else {
        0.0
    };

    let n1 = group1.len() as f64;
    let n2 = group2.len() as f64;

    // Welch's t-test calculation
    let s1_sq_n1 = var1 / n1;
    let s2_sq_n2 = var2 / n2;
    let se = (s1_sq_n1 + s2_sq_n2).sqrt();

    if se == 0.0 {
        return Ok(1.0);
    }

    let t_stat = (mean1 - mean2) / se;

    // Welch-Satterthwaite equation for degrees of freedom
    let df_num = (s1_sq_n1 + s2_sq_n2).powi(2);
    let df_denom = (s1_sq_n1.powi(2) / (n1 - 1.0)) + (s2_sq_n2.powi(2) / (n2 - 1.0));
    let df = if df_denom > 0.0 { df_num / df_denom } else { 1.0 };

    // Calculate p-value using t-distribution
    let t_dist = StudentsT::new(0.0, 1.0, df).map_err(|e| format!("Error creating t-distribution: {}", e))?;
    let p_value = 2.0 * (1.0 - t_dist.cdf(t_stat.abs()));

    Ok(p_value.max(1e-10).min(1.0)) // Bound p-value
}

fn apply_benjamini_hochberg_correction(
    mut results: Vec<DifferentialResult>,
    alpha: f64,
) -> Result<Vec<DifferentialResult>, Box<dyn Error>> {
    // Sort by p-value
    results.sort_by(|a, b| a.p_value.partial_cmp(&b.p_value).unwrap_or(std::cmp::Ordering::Equal));

    let m = results.len() as f64;
    let mut corrected_results = Vec::new();

    // Apply Benjamini-Hochberg correction
    for (i, mut result) in results.into_iter().enumerate() {
        let rank = (i + 1) as f64;
        result.adjusted_p_value = (result.p_value * m / rank).min(1.0);
        result.significant = result.adjusted_p_value < alpha;
        corrected_results.push(result);
    }

    // Ensure monotonicity of adjusted p-values
    for i in (1..corrected_results.len()).rev() {
        if corrected_results[i].adjusted_p_value < corrected_results[i-1].adjusted_p_value {
            corrected_results[i-1].adjusted_p_value = corrected_results[i].adjusted_p_value;
        }
    }

    // Re-sort by gene name for consistent output
    corrected_results.sort_by(|a, b| a.gene_id.cmp(&b.gene_id));

    Ok(corrected_results)
}

fn generate_analysis_stats(
    results: &[DifferentialResult],
    sample_groups: &HashMap<String, String>,
    control_group: &str,
    treatment_group: &str,
) -> Result<AnalysisStats, Box<dyn Error>> {
    let total_genes = results.len();
    let significant_genes = results.iter().filter(|r| r.significant).count();
    let upregulated_genes = results.iter().filter(|r| r.significant && r.log2_fold_change > 0.0).count();
    let downregulated_genes = results.iter().filter(|r| r.significant && r.log2_fold_change < 0.0).count();

    let control_samples = sample_groups.values().filter(|&g| g == control_group).count();
    let treatment_samples = sample_groups.values().filter(|&g| g == treatment_group).count();
    let total_samples = control_samples + treatment_samples;

    Ok(AnalysisStats {
        total_genes,
        total_samples,
        control_samples,
        treatment_samples,
        significant_genes,
        upregulated_genes,
        downregulated_genes,
    })
}

fn write_results(results: &[DifferentialResult], filename: &str) -> Result<(), Box<dyn Error>> {
    let mut file = File::create(filename)?;

    // Write header
    writeln!(file, "gene_id\tcontrol_mean\ttreatment_mean\tlog2_fold_change\tp_value\tadjusted_p_value\tsignificant")?;

    // Write results
    for result in results {
        writeln!(
            file,
            "{}\t{:.6}\t{:.6}\t{:.6}\t{:.2e}\t{:.2e}\t{}",
            result.gene_id,
            result.control_mean,
            result.treatment_mean,
            result.log2_fold_change,
            result.p_value,
            result.adjusted_p_value,
            result.significant
        )?;
    }

    Ok(())
}

fn write_analysis_stats(
    stats: &AnalysisStats,
    results: &[DifferentialResult],
    filename: &str,
) -> Result<(), Box<dyn Error>> {
    let mut file = File::create(filename)?;

    writeln!(file, "Differential Expression Analysis Statistics")?;
    writeln!(file, "==========================================")?;
    writeln!(file, "Total genes analyzed: {}", stats.total_genes)?;
    writeln!(file, "Total samples: {}", stats.total_samples)?;
    writeln!(file, "Control samples: {}", stats.control_samples)?;
    writeln!(file, "Treatment samples: {}", stats.treatment_samples)?;
    writeln!(file)?;
    writeln!(file, "Significant genes: {} ({:.1}%)", 
             stats.significant_genes, 
             stats.significant_genes as f64 / stats.total_genes as f64 * 100.0)?;
    writeln!(file, "Upregulated genes: {} ({:.1}%)", 
             stats.upregulated_genes,
             stats.upregulated_genes as f64 / stats.total_genes as f64 * 100.0)?;
    writeln!(file, "Downregulated genes: {} ({:.1}%)", 
             stats.downregulated_genes,
             stats.downregulated_genes as f64 / stats.total_genes as f64 * 100.0)?;
    writeln!(file)?;

    // Write top significant genes
    let mut sig_results: Vec<_> = results.iter().filter(|r| r.significant).collect();
    sig_results.sort_by(|a, b| a.adjusted_p_value.partial_cmp(&b.adjusted_p_value).unwrap());

    writeln!(file, "Top 10 Most Significant Genes:")?;
    writeln!(file, "gene_id\tlog2FC\tadj_p_value")?;
    for result in sig_results.iter().take(10) {
        writeln!(file, "{}\t{:.3}\t{:.2e}", 
                 result.gene_id, 
                 result.log2_fold_change, 
                 result.adjusted_p_value)?;
    }

    Ok(())
}