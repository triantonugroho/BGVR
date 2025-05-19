use anyhow::{Context, Result};
use clap::{Parser, ValueEnum};
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

// Custom error types
#[derive(thiserror::Error, Debug)]
pub enum RnaSeqError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
    #[error("Parse error: {0}")]
    Parse(#[from] std::num::ParseFloatError),
    #[error("Math error: {0}")]
    Math(String),
    #[error("Validation error: {0}")]
    Validation(String),
}

// Normalization methods
#[derive(Debug, Clone, ValueEnum, Serialize, Deserialize)]
pub enum NormalizationMethod {
    /// Trimmed Mean of M-values normalization
    Tmm,
    /// Counts Per Million normalization
    Cpm,
    /// Standard mean normalization
    Standard,
    /// DESeq2-style median-of-ratios normalization
    DeseqMor,
    /// Upper quartile normalization
    UpperQuartile,
}

impl std::fmt::Display for NormalizationMethod {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            NormalizationMethod::Tmm => write!(f, "TMM"),
            NormalizationMethod::Cpm => write!(f, "CPM"),
            NormalizationMethod::Standard => write!(f, "Standard"),
            NormalizationMethod::DeseqMor => write!(f, "DESeq-MOR"),
            NormalizationMethod::UpperQuartile => write!(f, "UpperQuartile"),
        }
    }
}

// Gene count structure
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeneCount {
    pub gene_id: String,
    pub raw_count: f64,
    pub normalized_count: Option<f64>,
}

impl GeneCount {
    pub fn new(gene_id: String, raw_count: f64) -> Result<Self> {
        if gene_id.trim().is_empty() {
            return Err(RnaSeqError::Validation("Gene ID cannot be empty".to_string()).into());
        }
        if raw_count < 0.0 || !raw_count.is_finite() {
            return Err(RnaSeqError::Validation(format!("Invalid count: {}", raw_count)).into());
        }
        
        Ok(Self {
            gene_id,
            raw_count,
            normalized_count: None,
        })
    }
}

// Sample information
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Sample {
    pub name: String,
    pub genes: Vec<GeneCount>,
}

// Main analyzer class
pub struct RnaSeqAnalyzer {
    pub samples: Vec<Sample>,
    pub normalization_method: NormalizationMethod,
    pub min_count: f64,
    pub min_samples: usize,
    pub verbose: bool,
}

impl RnaSeqAnalyzer {
    pub fn new(
        normalization_method: NormalizationMethod,
        min_count: f64,
        min_samples: usize,
        verbose: bool,
    ) -> Self {
        Self {
            samples: Vec::new(),
            normalization_method,
            min_count,
            min_samples,
            verbose,
        }
    }

    /// Load samples from count files
    pub fn load_samples_from_files(&mut self, count_files: Vec<PathBuf>) -> Result<()> {
        let pb = if self.verbose {
            let pb = ProgressBar::new(count_files.len() as u64);
            pb.set_style(
                ProgressStyle::default_bar()
                    .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({msg})")?
                    .progress_chars("#>-"),
            );
            Some(pb)
        } else {
            None
        };

        let samples: Result<Vec<_>> = count_files
            .par_iter()
            .map(|file_path| {
                let sample_name = file_path
                    .file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or("unknown")
                    .to_string();
                
                let genes = self.load_gene_counts(file_path)?;
                
                if let Some(pb) = &pb {
                    pb.set_message(sample_name.clone());
                    pb.inc(1);
                }
                
                Ok(Sample {
                    name: sample_name,
                    genes,
                })
            })
            .collect();

        if let Some(pb) = pb {
            pb.finish_with_message("Loaded all samples");
        }

        self.samples = samples?;
        Ok(())
    }

    /// Load gene counts from a single file
    fn load_gene_counts(&self, file_path: &Path) -> Result<Vec<GeneCount>> {
        let file = File::open(file_path)
            .with_context(|| format!("Failed to open file: {}", file_path.display()))?;
        
        let reader = BufReader::new(file);
        let mut genes = Vec::new();

        for (line_num, line) in reader.lines().enumerate() {
            let line = line?;
            
            // Skip empty lines, comments, and HTSeq statistics
            if line.trim().is_empty() || line.starts_with('#') || line.starts_with("__") {
                continue;
            }
            
            let parts: Vec<&str> = line.trim().split('\t').collect();
            if parts.len() < 2 {
                if self.verbose {
                    eprintln!("Warning: Skipping malformed line {} in {}: {}", 
                             line_num + 1, file_path.display(), line);
                }
                continue;
            }
            
            match parts[1].parse::<f64>() {
                Ok(count) => {
                    if let Ok(gene) = GeneCount::new(parts[0].to_string(), count) {
                        genes.push(gene);
                    }
                }
                Err(_) => {
                    if self.verbose {
                        eprintln!("Warning: Invalid count '{}' on line {} in {}", 
                                 parts[1], line_num + 1, file_path.display());
                    }
                }
            }
        }

        Ok(genes)
    }

    /// Filter genes based on minimum count and sample criteria
    pub fn filter_genes(&mut self) -> Result<()> {
        if self.samples.is_empty() {
            return Err(RnaSeqError::Validation("No samples loaded".to_string()).into());
        }

        // Get all unique gene IDs
        let mut all_gene_ids: Vec<String> = self.samples[0]
            .genes
            .iter()
            .map(|g| g.gene_id.clone())
            .collect();
        all_gene_ids.sort();
        all_gene_ids.dedup();

        // Filter genes that meet criteria
        let filtered_genes: Vec<String> = all_gene_ids
            .par_iter()
            .filter_map(|gene_id| {
                let mut samples_with_min_count = 0;
                
                for sample in &self.samples {
                    if let Some(gene) = sample.genes.iter().find(|g| g.gene_id == *gene_id) {
                        if gene.raw_count >= self.min_count {
                            samples_with_min_count += 1;
                        }
                    }
                }
                
                if samples_with_min_count >= self.min_samples {
                    Some(gene_id.clone())
                } else {
                    None
                }
            })
            .collect();

        // Filter samples to only include filtered genes
        for sample in &mut self.samples {
            sample.genes.retain(|gene| filtered_genes.contains(&gene.gene_id));
            sample.genes.sort_by(|a, b| a.gene_id.cmp(&b.gene_id));
        }

        if self.verbose {
            println!("Filtered {} genes (from {} total)", 
                    filtered_genes.len(), all_gene_ids.len());
        }

        Ok(())
    }

    /// Normalize counts using the specified method
    pub fn normalize(&mut self) -> Result<()> {
        if self.samples.is_empty() {
            return Err(RnaSeqError::Validation("No samples to normalize".to_string()).into());
        }

        match self.normalization_method {
            NormalizationMethod::Tmm => self.normalize_tmm(),
            NormalizationMethod::Cpm => self.normalize_cpm(),
            NormalizationMethod::Standard => self.normalize_standard(),
            NormalizationMethod::DeseqMor => self.normalize_deseq_mor(),
            NormalizationMethod::UpperQuartile => self.normalize_upper_quartile(),
        }
    }

    /// TMM normalization (simplified)
    fn normalize_tmm(&mut self) -> Result<()> {
        // Simplified TMM using library size normalization
        let library_sizes: Vec<f64> = self.samples
            .iter()
            .map(|sample| sample.genes.iter().map(|g| g.raw_count).sum())
            .collect();

        let reference_library_size = library_sizes.iter().cloned().fold(0.0, f64::max);
        
        for (i, sample) in self.samples.iter_mut().enumerate() {
            let normalization_factor = library_sizes[i] / reference_library_size;
            for gene in &mut sample.genes {
                gene.normalized_count = Some(gene.raw_count / normalization_factor);
            }
        }

        if self.verbose {
            println!("TMM normalization completed");
        }
        Ok(())
    }

    /// CPM normalization
    fn normalize_cpm(&mut self) -> Result<()> {
        for sample in &mut self.samples {
            let total_counts: f64 = sample.genes.iter().map(|g| g.raw_count).sum();
            
            if total_counts == 0.0 {
                return Err(RnaSeqError::Math("Total counts is zero".to_string()).into());
            }

            for gene in &mut sample.genes {
                gene.normalized_count = Some((gene.raw_count / total_counts) * 1_000_000.0);
            }
        }

        if self.verbose {
            println!("CPM normalization completed");
        }
        Ok(())
    }

    /// Standard mean normalization
    fn normalize_standard(&mut self) -> Result<()> {
        for sample in &mut self.samples {
            let raw_counts: Vec<f64> = sample.genes.iter().map(|g| g.raw_count).collect();
            let mean_count = raw_counts.iter().cloned().fold(0.0, |a, b| a + b) / raw_counts.len() as f64;
            
            if mean_count == 0.0 {
                return Err(RnaSeqError::Math("Mean count is zero".to_string()).into());
            }

            for gene in &mut sample.genes {
                gene.normalized_count = Some(gene.raw_count / mean_count);
            }
        }

        if self.verbose {
            println!("Standard normalization completed");
        }
        Ok(())
    }

    /// DESeq2-style median of ratios normalization
    fn normalize_deseq_mor(&mut self) -> Result<()> {
        // Calculate geometric mean for each gene across samples
        let mut gene_geometric_means = HashMap::new();
        
        if let Some(first_sample) = self.samples.first() {
            for gene in &first_sample.genes {
                let gene_id = &gene.gene_id;
                let mut counts = Vec::new();
                
                for sample in &self.samples {
                    if let Some(g) = sample.genes.iter().find(|g| g.gene_id == *gene_id) {
                        if g.raw_count > 0.0 {
                            counts.push(g.raw_count);
                        }
                    }
                }
                
                if !counts.is_empty() {
                    let log_sum: f64 = counts.iter().map(|&x| x.ln()).sum();
                    let geometric_mean = (log_sum / counts.len() as f64).exp();
                    gene_geometric_means.insert(gene_id.clone(), geometric_mean);
                }
            }
        }

        // Calculate normalization factors for each sample
        for sample in &mut self.samples {
            let mut ratios = Vec::new();
            
            for gene in &sample.genes {
                if let Some(&geometric_mean) = gene_geometric_means.get(&gene.gene_id) {
                    if geometric_mean > 0.0 && gene.raw_count > 0.0 {
                        ratios.push(gene.raw_count / geometric_mean);
                    }
                }
            }
            
            if !ratios.is_empty() {
                // Calculate median manually
                ratios.sort_by(|a, b| a.partial_cmp(b).unwrap());
                let size_factor = if ratios.len() % 2 == 0 {
                    let mid = ratios.len() / 2;
                    (ratios[mid - 1] + ratios[mid]) / 2.0
                } else {
                    ratios[ratios.len() / 2]
                };
                
                for gene in &mut sample.genes {
                    gene.normalized_count = Some(gene.raw_count / size_factor);
                }
            }
        }

        if self.verbose {
            println!("DESeq median-of-ratios normalization completed");
        }
        Ok(())
    }

    /// Upper quartile normalization
    fn normalize_upper_quartile(&mut self) -> Result<()> {
        for sample in &mut self.samples {
            let mut counts: Vec<f64> = sample.genes.iter()
                .map(|g| g.raw_count)
                .filter(|&x| x > 0.0)
                .collect();
            
            if counts.is_empty() {
                return Err(RnaSeqError::Math("No non-zero counts found".to_string()).into());
            }
            
            // Sort and calculate 75th percentile (upper quartile)
            counts.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let index = (counts.len() as f64 * 0.75) as usize;
            let upper_quartile = if index >= counts.len() {
                counts[counts.len() - 1]
            } else {
                counts[index]
            };
            
            if upper_quartile == 0.0 {
                return Err(RnaSeqError::Math("Upper quartile is zero".to_string()).into());
            }

            for gene in &mut sample.genes {
                gene.normalized_count = Some(gene.raw_count / upper_quartile * 1_000_000.0);
            }
        }

        if self.verbose {
            println!("Upper quartile normalization completed");
        }
        Ok(())
    }

    /// Write results to output files
    pub fn write_results(&self, output_prefix: &str) -> Result<()> {
        let counts_file = format!("{}_normalized_counts.tsv", output_prefix);
        let mut writer = BufWriter::new(File::create(&counts_file)?);
        
        // Write header
        write!(writer, "gene_id")?;
        for sample in &self.samples {
            write!(writer, "\t{}", sample.name)?;
        }
        writeln!(writer)?;

        // Write data
        if let Some(first_sample) = self.samples.first() {
            for (gene_idx, gene) in first_sample.genes.iter().enumerate() {
                write!(writer, "{}", gene.gene_id)?;
                for sample in &self.samples {
                    if let Some(normalized_count) = sample.genes[gene_idx].normalized_count {
                        write!(writer, "\t{:.6}", normalized_count)?;
                    } else {
                        write!(writer, "\tNA")?;
                    }
                }
                writeln!(writer)?;
            }
        }
        writer.flush()?;

        // Write summary statistics
        let summary_file = format!("{}_summary.json", output_prefix);
        let summary_stats = self.generate_summary_stats();
        let summary_json = serde_json::to_string_pretty(&summary_stats)?;
        std::fs::write(summary_file, summary_json)?;

        if self.verbose {
            println!("Results written to {} and {}_summary.json", counts_file, output_prefix);
        }

        Ok(())
    }

    /// Generate summary statistics
    fn generate_summary_stats(&self) -> serde_json::Value {
        let total_genes = if let Some(sample) = self.samples.first() {
            sample.genes.len()
        } else {
            0
        };

        let library_sizes: Vec<f64> = self.samples
            .iter()
            .map(|sample| sample.genes.iter().map(|g| g.raw_count).sum())
            .collect();

        let detected_genes_per_sample: Vec<usize> = self.samples
            .iter()
            .map(|sample| sample.genes.iter().filter(|g| g.raw_count > 0.0).count())
            .collect();

        serde_json::json!({
            "total_genes": total_genes,
            "total_samples": self.samples.len(),
            "library_sizes": library_sizes,
            "detected_genes_per_sample": detected_genes_per_sample,
            "normalization_method": self.normalization_method.to_string()
        })
    }
}

// CLI interface
#[derive(Parser)]
#[command(name = "rnaseq-analyzer")]
#[command(about = "A comprehensive RNA-seq count normalization and analysis tool")]
pub struct Cli {
    /// Input count files (tab-separated: gene_id count)
    #[arg(short, long, value_delimiter = ',')]
    pub input: Vec<PathBuf>,

    /// Output prefix for result files
    #[arg(short, long, default_value = "output")]
    pub output: String,

    /// Normalization method
    #[arg(short, long, default_value = "cpm")]
    pub method: NormalizationMethod,

    /// Minimum count threshold
    #[arg(long, default_value = "1.0")]
    pub min_count: f64,

    /// Minimum number of samples with min_count
    #[arg(long, default_value = "1")]
    pub min_samples: usize,

    /// Verbose output
    #[arg(short, long)]
    pub verbose: bool,

    /// Show summary statistics
    #[arg(short, long)]
    pub summary: bool,
}

fn main() -> Result<()> {
    env_logger::init();
    let cli = Cli::parse();

    if cli.input.is_empty() {
        eprintln!("Error: No input files specified");
        std::process::exit(1);
    }

    // Initialize analyzer
    let mut analyzer = RnaSeqAnalyzer::new(
        cli.method.clone(),
        cli.min_count,
        cli.min_samples,
        cli.verbose,
    );

    // Load samples
    if cli.verbose {
        println!("Loading {} count files...", cli.input.len());
    }
    analyzer.load_samples_from_files(cli.input)?;

    // Filter genes
    analyzer.filter_genes()?;

    // Normalize counts
    if cli.verbose {
        println!("Normalizing counts using {} method...", cli.method);
    }
    analyzer.normalize()?;

    // Write results
    analyzer.write_results(&cli.output)?;

    // Display summary if requested
    if cli.summary {
        let stats = analyzer.generate_summary_stats();
        println!("\n=== ANALYSIS SUMMARY ===");
        println!("Total genes: {}", stats["total_genes"]);
        println!("Total samples: {}", stats["total_samples"]);
        println!("Normalization method: {}", cli.method);
        
        println!("\nLibrary sizes:");
        if let Some(sizes) = stats["library_sizes"].as_array() {
            for (i, size) in sizes.iter().enumerate() {
                println!("  Sample {}: {:.0}", i + 1, size.as_f64().unwrap_or(0.0));
            }
        }
        
        println!("\nDetected genes per sample:");
        if let Some(detected) = stats["detected_genes_per_sample"].as_array() {
            for (i, count) in detected.iter().enumerate() {
                println!("  Sample {}: {}", i + 1, count.as_u64().unwrap_or(0));
            }
        }
    }

    if cli.verbose {
        println!("\nAnalysis completed successfully!");
    }

    Ok(())
}