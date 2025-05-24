use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::time::Instant;

use anyhow::{Context, Result};
use clap::Parser;
use dashmap::DashMap;
use flate2::read::GzDecoder;
use indicatif::{ProgressBar, ProgressStyle};
use log::{info, warn, error};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use thiserror::Error;

#[derive(Error, Debug)]
pub enum PseudoAlignError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
    #[error("JSON parsing error: {0}")]
    Json(#[from] serde_json::Error),
    #[error("Invalid k-mer length: expected {expected}, got {actual}")]
    InvalidKmerLength { expected: usize, actual: usize },
    #[error("Empty index file")]
    EmptyIndex,
    #[error("No reads found in input file")]
    NoReads,
}

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// K-mer index file (JSON format)
    #[arg(short, long, default_value = "data/kmer_index.json")]
    index: String,

    /// Input reads file (FASTQ format, optionally gzipped)
    #[arg(short, long, default_value = "data/reads.fastq")]
    reads: String,

    /// Output quantification file
    #[arg(short, long, default_value = "results/quantification.tsv")]
    output: String,

    /// K-mer length
    #[arg(short, long, default_value_t = 31)]
    kmer_length: usize,

    /// Number of threads to use
    #[arg(short, long, default_value_t = 4)]
    threads: usize,

    /// Minimum read length to process
    #[arg(long, default_value_t = 50)]
    min_read_length: usize,

    /// Verbose logging
    #[arg(short, long)]
    verbose: bool,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct KmerIndex {
    pub kmer: String,
    pub transcripts: Vec<String>,
    pub transcript_positions: Option<Vec<usize>>,
}

#[derive(Debug, Clone)]
pub struct TranscriptQuantification {
    pub transcript_id: String,
    pub count: f64,
    pub tpm: f64,
    pub effective_length: f64,
}

#[derive(Debug)]
pub struct QuantificationResults {
    pub transcripts: Vec<TranscriptQuantification>,
    pub total_reads: usize,
    pub aligned_reads: usize,
    pub processing_time: std::time::Duration,
}

pub struct PseudoAligner {
    kmer_index: HashMap<String, Vec<String>>,
    transcript_lengths: HashMap<String, usize>,
    kmer_length: usize,
    min_read_length: usize,
}

impl PseudoAligner {
    pub fn new(index_path: &str, kmer_length: usize, min_read_length: usize) -> Result<Self> {
        let start = Instant::now();
        info!("Loading k-mer index from: {}", index_path);
        
        let file = File::open(index_path)
            .with_context(|| format!("Failed to open index file: {}", index_path))?;
        let reader = BufReader::new(file);
        let kmer_entries: Vec<KmerIndex> = serde_json::from_reader(reader)
            .context("Failed to parse k-mer index JSON")?;

        if kmer_entries.is_empty() {
            return Err(PseudoAlignError::EmptyIndex.into());
        }

        let mut kmer_index = HashMap::new();
        let mut transcript_lengths = HashMap::new();

        for entry in kmer_entries {
            if entry.kmer.len() != kmer_length {
                warn!(
                    "K-mer length mismatch: expected {}, got {} for k-mer {}",
                    kmer_length,
                    entry.kmer.len(),
                    entry.kmer
                );
                continue;
            }
            kmer_index.insert(entry.kmer, entry.transcripts.clone());
            
            for transcript in &entry.transcripts {
                transcript_lengths.entry(transcript.clone()).or_insert(1000);
            }
        }

        let load_time = start.elapsed();
        info!(
            "Loaded {} k-mers for {} transcripts in {:?}",
            kmer_index.len(),
            transcript_lengths.len(),
            load_time
        );

        Ok(PseudoAligner {
            kmer_index,
            transcript_lengths,
            kmer_length,
            min_read_length,
        })
    }

    pub fn quantify_reads(&self, reads_path: &str, num_threads: usize) -> Result<QuantificationResults> {
        let start = Instant::now();
        rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build_global()
            .context("Failed to initialize thread pool")?;

        let transcript_counts: Arc<DashMap<String, AtomicUsize>> = Arc::new(DashMap::new());
        let total_reads = Arc::new(AtomicUsize::new(0));

        for transcript in self.transcript_lengths.keys() {
            transcript_counts.insert(transcript.clone(), AtomicUsize::new(0));
        }

        info!("Processing reads from: {}", reads_path);
        let pb = ProgressBar::new_spinner();
        pb.set_style(
            ProgressStyle::default_spinner()
                .template("{spinner:.green} [{elapsed_precise}] {msg}")
                .unwrap(),
        );

        let reader: Box<dyn BufRead> = if reads_path.ends_with(".gz") {
            Box::new(BufReader::new(GzDecoder::new(File::open(reads_path)?)))
        } else {
            Box::new(BufReader::new(File::open(reads_path)?))
        };

        let mut lines = reader.lines();
        let mut read_count = 0;
        let mut sequences = Vec::new();

        while let Some(header_line) = lines.next() {
            let _header = header_line?;
            
            if let Some(sequence_line) = lines.next() {
                let sequence = sequence_line?;
                total_reads.fetch_add(1, Ordering::Relaxed);
                read_count += 1;

                if read_count % 10000 == 0 {
                    pb.set_message(format!("Reading {} sequences", read_count));
                }

                if sequence.len() >= self.min_read_length {
                    sequences.push(sequence);
                }
            }

            if lines.next().is_some() {
                lines.next();
            }
        }

        pb.finish_with_message(format!("Read {} sequences", sequences.len()));

        let pb2 = ProgressBar::new(sequences.len() as u64);
        pb2.set_style(
            ProgressStyle::default_bar()
                .template("{bar:40.cyan/blue} {pos:>7}/{len:7} {msg}")
                .unwrap(),
        );

        let kmer_index = &self.kmer_index;
        let kmer_length = self.kmer_length;
        let aligned_count = AtomicUsize::new(0);

        sequences.par_iter().for_each(|sequence| {
            let mut transcript_hits: HashMap<String, usize> = HashMap::new();
            
            for i in 0..=sequence.len().saturating_sub(kmer_length) {
                let kmer = &sequence[i..i + kmer_length];
                if let Some(transcripts) = kmer_index.get(kmer) {
                    for transcript in transcripts {
                        *transcript_hits.entry(transcript.clone()).or_insert(0) += 1;
                    }
                }
            }

            if !transcript_hits.is_empty() {
                aligned_count.fetch_add(1, Ordering::Relaxed);
                
                if let Some((best_transcript, _)) = transcript_hits
                    .iter()
                    .max_by_key(|(_, &count)| count)
                {
                    if let Some(counter) = transcript_counts.get(best_transcript) {
                        counter.fetch_add(1, Ordering::Relaxed);
                    }
                }
            }

            pb2.inc(1);
        });

        pb2.finish_with_message("Processing complete");

        let mut transcripts = Vec::new();
        let total_counts: f64 = transcript_counts
            .iter()
            .map(|entry| entry.value().load(Ordering::Relaxed) as f64)
            .sum();

        for entry in transcript_counts.iter() {
            let transcript_id = entry.key().clone();
            let count = entry.value().load(Ordering::Relaxed) as f64;
            let effective_length = *self.transcript_lengths.get(&transcript_id).unwrap_or(&1000) as f64;
            
            let tpm = if total_counts > 0.0 {
                (count / effective_length) * 1_000_000.0 / 
                (transcript_counts
                    .iter()
                    .map(|e| e.value().load(Ordering::Relaxed) as f64 / 
                         *self.transcript_lengths.get(e.key()).unwrap_or(&1000) as f64)
                    .sum::<f64>())
            } else {
                0.0
            };

            transcripts.push(TranscriptQuantification {
                transcript_id,
                count,
                tpm,
                effective_length,
            });
        }

        transcripts.sort_by(|a, b| b.count.partial_cmp(&a.count).unwrap());

        Ok(QuantificationResults {
            transcripts,
            total_reads: total_reads.load(Ordering::Relaxed),
            aligned_reads: aligned_count.load(Ordering::Relaxed),
            processing_time: start.elapsed(),
        })
    }

    pub fn write_results(&self, results: &QuantificationResults, output_path: &str) -> Result<()> {
        if let Some(parent) = Path::new(output_path).parent() {
            std::fs::create_dir_all(parent)?;
        }

        let file = File::create(output_path)
            .with_context(|| format!("Failed to create output file: {}", output_path))?;
        let mut writer = BufWriter::new(file);

        writeln!(writer, "transcript_id\tcount\ttpm\teffective_length")?;

        for quant in &results.transcripts {
            writeln!(
                writer,
                "{}\t{:.2}\t{:.6}\t{:.0}",
                quant.transcript_id, quant.count, quant.tpm, quant.effective_length
            )?;
        }

        writer.flush()?;
        info!("Results written to: {}", output_path);

        println!("\n=== Quantification Summary ===");
        println!("Total reads processed: {}", results.total_reads);
        println!("Reads aligned: {}", results.aligned_reads);
        println!(
            "Alignment rate: {:.2}%",
            (results.aligned_reads as f64 / results.total_reads as f64) * 100.0
        );
        println!("Processing time: {:?}", results.processing_time);
        println!(
            "Transcripts with non-zero counts: {}",
            results.transcripts.iter().filter(|t| t.count > 0.0).count()
        );

        Ok(())
    }
}

fn main() -> Result<()> {
    let args = Args::parse();

    let log_level = if args.verbose { "debug" } else { "info" };
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(log_level)).init();

    info!("Starting pseudo-alignment with parameters:");
    info!("  Index file: {}", args.index);
    info!("  Reads file: {}", args.reads);
    info!("  Output file: {}", args.output);
    info!("  K-mer length: {}", args.kmer_length);
    info!("  Threads: {}", args.threads);
    info!("  Min read length: {}", args.min_read_length);

    if !Path::new(&args.index).exists() {
        error!("Index file does not exist: {}", args.index);
        std::process::exit(1);
    }
    if !Path::new(&args.reads).exists() {
        error!("Reads file does not exist: {}", args.reads);
        std::process::exit(1);
    }

    let aligner = PseudoAligner::new(&args.index, args.kmer_length, args.min_read_length)
        .context("Failed to initialize pseudo-aligner")?;

    let results = aligner
        .quantify_reads(&args.reads, args.threads)
        .context("Failed to quantify reads")?;

    aligner
        .write_results(&results, &args.output)
        .context("Failed to write results")?;

    info!("Pseudo-alignment completed successfully!");
    Ok(())
}
