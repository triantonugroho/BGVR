use rayon::prelude::*;
use needletail::{parse_fastx_file};
use serde::{Serialize, Deserialize};
use anyhow::{Result, Context};
use std::collections::HashMap;
use std::fs::File;
use clap::Parser;

/// A histogram of read lengths.
#[derive(Debug, Serialize, Deserialize)]
struct ReadLengthHistogram {
    counts: HashMap<usize, usize>,
}

/// Struct for storing partial results (histogram, total reads, total bases).
#[derive(Debug, Serialize, Deserialize)]
struct PartialOutput {
    histogram: ReadLengthHistogram,
    total_reads: usize,
    total_bases: usize,
}

/// Command-line arguments for specifying input/output files.
#[derive(Parser, Debug)]
#[command(name = "hpc_fastq_stats")]
struct Args {
    /// Path to the input FASTQ file (supports gzipped FASTQ).
    #[arg(long)]
    input: String,

    /// Name of the output JSON file for partial results.
    #[arg(long, default_value = "partial_output.json")]
    output: String,
}

fn main() -> Result<()> {
    // Parse command-line arguments.
    let args = Args::parse();

    // Open the FASTQ file with needletail, providing an error context if it fails.
    let mut reader = parse_fastx_file(&args.input)
        .with_context(|| format!("Failed to open FASTQ file at {}", args.input))?;

    // Initialize stats containers.
    let mut histogram = HashMap::new();
    let mut total_reads = 0usize;
    let mut total_bases = 0usize;

    // Read each record in a streaming fashion.
    while let Some(record) = reader.next() {
        let seqrec = record?;
        let len = seqrec.seq().len();
        *histogram.entry(len).or_insert(0) += 1;
        total_reads += 1;
        total_bases += len;
    }

    // Package results into a PartialOutput struct.
    let partial = PartialOutput {
        histogram: ReadLengthHistogram { counts: histogram },
        total_reads,
        total_bases,
    };

    // Write partial results to a JSON file.
    serde_json::to_writer(
        File::create(&args.output)
            .with_context(|| format!("Failed to create output file at {}", args.output))?,
        &partial,
    )?;

    println!("Wrote partial output to {}", args.output);
    println!("Reads: {}, Bases: {}", total_reads, total_bases);
    Ok(())
}
