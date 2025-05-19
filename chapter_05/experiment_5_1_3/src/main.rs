use rayon::prelude::*;
use serde::{Serialize, Deserialize};
use anyhow::{Result, Context};
use std::collections::HashMap;
use std::fs;
use std::io::BufReader;
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

/// Command-line arguments for merging multiple JSON files.
#[derive(Parser, Debug)]
#[command(name = "merge_partial_outputs")]
struct Args {
    /// Paths to one or more partial JSON files to merge.
    #[arg(long, num_args = 1..)]
    inputs: Vec<String>,

    /// Name of the merged output JSON file.
    #[arg(long, default_value = "merged_output.json")]
    output: String,
}

fn main() -> Result<()> {
    let args = Args::parse();

    // Load all partial JSON files in parallel.
    let partials: Vec<PartialOutput> = args.inputs.par_iter()
        .map(|path| {
            let file = fs::File::open(path)
                .with_context(|| format!("Failed to open partial file {}", path))?;
            let reader = BufReader::new(file);
            let partial: PartialOutput = serde_json::from_reader(reader)
                .with_context(|| format!("Failed to parse JSON from file {}", path))?;
            Ok::<_, anyhow::Error>(partial)
        })
        .collect::<Result<_>>()?;

    // Merge histograms, total reads, and total bases.
    let mut merged_histogram = HashMap::new();
    let mut total_reads = 0usize;
    let mut total_bases = 0usize;

    for partial in partials {
        total_reads += partial.total_reads;
        total_bases += partial.total_bases;
        for (length, count) in partial.histogram.counts {
            *merged_histogram.entry(length).or_insert(0) += count;
        }
    }

    // Create a single PartialOutput for the merged data.
    let merged_output = PartialOutput {
        histogram: ReadLengthHistogram { counts: merged_histogram },
        total_reads,
        total_bases,
    };

    // Write the merged results to a JSON file.
    serde_json::to_writer(
        fs::File::create(&args.output)
            .with_context(|| format!("Failed to create merged output file {}", args.output))?,
        &merged_output,
    )?;

    println!("Merged {} partial files.", args.inputs.len());
    println!("Total reads: {}, total bases: {}", total_reads, total_bases);
    println!("Wrote merged output to {}", args.output);

    Ok(())
}
