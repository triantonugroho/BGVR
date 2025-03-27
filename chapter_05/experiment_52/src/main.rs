use rayon::prelude::*;
use needletail::{parse_fastx_file, FastxReader, Sequence};
use anyhow::{Result, Context};
use serde::{Serialize, Deserialize};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use clap::Parser;

/// Command-line arguments for chunk-based trimming.
#[derive(Parser, Debug)]
#[command(name = "sliding_window_trimmer")]
#[command(about = "Sliding window trimmer with chunked parallel processing")]
struct Args {
    /// Input FASTQ file (gzipped or plain).
    #[arg(long)]
    input: PathBuf,

    /// Prefix for chunk-based partial output. For each chunk, a file will be created with this prefix.
    #[arg(long, default_value = "partial_trim_output")]
    output_prefix: String,

    /// Maximum number of reads to process in one chunk.
    #[arg(long, default_value_t = 10_000)]
    chunk_size: usize,

    /// Sliding window size for trimming.
    #[arg(long, default_value_t = 4)]
    window_size: usize,

    /// Average phred-score threshold within the window.
    #[arg(long, default_value_t = 20)]
    quality_threshold: u8,
}

/// Stores the trimmed result for each read.
#[derive(Serialize, Deserialize, Debug)]
struct TrimResult {
    read_id: String,
    trimmed_seq: Vec<u8>,
    trimmed_quals: Vec<u8>,
}

/// Performs a simple sliding window trim from both ends.
fn sliding_window_trim(
    seq: &[u8],
    quals: &[u8],
    window_size: usize,
    threshold: u8,
) -> (Vec<u8>, Vec<u8>) {
    let mut start = 0;
    let mut end = seq.len();
    let mut current_sum = 0u32;

    // Trim from the left side.
    for i in 0..seq.len() {
        current_sum += quals[i] as u32;
        if i >= window_size {
            current_sum -= quals[i - window_size] as u32;
        }
        if i >= window_size {
            let avg = current_sum as f32 / window_size as f32;
            if avg < threshold as f32 {
                start = i + 1;
            } else {
                break;
            }
        }
    }

    // Reset and trim from the right side.
    current_sum = 0;
    for i in (0..seq.len()).rev() {
        current_sum += quals[i] as u32;
        if (seq.len() - 1 - i) >= window_size {
            current_sum -= quals[i + window_size] as u32;
        }
        if (seq.len() - 1 - i) >= window_size {
            let avg = current_sum as f32 / window_size as f32;
            if avg < threshold as f32 {
                end = i;
            } else {
                break;
            }
        }
    }

    (seq[start..end].to_vec(), quals[start..end].to_vec())
}

fn main() -> Result<()> {
    let args = Args::parse();

    // Open the FASTQ (possibly gzipped) and create a reader.
    let mut reader = parse_fastx_file(&args.input)
        .with_context(|| format!("Failed to open FASTQ file at {:?}", args.input))?;

    // We'll read in chunks of records to avoid storing the entire file in memory at once.
    let mut chunk_counter = 0;
    let mut records_buffer = Vec::with_capacity(args.chunk_size);

    loop {
        // Read up to chunk_size records or until EOF.
        records_buffer.clear();
        for _ in 0..args.chunk_size {
            if let Some(record) = reader.next() {
                records_buffer.push(record?);
            } else {
                break;
            }
        }

        // If no records were loaded, we've reached EOF.
        if records_buffer.is_empty() {
            break;
        }

        // Process the chunk in parallel.
        let trim_results: Vec<TrimResult> = records_buffer
            .par_iter()
            .map(|rec| {
                let seq = rec.seq();
                let quals = rec.qual().unwrap_or(&[]);
                let (trimmed_seq, trimmed_quals) =
                    sliding_window_trim(seq, quals, args.window_size, args.quality_threshold);

                TrimResult {
                    read_id: rec.id().to_string(),
                    trimmed_seq,
                    trimmed_quals,
                }
            })
            .collect();

        // Write out partial results as JSON.
        // Example file name: partial_trim_output_chunk_0.json, partial_trim_output_chunk_1.json, etc.
        let output_file = format!("{}_chunk_{}.json", args.output_prefix, chunk_counter);
        let file = File::create(&output_file)
            .with_context(|| format!("Failed to create output file {}", output_file))?;
        let mut writer = BufWriter::new(file);
        serde_json::to_writer(&mut writer, &trim_results)
            .with_context(|| format!("Failed to serialize JSON to {}", output_file))?;
        writer.flush()?;

        println!(
            "Processed chunk {} ({} reads). Partial output: {}",
            chunk_counter,
            trim_results.len(),
            output_file
        );

        chunk_counter += 1;
    }

    println!("Completed trimming all chunks from {:?}", args.input);
    println!("Generated {} partial output files.", chunk_counter);

    Ok(())
}
