use anyhow::{Context, Result};
use rayon::prelude::*;
use serde::{Serialize, Deserialize};
use std::collections::HashMap;
use std::fs::{File, create_dir_all};
use std::io::{BufWriter, BufReader, Write};
use std::path::PathBuf;
use clap::Parser;

/// Represents a segment of a read alignment, including minimal CIGAR data.
#[derive(Debug, Clone, Serialize, Deserialize)]
struct AlignmentSegment {
    read_id: String,
    chrom: String,
    start: u64,
    cigar: String,
    orientation: char, // '+', '-' or more complex
}

/// Represents a detected breakpoint in naive structural variant analysis.
#[derive(Debug, Serialize, Deserialize)]
struct Breakpoint {
    read_id: String,
    chrom: String,
    pos: u64,
    sv_type: String,
}

/// Command-line arguments for chunk-based read parsing and breakpoint detection.
#[derive(Parser, Debug)]
#[command(name = "sv_detector")]
#[command(about = "Naive split-read based structural variant detector in Rust")]
struct Args {
    /// Input file containing alignment segments in JSON format.
    #[arg(long)]
    alignment_input: PathBuf,

    /// Number of alignment segments to load per chunk.
    #[arg(long, default_value_t = 5_000)]
    chunk_size: usize,

    /// Directory to which partial breakpoint results will be written.
    #[arg(long, default_value = "partial_breakpoints")]
    partial_output_dir: PathBuf,

    /// Final merged breakpoint results JSON file.
    #[arg(long, default_value = "merged_breakpoints.json")]
    merged_output: PathBuf,
}

/// A partial container holding a batch of breakpoints.
#[derive(Debug, Serialize, Deserialize)]
struct PartialBreakpoints {
    breakpoints: Vec<Breakpoint>,
}

/// Naive parser for the numeric length in a CIGAR string (e.g., "50M" => 50).
fn parse_cigar_len(cigar: &str) -> u64 {
    let len_str = cigar.trim_end_matches(|c: char| !c.is_numeric());
    len_str.parse::<u64>().unwrap_or(0)
}

/// Detects breakpoints by examining consecutive segments from the same read.
fn detect_breakpoints(segments: &[AlignmentSegment]) -> Vec<Breakpoint> {
    // If a single read has multiple alignment segments, we consider naive breakpoints.
    if segments.len() < 2 {
        return Vec::new();
    }
    let mut bps = Vec::new();
    for window in segments.windows(2) {
        let first = &window[0];
        let second = &window[1];
        if first.chrom == second.chrom {
            // Potential 'intra-chromosomal' breakpoint where the first alignment ends
            let breakpos = first.start + parse_cigar_len(&first.cigar);
            bps.push(Breakpoint {
                read_id: first.read_id.clone(),
                chrom: first.chrom.clone(),
                pos: breakpos,
                sv_type: "intra-chr".to_string(),
            });
        } else {
            // If on different chromosomes, a naive translocation guess
            bps.push(Breakpoint {
                read_id: first.read_id.clone(),
                chrom: first.chrom.clone(),
                pos: first.start,
                sv_type: "translocation".to_string(),
            });
        }
    }
    bps
}

/// Merges two vectors of breakpoints by concatenation.
fn merge_breakpoints(mut global: Vec<Breakpoint>, mut local: Vec<Breakpoint>) -> Vec<Breakpoint> {
    global.append(&mut local);
    global
}

/// Reads a specified number of AlignmentSegment records from a JSON source.
fn read_chunk<R: std::io::Read>(
    reader: &mut serde_json::StreamDeserializer<'_, serde_json::de::IoRead<R>, AlignmentSegment>,
    chunk_size: usize,
) -> Result<Vec<AlignmentSegment>> {
    let mut chunk = Vec::with_capacity(chunk_size);
    for _ in 0..chunk_size {
        if let Some(res) = reader.next() {
            let record = res?;
            chunk.push(record);
        } else {
            break;
        }
    }
    Ok(chunk)
}

fn main() -> Result<()> {
    let args = Args::parse();

    // Ensure the directory for partial outputs exists.
    create_dir_all(&args.partial_output_dir)
        .with_context(|| format!("Failed to create partial output directory {:?}", args.partial_output_dir))?;

    // Open the alignment JSON file and prepare a streaming deserializer.
    let file = File::open(&args.alignment_input)
        .with_context(|| format!("Failed to open alignment input file {:?}", args.alignment_input))?;
    let reader = BufReader::new(file);
    let mut stream = serde_json::Deserializer::from_reader(reader).into_iter::<AlignmentSegment>();

    let mut chunk_index = 0usize;
    loop {
        // Read a chunk of alignment records.
        let records = read_chunk(&mut stream, args.chunk_size)?;
        if records.is_empty() {
            break;
        }

        // Group alignments by read ID.
        let mut read_map: HashMap<String, Vec<AlignmentSegment>> = HashMap::new();
        for seg in records {
            read_map.entry(seg.read_id.clone()).or_default().push(seg);
        }

        // Detect breakpoints in parallel.
        let partial_breakpoints: Vec<Breakpoint> = read_map
            .par_iter()
            .flat_map(|(_read_id, segs)| {
                // Sort segments by start for consistent ordering
                let mut ordered = segs.clone();
                ordered.sort_by_key(|s| s.start);
                detect_breakpoints(&ordered)
            })
            .collect();

        // Serialize these partial breakpoints to disk.
        let partial_res = PartialBreakpoints { breakpoints: partial_breakpoints };
        let chunk_path = args.partial_output_dir.join(format!("partial_breakpoints_{}.json", chunk_index));
        let chunk_file = File::create(&chunk_path)
            .with_context(|| format!("Failed to create partial breakpoint file {:?}", chunk_path))?;
        serde_json::to_writer(BufWriter::new(chunk_file), &partial_res)
            .with_context(|| format!("Failed to write partial breakpoint data to {:?}", chunk_path))?;

        println!(
            "Processed chunk {} ({} read groups). Wrote partial results to {:?}",
            chunk_index, read_map.len(), chunk_path
        );
        chunk_index += 1;
    }

    // Now merge all partial outputs into a single file.
    let dir_entries = std::fs::read_dir(&args.partial_output_dir)
        .with_context(|| format!("Failed to read partial output directory {:?}", args.partial_output_dir))?;

    let mut merged = Vec::new();
    for entry in dir_entries {
        let path = entry?.path();
        if path.file_name().map_or(false, |p| p.to_string_lossy().starts_with("partial_breakpoints_")) {
            let file = File::open(&path)
                .with_context(|| format!("Failed to open partial breakpoints file {:?}", path))?;
            let partial: PartialBreakpoints = serde_json::from_reader(BufReader::new(file))
                .with_context(|| format!("Failed to parse partial breakpoints from {:?}", path))?;
            merged = merge_breakpoints(merged, partial.breakpoints);
        }
    }

    let merged_file = File::create(&args.merged_output)
        .with_context(|| format!("Failed to create merged output file {:?}", args.merged_output))?;
    serde_json::to_writer(BufWriter::new(merged_file), &merged)
        .with_context(|| format!("Failed to write merged breakpoints to {:?}", args.merged_output))?;

    println!(
        "Merged {} total breakpoints into {:?}.",
        merged.len(),
        args.merged_output
    );

    Ok(())
}
