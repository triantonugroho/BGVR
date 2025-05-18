use rayon::prelude::*;
use needletail::{parse_fastx_file};
use dashmap::DashMap;
use serde::{Serialize, Deserialize};
use std::sync::Arc;
use std::path::PathBuf;
use std::fs;
use std::error::Error;

/// A record capturing the final k-mer counts for a specific k.
#[derive(Debug, Serialize, Deserialize)]
struct KmerCountResult {
    k: usize,
    counts: Vec<(String, usize)>,
}

fn main() -> Result<(), Box<dyn Error>> {
    // 1) Read command-line arguments: path to FASTQ/FASTA and desired k.
    let input_path = std::env::args().nth(1).unwrap_or_else(|| "reads.fq".to_string());
    let k: usize = std::env::args()
        .nth(2)
        .unwrap_or_else(|| "31".to_string())
        .parse()
        .map_err(|e| format!("Failed to parse k: {}", e))?;

    // 2) Create a concurrency-friendly map for final merges.
    //    We'll store (k-mer -> count).
    let kmer_map = Arc::new(DashMap::<String, usize>::new());

    // 3) Attempt to read the FASTA/FASTQ file.
    let mut reader = parse_fastx_file(PathBuf::from(&input_path))
        .map_err(|e| format!("Error parsing input file: {}", e))?;
    
    // 4) Store reads in-memory (for demonstration). In HPC usage,
    //    you might chunk or stream to handle large datasets.
    let mut seq_store = Vec::new();
    while let Some(record) = reader.next() {
        let rec = record?;
        seq_store.push(rec.seq().to_vec());
    }

    // 5) For each read in parallel, build a local HashMap of k-mer counts,
    //    then merge into the DashMap. This approach avoids partial locking
    //    on every k-mer insert, boosting throughput for large data.
    seq_store.par_iter().for_each(|seq| {
        let mut local_counts = std::collections::HashMap::<String, usize>::new();
        
        // Generate all k-mers from the read.
        for i in 0..=(seq.len().saturating_sub(k)) {
            let fragment = &seq[i..i + k];
            let kmer_str = String::from_utf8_lossy(fragment).to_string();
            *local_counts.entry(kmer_str).or_insert(0) += 1;
        }

        // Merge local counts into the global DashMap
        for (kmer, count) in local_counts {
            *kmer_map.entry(kmer).or_insert(0) += count;
        }
    });

    // 6) Extract data from the DashMap and sort by k-mer lex order.
    //    This step is single-threaded, but typically only needed once.
    let mut results: Vec<_> = kmer_map
        .iter()
        .map(|entry| (entry.key().clone(), *entry.value()))
        .collect();
    results.sort_by(|a, b| a.0.cmp(&b.0));

    // 7) Serialize k-mer counts as JSON, writing to a file for downstream use.
    let out = KmerCountResult { k, counts: results };
    fs::write("kmer_counts.json", serde_json::to_string_pretty(&out)?)?;

    println!("Wrote k-mer counts to kmer_counts.json. Done!");
    Ok(())
}