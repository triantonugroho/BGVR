use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::Mutex;
use std::fs::File;
use std::io::{Write, BufWriter};

use anyhow::{Context, Result};
use bio::io::fasta;
use rayon::prelude::*;
use log::{info, warn};
use env_logger::{Builder, Target};

#[derive(Debug)]
struct DeBruijnGraph {
    edges: HashMap<String, Vec<String>>,
}

impl DeBruijnGraph {
    fn new() -> Self {
        DeBruijnGraph { edges: HashMap::new() }
    }
    fn add_edge(&mut self, prefix: &str, suffix: &str) {
        self.edges
            .entry(prefix.to_string())
            .or_insert_with(Vec::new)
            .push(suffix.to_string());
    }
    fn merge(&mut self, other: DeBruijnGraph) {
        for (pref, suff_list) in other.edges {
            self.edges
                .entry(pref)
                .or_insert_with(Vec::new)
                .extend(suff_list);
        }
    }
}

fn build_partial_graph(reads: &[Vec<u8>], k: usize) -> DeBruijnGraph {
    let mut local_graph = DeBruijnGraph::new();
    for seq in reads {
        let s = String::from_utf8_lossy(seq);
        if s.len() < k {
            continue;
        }
        for i in 0..=s.len() - k {
            let prefix = &s[i..i + k - 1];
            let suffix = &s[i + 1..i + k];
            local_graph.add_edge(prefix, suffix);
        }
    }
    local_graph
}

fn main() -> Result<()> {
    // Open a file to log output (output_new.txt)
    let log_file = File::create("output_new.txt").with_context(|| "Failed to create log file")?;
    let mut log_writer = BufWriter::new(log_file); // Using BufWriter for more efficient file writing

    // Set up the logger to log to both stdout and the log file, without timestamp
    Builder::new()
        .format(|buf, record| {
            writeln!(buf, "{}", record.args()) // Only log the message (no timestamp, level, etc.)
        })
        .filter_level(log::LevelFilter::Info) // Ensure logging is at 'Info' level
        .target(Target::Stdout)     // Log to the console
        .init();

    // Hardcoded file path for input FASTA
    let input_path = PathBuf::from(r"C:\Users\trian\BGVR\chapter_02\experiment_23_2\src\reads.fasta");
    let k = 21; // Default k-mer size

    // Log the k-mer size for confirmation
    info!("Building de Bruijn graph with k-mer size = {}", k);

    // Try to read the FASTA file
    let reader = fasta::Reader::from_file(&input_path)
        .with_context(|| format!("Could not open FASTA file: {:?}", input_path))?;

    // Collect the result of reading records into a Result<Vec<Vec<u8>>>
    let reads: Vec<Vec<u8>> = reader
        .records()
        .map(|r| {
            let record = r.context("Error reading FASTA record.")?;
            info!("Read record of length: {}", record.seq().len()); // Log the length of each sequence
            Ok(record.seq().to_owned())
        })
        .collect::<Result<Vec<Vec<u8>>, anyhow::Error>>()?;

    info!("Total reads found: {}", reads.len());

    if reads.is_empty() {
        warn!("No reads found; exiting...");
        return Ok(());
    }

    let final_graph = Mutex::new(DeBruijnGraph::new());
    let chunk_size = 1000;

    // Process reads in parallel
    reads
        .par_chunks(chunk_size)
        .map(|chunk| build_partial_graph(chunk, k))
        .for_each(|partial| {
            let mut global = final_graph.lock().unwrap();
            global.merge(partial);
        });

    let final_graph = final_graph.into_inner().unwrap();
    let graph_info = format!("Built graph with {} (k-1)-mer prefixes", final_graph.edges.len());

    // Log the graph info to both stdout and output_new.txt
    info!("{}", graph_info);
    writeln!(log_writer, "{}", graph_info)?;

    Ok(())
}
