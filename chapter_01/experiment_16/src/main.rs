use rayon::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader, Result, Write};

/// A helper function to load a large genomic sequence from a FASTA file into a String.
fn load_genome(path: &str) -> Result<String> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut genome = String::new();
    for line in reader.lines() {
        let line = line?;
        if !line.starts_with('>') {
            genome.push_str(&line.trim());
        }
    }
    Ok(genome)
}

/// Constructs a naive suffix array for a given &str by sorting suffixes.
fn build_suffix_array(seq: &str) -> Vec<usize> {
    let mut suffixes: Vec<usize> = (0..seq.len()).collect();
    suffixes.sort_by(|&a, &b| seq[a..].cmp(&seq[b..]));
    suffixes
}

/// Splits the genome into `num_chunks` pieces, builds partial suffix arrays in parallel,
/// then merges them and does a final sort for a global suffix array.
fn build_parallel_suffix_array(genome: &str, num_chunks: usize) -> Vec<usize> {
    let chunk_size = (genome.len().max(1)) / num_chunks;
    let chunks: Vec<(usize, &str)> = (0..num_chunks)
        .map(|i| {
            let start = i * chunk_size;
            let end = if i == num_chunks - 1 {
                genome.len()
            } else {
                (i + 1) * chunk_size
            };
            (start, &genome[start..end])
        })
        .collect();

    let mut partial_results: Vec<(usize, Vec<usize>)> = chunks
        .par_iter()
        .map(|(offset, chunk)| {
            let sa = build_suffix_array(chunk);
            let adjusted_sa: Vec<usize> = sa.into_iter().map(|i| i + offset).collect();
            (*offset, adjusted_sa)
        })
        .collect();

    let mut combined: Vec<usize> = partial_results
        .iter_mut()
        .flat_map(|(_, sa)| sa.drain(..))
        .collect();

    combined.sort_by(|&a, &b| genome[a..].cmp(&genome[b..]));
    combined
}

fn main() -> Result<()> {
    // Define the input FASTA file and output file
    let fasta_path = "C:\\Users\\trian\\BGVR\\chapter_01\\experiment_16\\src\\example_genome.fasta";
    let output_path = "C:\\Users\\trian\\BGVR\\chapter_01\\experiment_16\\src\\output.txt";

    // Load the genomic sequence
    let genome = load_genome(fasta_path)?;

    // Build suffix array in parallel
    let num_chunks = 8;
    let suffix_array = build_parallel_suffix_array(&genome, num_chunks);

    // Save output to file
    let mut output_file = File::create(output_path)?;
    writeln!(output_file, "Genome length: {}", genome.len())?;
    writeln!(output_file, "Suffix array length: {}", suffix_array.len())?;
    writeln!(output_file, "First 10 entries in suffix array: {:?}", &suffix_array[..10.min(suffix_array.len())])?;

    println!("Processing complete! Results saved to {}", output_path);
    
    Ok(())
}
