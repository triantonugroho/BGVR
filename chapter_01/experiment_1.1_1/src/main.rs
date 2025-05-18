use bio::io::fasta;
use rayon::prelude::*;
use std::error::Error;
use std::fs::File;
use std::io::Write;

/// Reads a FASTA file and returns a vector of sequences (owned `Vec<String>`).
/// No references to local data are returned, sidestepping E0515 errors.
fn load_sequences() -> Result<Vec<String>, Box<dyn Error>> {
    let reader = fasta::Reader::from_file("example.fasta")?;
    let mut sequences = Vec::new();

    // Read each record and push into a local vector (no Mutex needed if single-thread read)
    for record_result in reader.records() {
        let record = record_result?;
        let seq_str = String::from_utf8_lossy(record.seq()).to_string();
        sequences.push(seq_str);
    }

    Ok(sequences)
}

fn main() -> Result<(), Box<dyn Error>> {
    // 1) Load all sequences from the FASTA file into an owned vector
    let seqs = load_sequences()?;

    // 2) Open file for writing
    let mut file = File::create("output.txt")?;

    // 3) Process sequences in parallel and write to file
    let results: Vec<String> = seqs.par_iter()
        .map(|seq| {
            let gc_count = seq.chars().filter(|&c| c == 'G' || c == 'C').count();
            format!("Length: {}, GC: {}", seq.len(), gc_count)
        })
        .collect();

    // Write results to file
    for result in results {
        writeln!(file, "{}", result)?;
    }

    writeln!(file, "Successfully processed {} sequences.", seqs.len())?;
    println!("Output written to output.txt");
    
    Ok(())
}
