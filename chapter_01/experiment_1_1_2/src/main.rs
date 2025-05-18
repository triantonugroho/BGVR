use bio::io::fasta;
use std::error::Error;
use std::fs::File;
use std::io::Write;

fn main() -> Result<(), Box<dyn Error>> {
    // Open the FASTA file with the `bio::io::fasta` crate
    let reader = fasta::Reader::from_file("example.fasta")?;

    // Process sequences: filter and compute total GC content
    let total_gc: usize = reader
        .records()
        .map(|rec_res| {
            let record = rec_res.unwrap();
            String::from_utf8_lossy(record.seq()).to_string()
        })
        .filter(|seq| seq.len() >= 50)
        .map(|seq| seq.chars().filter(|&c| c == 'G' || c == 'C').count())
        .sum();

    // Write the result to output.txt
    let mut file = File::create("output.txt")?;
    writeln!(file, "Total GC content in sequences >= 50 nt: {}", total_gc)?;

    Ok(())
}
