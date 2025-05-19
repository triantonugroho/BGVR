use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;
use clap::{Parser};

/// A simple program to preprocess FASTQ reads
#[derive(Parser)]
#[command(name = "rust_preprocess")]
#[command(version = "1.0")]
#[command(about = "Preprocesses FASTQ reads", long_about = None)]
struct Cli {
    /// Input FASTQ file
    #[arg(short, long)]
    input: PathBuf,

    /// Output FASTQ file
    #[arg(short, long)]
    output: PathBuf,
}

fn main() -> std::io::Result<()> {
    let args = Cli::parse();

    // Open input file
    let input_file = File::open(&args.input)?;
    let reader = BufReader::new(input_file);

    // Create output file
    let output_file = File::create(&args.output)?;
    let mut writer = BufWriter::new(output_file);

    // Dummy passthrough: reads are copied as-is.
    // Insert your filtering logic here.
    for line in reader.lines() {
        let line = line?; // Unwrap line safely
        writeln!(writer, "{}", line)?;
    }

    println!("Preprocessing complete. Output saved to {:?}", args.output);
    Ok(())
}
