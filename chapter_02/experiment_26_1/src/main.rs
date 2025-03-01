use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = std::env::args().collect();

    // ✅ Ensure at least 5 arguments are provided
    if args.len() < 5 {
        eprintln!("Usage: {} --input <input.fastq> --output <output.fastq>", args[0]);
        std::process::exit(1);
    }

    let input_path = &args[2];
    let output_path = &args[4];

    let infile = File::open(input_path)?;
    let mut outfile = File::create(output_path)?;
    let reader = BufReader::new(infile);
    let mut lines = reader.lines();

    while let Some(Ok(header)) = lines.next() {
        let seq = lines.next().unwrap_or(Ok(String::new()))?;
        let plus_line = lines.next().unwrap_or(Ok(String::new()))?;
        let qual = lines.next().unwrap_or(Ok(String::new()))?;

        let trimmed_seq = trim_bases(&seq, 5);
        let trimmed_qual = trim_bases(&qual, 5);

        writeln!(outfile, "{}", header)?;
        writeln!(outfile, "{}", trimmed_seq)?;
        writeln!(outfile, "{}", plus_line)?;
        writeln!(outfile, "{}", trimmed_qual)?;
    }

    println!("Trimming complete: {} -> {}", input_path, output_path);
    Ok(())
}

fn trim_bases(s: &str, n: usize) -> String {
    if s.len() < 2 * n {
        s.to_string()
    } else {
        s[n..(s.len() - n)].to_string()
    }
}
