use std::process::{Command, exit};
use std::fs::{File, rename};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path};
use std::collections::HashMap;

fn main() {
    // Set the input and output paths directly inside the code
    let input_path = "SRR11192680".to_string(); // SRA ID
    let output_path = "C:\\Users\\trian\\BGVR\\chapter_01\\experiment_18_1_2\\src\\de_bruijn_graph.txt".to_string(); // Output file path

    // Set k and sra_toolkit_path as well
    let k = 21;
    let sra_toolkit_path = "C:\\SRAToolkit\\bin\\fasterq-dump.exe".to_string(); // Default path

    // Print input and output paths to verify
    println!("Input Path: {}", input_path);
    println!("Output Path: {}", output_path);

    // Download the FASTQ file using fasterq-dump
    println!("Downloading SRA file {}...", input_path);
    let output_dir = Path::new(&output_path).parent().unwrap_or_else(|| Path::new("."));
    let sra_toolkit_path = Path::new(&sra_toolkit_path);

    let status = Command::new(sra_toolkit_path)
        .arg(&input_path)
        .arg("--outdir")
        .arg(output_dir)
        .arg("--split-files")
        .status()
        .expect("Failed to execute fasterq-dump");

    if !status.success() {
        eprintln!("Error downloading FASTQ file for SRA ID {}", input_path);
        exit(1);
    }

    // Assuming the file is downloaded successfully
    let fastq_file = output_dir.join(format!("{}_1.fastq", input_path));
    let fastq_file_renamed = "reads.fastq"; // Rename to a standard name for further processing

    // Rename the downloaded file (if needed)
    rename(fastq_file, fastq_file_renamed).unwrap_or_else(|_| {
        eprintln!("Failed to rename downloaded FASTQ file.");
        exit(1);
    });

    println!("Download complete, file saved as {}", fastq_file_renamed);

    // Read the FASTQ file and gather sequences
    let file = File::open(fastq_file_renamed).unwrap_or_else(|_| {
        eprintln!("Error opening FASTQ file: {}", fastq_file_renamed);
        exit(1);
    });
    let reader = BufReader::new(file);

    let mut sequences = Vec::new();
    let mut line_count = 0;

    // Read FASTQ sequences
    for line_res in reader.lines() {
        let line = line_res.unwrap();
        line_count += 1;
        if line_count % 4 == 2 {
            sequences.push(line);
        }
    }

    // Build De Bruijn graph
    let graph = build_de_bruijn_graph(k, &sequences);

    // Write the graph to an output file
    let mut outfile = File::create(&output_path).unwrap_or_else(|_| {
        eprintln!("Error creating output file: {}", output_path);
        exit(1);
    });

    writeln!(outfile, "De Bruijn Graph (k={})", k).unwrap();
    writeln!(outfile, "Number of nodes: {}", graph.len()).unwrap();

    for (node, edges) in graph {
        writeln!(outfile, "{} => {:?}", node, edges).unwrap();
    }

    println!("De Bruijn graph construction complete. Output written to {}", output_path);
}

fn build_de_bruijn_graph(k: usize, sequences: &[String]) -> HashMap<String, Vec<String>> {
    let mut graph: HashMap<String, Vec<String>> = HashMap::new();

    for seq in sequences {
        if seq.len() < k {
            continue;
        }
        for i in 0..seq.len().saturating_sub(k) {
            let node = &seq[i..i + k];
            let edge = &seq[i + 1..i + k + 1];
            let entry = graph.entry(node.to_string()).or_insert_with(Vec::new);
            entry.push(edge.to_string());
        }
    }

    graph
}
