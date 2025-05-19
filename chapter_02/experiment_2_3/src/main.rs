use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use bio::io::fasta;
use ndarray::Array2;
use nalgebra::DMatrix;
use rayon::prelude::*;

fn build_de_bruijn(k: usize, sequences: &[String]) -> HashMap<String, Vec<String>> {
    sequences
        .par_iter()
        .map(|seq| {
            let mut local_map = HashMap::new();
            for i in 0..seq.len().saturating_sub(k) {
                let node = &seq[i..i + k];
                let edge = &seq[i + 1..i + k + 1];
                local_map
                    .entry(node.to_string())
                    .or_insert_with(Vec::new)
                    .push(edge.to_string());
            }
            local_map
        })
        .reduce(
            || HashMap::new(),
            |mut acc, local_map| {
                for (key, mut edges) in local_map {
                    acc.entry(key).or_insert_with(Vec::new).append(&mut edges);
                }
                acc
            },
        )
}

fn main() {
    let reader = fasta::Reader::from_file("src/reads.fasta")
        .expect("Cannot open FASTA file in 'src' directory");
    let sequences: Vec<String> = reader
        .records()
        .map(|r| {
            let record = r.expect("Invalid FASTA record");
            String::from_utf8(record.seq().to_vec()).expect("Sequence is not valid UTF-8")
        })
        .collect();

    let k = 21;
    let graph = build_de_bruijn(k, &sequences);

    let test_matrix = DMatrix::<f32>::from_element(5, 5, 1.0);
    let nd_array = Array2::<f32>::ones((5, 5));

    let output = format!(
        "Constructed De Bruijn graph with {} nodes.\n\
        nalgebra matrix: {} x {}\n\
        ndarray shape: {} x {}\n",
        graph.len(),
        test_matrix.nrows(),
        test_matrix.ncols(),
        nd_array.nrows(),
        nd_array.ncols()
    );

    // Print to terminal
    println!("{}", output);

    // Save output to a text file
    let mut file = File::create("output.txt").expect("Could not create output.txt");
    file.write_all(output.as_bytes()).expect("Could not write to output.txt");
}
