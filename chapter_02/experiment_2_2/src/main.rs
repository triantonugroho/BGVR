use std::collections::HashMap;
use std::fs::File;
use std::io::{Write, BufWriter};
use bio::io::fasta;
use ndarray::Array2;
use nalgebra::DMatrix;

#[derive(Debug)]
struct MRFEdge {
    potential: f32,
}

#[derive(Debug)]
struct MRF {
    adjacency: HashMap<(usize, usize), Vec<((usize, usize), MRFEdge)>>,
}

fn build_mrf_from_sequences(sequences: Vec<String>) -> MRF {
    let mut adjacency: HashMap<(usize, usize), Vec<((usize, usize), MRFEdge)>> = HashMap::new();

    for (seq_id, seq) in sequences.iter().enumerate() {
        for i in 0..seq.len().saturating_sub(1) {
            let node_a = (seq_id, i);
            let node_b = (seq_id, i + 1);
            let edge = MRFEdge { potential: 1.0 };

            adjacency
                .entry(node_a)
                .or_insert_with(Vec::new)
                .push((node_b, edge));
        }
    }

    MRF { adjacency }
}

fn main() {
    let reader = fasta::Reader::from_file("src/reads.fasta")
        .expect("Error opening FASTA in the 'src' directory");

    let seqs: Vec<String> = reader
        .records()
        .map(|rec| {
            let r = rec.expect("Invalid record");
            String::from_utf8(r.seq().to_vec()).expect("Invalid UTF-8 sequence")
        })
        .collect();

    let mrf = build_mrf_from_sequences(seqs);

    let matrix_nnalgebra = DMatrix::<f32>::from_element(10, 10, 1.0);
    let matrix_ndarray = Array2::<f32>::ones((10, 10));

    // Open the output file for writing
    let file = File::create("output.txt").expect("Unable to create file");
    let mut writer = BufWriter::new(file);

    writeln!(writer, "Constructed MRF with {} nodes", mrf.adjacency.len()).unwrap();
    writeln!(writer, "nalgebra matrix dimensions: {} x {}", matrix_nnalgebra.nrows(), matrix_nnalgebra.ncols()).unwrap();
    writeln!(writer, "ndarray matrix dimensions: {} x {}", matrix_ndarray.nrows(), matrix_ndarray.ncols()).unwrap();

    for (node, edges) in mrf.adjacency.iter().take(5) {
        writeln!(writer, "Node {:?} has edges:", node).unwrap();
        for (neighbor, edge) in edges {
            writeln!(writer, "  -> {:?} with potential = {}", neighbor, edge.potential).unwrap();
        }
    }

    println!("Output written to output.txt");
}
