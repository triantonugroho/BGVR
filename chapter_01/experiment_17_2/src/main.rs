use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use rand::{Rng, rng};
use nalgebra::DMatrix;
use rayon::prelude::*;

fn generate_adjacency_list(n: usize, edges_per_node: usize) -> HashMap<usize, Vec<usize>> {
    let mut adjacency: HashMap<usize, Vec<usize>> = HashMap::new();
    let mut rng = rng(); // Updated

    for i in 0..n {
        let mut neighbors = Vec::new();
        while neighbors.len() < edges_per_node {
            let neighbor = rng.random_range(0..n); // Updated
            if neighbor != i && !neighbors.contains(&neighbor) {
                neighbors.push(neighbor);
            }
        }
        adjacency.insert(i, neighbors);
    }

    adjacency
}

fn gene_gnn_iteration(
    adjacency: &HashMap<usize, Vec<usize>>,
    node_features: &DMatrix<f32>,
) -> DMatrix<f32> {
    let n = node_features.nrows();
    let dim = node_features.ncols();
    let mut updated = DMatrix::<f32>::zeros(n, dim);

    updated
        .as_mut_slice()
        .par_chunks_mut(dim)
        .enumerate()
        .for_each(|(i, row_buf)| {
            if let Some(neighbors) = adjacency.get(&i) {
                if neighbors.is_empty() {
                    // No neighbors, retain original features
                    row_buf.copy_from_slice(node_features.row(i).clone_owned().as_slice());
                } else {
                    let mut sum_features = vec![0.0; dim];
                    for &nbr in neighbors {
                        for d in 0..dim {
                            sum_features[d] += node_features[(nbr, d)];
                        }
                    }
                    let degree = neighbors.len() as f32;
                    for d in 0..dim {
                        row_buf[d] = sum_features[d] / degree;
                    }
                }
            } else {
                // If no entry in adjacency list, retain original features
                row_buf.copy_from_slice(node_features.row(i).clone_owned().as_slice());
            }
        });

    updated
}

fn main() -> std::io::Result<()> {
    let output_path = "C:\\Users\\trian\\BGVR\\chapter_01\\experiment_17_2\\src\\output.txt";
    let file = File::create(output_path)?;
    let mut writer = BufWriter::new(file);

    let num_nodes = 100;
    let feature_dim = 10;
    let edges_per_node = 5;

    let adjacency = generate_adjacency_list(num_nodes, edges_per_node);
    
    // Randomly initialize node features
    let mut rng = rng(); // Updated
    let expression_data = DMatrix::<f32>::from_fn(num_nodes, feature_dim, |_, _| rng.random()); // Updated

    writeln!(
        writer,
        "Expression data shape: {} x {}",
        expression_data.nrows(),
        expression_data.ncols()
    )?;

    let updated = gene_gnn_iteration(&adjacency, &expression_data);

    writeln!(
        writer,
        "Updated node features shape: {} x {}",
        updated.nrows(),
        updated.ncols()
    )?;

    for i in 0..updated.nrows() {
        let row: Vec<String> = (0..updated.ncols())
            .map(|j| format!("{:.6}", updated[(i, j)]))
            .collect();
        writeln!(writer, "{}", row.join("\t"))?;
    }

    println!("Output saved to {}", output_path);
    Ok(())
}
