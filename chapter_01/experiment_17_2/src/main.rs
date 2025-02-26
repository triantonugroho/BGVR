use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use ndarray::Array2; // 's' removed, we keep 'Array2' just for demonstration.
use nalgebra::DMatrix;
use rayon::prelude::*;

fn gene_gnn_iteration(
    adjacency: &HashMap<usize, Vec<usize>>,
    node_features: &DMatrix<f32>,
) -> DMatrix<f32> {
    let n = node_features.nrows();
    let dim = node_features.ncols();
    let mut updated = DMatrix::<f32>::zeros(n, dim);

    // Instead of `par_row_iter_mut()`, we operate on the underlying slice
    // and process chunks of size 'dim' in parallel.
    updated
        .as_mut_slice()
        .par_chunks_mut(dim)
        .enumerate()
        .for_each(|(i, row_buf)| {
            if let Some(neighbors) = adjacency.get(&i) {
                // Accumulate neighbor features
                let sum_features = neighbors.iter().fold(vec![0.0; dim], |mut acc, &nbr| {
                    for d in 0..dim {
                        acc[d] += node_features[(nbr, d)];
                    }
                    acc
                });

                // Average neighbor features into the row buffer
                for d in 0..dim {
                    row_buf[d] = sum_features[d] / (neighbors.len() as f32 + 1e-8);
                }
            }
        });

    updated
}

fn main() -> std::io::Result<()> {
    let output_path = "C:\\Users\\trian\\BGVR\\chapter_01\\experiment_17_2\\src\\output.txt";
    let mut file = File::create(output_path)?;

    // Example adjacency and expression data
    let adjacency: HashMap<usize, Vec<usize>> = HashMap::new();
    let expression_data = DMatrix::<f32>::from_element(100, 10, 1.0);
    let adjacency_nd = Array2::<f32>::ones((100, 100));

    writeln!(
        file,
        "Expression data shape: {} x {}",
        expression_data.nrows(),
        expression_data.ncols()
    )?;
    writeln!(
        file,
        "Adjacency array shape: {} x {}",
        adjacency_nd.nrows(),
        adjacency_nd.ncols()
    )?;

    let updated = gene_gnn_iteration(&adjacency, &expression_data);

    writeln!(
        file,
        "Updated node features shape: {} x {}",
        updated.nrows(),
        updated.ncols()
    )?;

    // Save updated node features to file
    for i in 0..updated.nrows() {
        let row: Vec<String> = (0..updated.ncols())
            .map(|j| format!("{:.6}", updated[(i, j)]))
            .collect();
        writeln!(file, "{}", row.join("\t"))?;
    }

    println!("Output saved to {}", output_path);

    Ok(())
}
