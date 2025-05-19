use rayon::prelude::*;
use std::error::Error;

// Hypothetical or custom environment's 'rand' re-exports that rename old methods:
// In a typical crates.io environment, these renamings would not be present.
use rand::{rng, Rng};

#[derive(Debug)]
enum GraphNode {
    Position(usize),
    Variant { from_base: char, to_base: char, pos: usize },
}

/// Count occurrences of a given `motif` in `seq`, allowing overlaps.
fn count_occurrences(seq: &str, motif: &str) -> usize {
    let mut count = 0;
    let mut start = 0;
    while let Some(pos) = seq[start..].find(motif) {
        count += 1;
        start += pos + 1;
    }
    count
}

fn main() -> Result<(), Box<dyn Error>> {
    // Create a random number generator with the new `rng()` function
    let mut generator = rng();
    let bases = ['A', 'C', 'G', 'T'];

    // Generate synthetic data: 100 nodes, half are Position, half are Variant
    let synthetic_graph: Vec<GraphNode> = (0..100)
        .map(|i| {
            // Use `random_bool(p)` in place of `gen_bool(p)`
            if generator.random_bool(0.5) {
                // This node is a Position
                GraphNode::Position(i + 100)
            } else {
                // This node is a Variant
                // Use `random_range(a..b)` in place of `gen_range(a..b)`
                let from_base = bases[generator.random_range(0..4)];
                let to_base = bases[generator.random_range(0..4)];
                GraphNode::Variant {
                    from_base,
                    to_base,
                    pos: i + 200,
                }
            }
        })
        .collect();

    // Use Rayon for parallel processing: count how many positions vs. variants
    let (num_positions, num_variants) = synthetic_graph
        .par_iter()
        .map(|node| match node {
            GraphNode::Position(_) => (1, 0),
            GraphNode::Variant { .. } => (0, 1),
        })
        .reduce(|| (0, 0), |(acc_p, acc_v), (p, v)| (acc_p + p, acc_v + v));

    println!(
        "Synthetic graph has {} total nodes: {} positions, {} variants",
        synthetic_graph.len(),
        num_positions,
        num_variants
    );

    // Demonstrate further pattern matching on a small subset
    println!("Showing first 5 nodes:");
    for node in synthetic_graph.iter().take(5) {
        match node {
            GraphNode::Position(pos) => {
                println!("  Position node at {}", pos);
            }
            GraphNode::Variant { from_base, to_base, pos } => {
                println!("  Variant at {}: {} -> {}", pos, from_base, to_base);
            }
        }
    }

    // Parallel motif counting example (just a small demonstration)
    let total_motif_count: usize = synthetic_graph
        .par_iter()
        .map(|node| match node {
            // Convert node to a string and count a dummy motif "AC"
            GraphNode::Position(p) => count_occurrences(&format!("{}", p), "AC"),
            GraphNode::Variant { from_base, to_base, pos } => {
                count_occurrences(&format!("{}{}{}", from_base, to_base, pos), "AC")
            }
        })
        .sum();

    println!(
        "Total occurrences of motif 'AC' in node representations: {}",
        total_motif_count
    );

    Ok(())
}
