# 1.2. Markov Random Fields

Markov Random Fields (MRFs) are graphical models that capture dependencies among variables—in this case, the nucleotides at each position in a sequence. Even a simple MRF linking consecutive nucleotides can reveal how Rust helps construct and store adjacency structures that represent these relationships. By extending this approach, you can design sophisticated models or run inference algorithms (such as belief propagation) to analyze motifs omarkor interactions in genomic data.

In the experiment, several Rust crates facilitate different tasks. The “bio” crate provides utilities for reading and writing FASTA files. The “nalgebra” and “ndarray” crates assist with linear algebra and multidimensional array operations, which are common needs for high-performance computing tasks or for computationally intensive algorithms. The “rayon” crate, known for its data-parallel functionality, can help distribute workload across multiple threads. The code will look for a FASTA file named reads.fasta in the src directory.

This code reads FASTA sequences from a file using the Rust “bio” crate, then constructs a simple Markov Random Field (MRF) by treating each position in each sequence as a graph node identified by (sequence_id, position). It links consecutive positions within a sequence (e.g., index i to index i+1) by creating edges that include a numeric “potential.” The code stores these connections in a hash map mapping each node to its neighbors. After the MRF is built, it demonstrates how to use Rust’s linear algebra crates (“nalgebra” and “ndarray”) to create and display example matrices, and prints out basic information about the constructed MRF, including a few example nodes and their edges.

Reading sequences is handled by the “bio::io::fasta” crate, which provides a convenient interface for parsing FASTA files into Rust strings. Rust enforces UTF-8 encoding, preventing many common errors that might otherwise appear when handling large sets of biological data. While the final example uses a single-threaded approach for clarity, the “rayon” crate could be used to parallelize the MRF-building loop, aggregating partial results in a concurrency-safe manner. From here, you can refine the potentials to capture more biologically realistic dependencies, apply inference algorithms for motif detection or sequence analysis, or even expand the MRF to incorporate higher-order or long-range interactions.

Files contents:
* main.rs (rust script)
* main.nf (nextflow script)
* reads.fasta (fasta file)

How to run:
cargo run main.nf (run the nextflow script that will run the main.rs and save the output in output.txt)
  
[dependencies]
bio = "2.0.3"
nalgebra = "0.33.2"
ndarray = "0.16.1"
rayon = "1.10.0"


