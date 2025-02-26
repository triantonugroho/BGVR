# 1.2. Markov Random Fields

Markov Random Fields (MRFs) are graphical models that capture dependencies among variables—in this case, the nucleotides at each position in a sequence. Even a simple MRF linking consecutive nucleotides can reveal how Rust helps construct and store adjacency structures that represent these relationships. By extending this approach, you can design sophisticated models or run inference algorithms (such as belief propagation) to analyze motifs omarkor interactions in genomic data.

In the experiment, several Rust crates facilitate different tasks. The “bio” crate provides utilities for reading and writing FASTA files. The “nalgebra” and “ndarray” crates assist with linear algebra and multidimensional array operations, which are common needs for high-performance computing tasks or for computationally intensive algorithms. The “rayon” crate, known for its data-parallel functionality, can help distribute workload across multiple threads. The code will look for a FASTA file named reads.fasta in the src directory.

This code reads FASTA sequences from a file using the Rust “bio” crate, then constructs a simple Markov Random Field (MRF) by treating each position in each sequence as a graph node identified by (sequence_id, position). It links consecutive positions within a sequence (e.g., index i to index i+1) by creating edges that include a numeric “potential.” The code stores these connections in a hash map mapping each node to its neighbors. After the MRF is built, it demonstrates how to use Rust’s linear algebra crates (“nalgebra” and “ndarray”) to create and display example matrices, and prints out basic information about the constructed MRF, including a few example nodes and their edges.

Reading sequences is handled by the “bio::io::fasta” crate, which provides a convenient interface for parsing FASTA files into Rust strings. Rust enforces UTF-8 encoding, preventing many common errors that might otherwise appear when handling large sets of biological data. While the final example uses a single-threaded approach for clarity, the “rayon” crate could be used to parallelize the MRF-building loop, aggregating partial results in a concurrency-safe manner. From here, you can refine the potentials to capture more biologically realistic dependencies, apply inference algorithms for motif detection or sequence analysis, or even expand the MRF to incorporate higher-order or long-range interactions.

## Files contents:
* main.rs (rust script)
* main.nf (nextflow script)
* reads.fasta (fasta file)

## How to run:

cargo run main.nf (run the nextflow script that will run the main.rs and save the output in output.txt)

## Explanation of the Output:

The output provides details about the constructed Markov Random Field (MRF) based on the input sequences from the FASTA file. Here's a breakdown of the results:

1. Total Number of Nodes

Constructed MRF with 251 nodes
The MRF consists of 251 nodes, which represent positions in the sequences from the FASTA file.
Nodes are indexed using a (sequence_id, position) format, where:
sequence_id is the index of the sequence in the input list.
position represents a specific character's index within that sequence.

2. Matrix Dimensions

nalgebra matrix dimensions: 10 x 10
ndarray matrix dimensions: 10 x 10

* Two matrices of size 10 × 10 were created:
  * Nalgebra's DMatrix
  * Ndarray's Array2
* These matrices are initialized with values, but their purpose is not directly related to the MRF construction.

3. Sample Nodes and Their Connections

Node (4, 41) has edges:
  -> (4, 42) with potential = 1
Node (0, 14) has edges:
  -> (0, 15) with potential = 1
Node (0, 5) has edges:
  -> (0, 6) with potential = 1
Node (0, 46) has edges:
  -> (0, 47) with potential = 1
Node (4, 35) has edges:
  -> (4, 36) with potential = 1

* Each node is a (sequence_id, position) tuple.
* The output shows 5 nodes and their edges.
* Each node (i, j) connects to (i, j+1), forming a simple chain structure.
* The potential value of 1.0 suggests a uniform edge weight, indicating equal transition probabilities between consecutive positions in the sequences.

## Conclusion

The program successfully builds an MRF from the input FASTA sequences, representing positional relationships in the data. The adjacency structure follows a sequential connection pattern, ensuring each character in a sequence is linked to its next character. The generated output file (output.txt) records this structure along with matrix dimensions from nalgebra and ndarray.

## [dependencies]

bio = "2.0.3"

nalgebra = "0.33.2"

ndarray = "0.16.1"

rayon = "1.10.0"


