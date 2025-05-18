## 2.3. De Bruijn Graphs

### experiment_2.3

Below is a Rust code that demonstrates building a simple De Bruijn graph from a FASTA file by extracting overlapping k-mers and linking them, thereby avoiding pairwise alignments for each read. De Bruijn graphs are commonly employed in genome assembly because they can efficiently represent overlapping k-mers, simplifying reconstruction. By keeping memory usage explicit and parallelizing iteration, the program remains robust and scalable to large data sets containing millions of reads. It leverages Rust’s strong memory safety and concurrency features, along with crates like “rayon” for parallel processing and “nalgebra”/“ndarray” for numeric tasks, ensuring efficient performance in high-throughput sequencing scenarios.

This program imports crates for bioinformatics (bio), linear algebra (nalgebra, ndarray), and parallel computation (rayon). In the build_de_bruijn function, each sequence is processed by extracting overlapping k-mers of length k (referred to as node), along with the subsequent overlapping k-mer (edge). Rather than modifying a single shared hash map in multiple threads, each thread accumulates its own local map of k-mer relationships. Rayon’s .par_iter() enables concurrent iteration over the sequences, and a reduce operation merges these local maps into a single global HashMap<String, Vec<String>>.

After reading a FASTA file named reads.fasta from the src directory, the main function collects all sequences into a vector and invokes build_de_bruijn with a chosen k-mer size (k = 21). Once the De Bruijn graph is built, the code demonstrates HPC-oriented crates by creating example matrices with “nalgebra” and “ndarray.” Finally, it prints out details about the graph’s size and the matrices, confirming that the parallel construction and supporting data structures have been set up correctly.

Several success stories highlight substantial reductions in runtime and cost once legacy Python or Java components are rewritten in Rust, particularly for k-mer counting or parallel motif searches. By combining HPC scheduling with containerized Rust executables, these organizations accelerate the pace of biomarker discovery and gene therapy research while preserving the reproducibility needed for regulatory compliance and collaborative studies.

#### Project Structure:

```plaintext
experiment_2.3/
├── Cargo.toml                     # Rust project configuration and dependencies
└── src/
    ├── main.rs                    # Main Rust script containing program logic
    ├── reads.fasta                # FASTA file containing sequence reads
    └── output.txt                 # Output file
```

#### Cargo.toml

```toml
[package]
name = "experiment_2.3"
version = "0.1.0"
edition = "2021"

[dependencies]
bio = "2.0.3"
nalgebra = "0.33.2"
ndarray = "0.16.1"
rayon = "1.10.0"
```

#### How to run

```powershell
cargo run | tee output.txt
```

(run main.rs and save the output in output.txt)

#### Explanation of the Output

The output describes the results of constructing a De Bruijn graph from the input sequences stored in a FASTA file. Below is a breakdown of each part of the output:

##### 1. De Bruijn Graph Construction

* Constructed De Bruijn graph with 150 nodes.
* The program reads FASTA sequences and constructs a De Bruijn graph using a k-mer length of 21 (k = 21).
* A total of 150 nodes were created.
* Each node represents a k-mer (a substring of length k from the sequences).
* Edges represent transitions between overlapping k-mers.

##### 2. Matrix Dimensions

nalgebra matrix: 5 x 5

ndarray shape: 5 x 5

* Two matrices of size 5 × 5 were created:
  * nalgebra::DMatrix (a dynamically sized matrix from the nalgebra library).
  * ndarray::Array2 (a 2D array from the ndarray library).
* Both matrices are initialized with values (all elements are 1.0), but they are not directly related to the graph construction.

##### 3. Output Storage

* The program prints the results to the terminal.
* It also saves the output to a text file (output.txt) for further analysis.

#### Conclusion

The program successfully constructs a parallelized De Bruijn graph using Rayon for efficiency. The graph consists of 150 nodes, each representing a k-mer. Additionally, it creates two example 5 × 5 matrices for demonstration purposes. The results are both displayed on the terminal and stored in output.txt.

