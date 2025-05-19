## 3.1. Introduction to Data Structures and Algorithms

### experiment_3_1_3

This code demonstrates a scenario where different types of genomic data—such as microbial, eukaryotic, pangenome, single-cell transcriptomic, and Hi-C assays—require different indexing approaches. It defines a trait-based framework that can build either linear or graph-based indexes depending on the genomic context, then uses Rayon to construct each index in parallel. The result is a scalable infrastructure where large datasets, each with unique characteristics, can be assigned the most suitable index structure and processed concurrently.

Internally, the code maintains a vector of tasks that specify both the genome type and an approximate data size. For each task, it optionally simulates a heavy operation (representing actual data loading or transformation) and then constructs the appropriate index. The chosen index is returned as a boxed trait object, ensuring that both linear and graph-based indexes implement the same interface. By leveraging parallel iterators, the creation of multiple indexes can be distributed across CPU cores, helping manage computational load when building large or complex data structures.

#### Project Structure:

```plaintext
experiment_3_1_3/
└── Cargo.toml                     # Rust project configuration and dependencies
src/
├── main.rs                        # Main Rust script containing program logic
└── output.txt                     # Text output file
```

#### Cargo.toml

```toml
[package]
name = "experiment_3_1_2"
version = "0.1.0"
edition = "2021"

[dependencies]
rand = "0.9.0"
rayon = "1.10.0"
```

#### How to run:

run in powershell:

```powershell
cargo run | tee output.txt
```

(run main.rs and save the output in output.txt)
  

#### Explanation of the Output:
This program simulates the construction of genomic indexes for different types of genomic data using Rayon for parallel execution. Here’s a breakdown of the program and the output:

##### Key Concepts:

###### 1. GenomeType Enum:

* This enum represents different genomic data scenarios. It includes:
* Microbial: Small genomes (e.g., bacterial genomes, in the megabase range).
* Eukaryotic: Larger genomes (e.g., human genomes, in the gigabase range).
* Pangenome: A reference to multiple genomes (typically graph-based indexing).
* SingleCellTranscriptomic: Data from single-cell RNA sequencing (typically graph-based).
* HiCAssay: Data for 3D genome structures or contact maps (also graph-based).

###### 2. GenomicIndex Trait:

* This is a common interface that different index types implement. It has a describe() method to provide a textual description of the index.

###### 3. Index Types:

* LinearIndex: Used for smaller genomes (e.g., microbial or eukaryotic genomes).
* GraphIndex: Used for larger, more complex data such as pangenomes or single-cell transcriptomics.

###### 4. build_index Function:

* This function takes a GenomeType and an approximate genome size and returns an appropriate index (LinearIndex or GraphIndex).

###### 5. Parallelism with Rayon:

* The program leverages Rayon to parallelize the creation of indexes and their subsequent descriptions. This is particularly useful for handling large datasets or computationally expensive tasks.

###### 6. Simulating Heavy Operations:

* simulate_heavy_operation is a placeholder for operations that might be computationally expensive, such as I/O, compression, or complex graph-building tasks.

##### Output Analysis:
The program generates a list of tasks, where each task corresponds to a specific type of genomic data with an associated genome size. Then, it builds indexes in parallel for each task. The final result prints the descriptions of the generated indexes:

```rust
LinearIndex for genome of size: 2000000
LinearIndex for genome of size: 3000000000
GraphIndex with 12000000 nodes, colored = true
GraphIndex with 250000 nodes, colored = false
GraphIndex with 2000000 nodes, colored = false
```

* First Output (LinearIndex for microbial genome of size 2,000,000):
  * The first genome type is Microbial with a genome size of 2,000,000 bases (approximately 2MB), so it uses a LinearIndex.
* Second Output (LinearIndex for eukaryotic genome of size 3,000,000,000):
  * The second genome type is Eukaryotic with a genome size of 3,000,000,000 bases (approximately 3GB), so it also uses a LinearIndex.
* Third Output (GraphIndex for pangenome of size 12,000,000,000):
  * The third genome type is Pangenome with a genome size of 12,000,000,000 bases (approximately 12GB), which requires a more complex GraphIndex. The number of nodes is approximately 12,000,000, and it is colored (indicating multiple samples).
* Fourth Output (GraphIndex for single-cell transcriptomics of size 500,000,000):
  * The fourth genome type is SingleCellTranscriptomic with a genome size of 500,000,000 bases (approximately 500MB), requiring a non-colored GraphIndex. The number of nodes is approximately 250,000.
* Fifth Output (GraphIndex for Hi-C assay of size 1,000,000,000):
  * The fifth genome type is HiCAssay with a genome size of 1,000,000,000 bases (approximately 1GB), requiring a non-colored GraphIndex. The number of nodes is approximately 2,000,000.

#### Conclusion:
* The program efficiently parallelizes the construction of genomic indexes based on different genomic data types, using Rayon to speed up the process.
* LinearIndex is suitable for smaller genomes, while GraphIndex is used for larger, more complex genomic data.
* This approach scales well for different types of genomic data, enabling the handling of large datasets in parallel.
In practice, the program could be expanded with more sophisticated data structures (e.g., FM-index, suffix array) and additional functionalities for querying and manipulating these indexes.



