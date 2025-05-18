## 4.2. Graph-Based Models for Gene Regulatory Networks (GRNs)

### experiment_42

#### 1. Nextflow 
The following scenario demonstrates how to compute a correlation-based adjacency matrix from a synthetic gene expression dataset using both Rust and Nextflow. Rust (with the rayon and ndarray crates) performs high-performance, parallelized computation of pairwise Pearson correlations, while Nextflow orchestrates the workflow in a reproducible and scalable manner. By combining these technologies, you can easily transition from local development to large-scale HPC or cloud environments without changing your core Rust code.

The main.nf script defines a single process, buildAdjacencyMatrix, which compiles and executes the Rust program. In Nextflow, processes are self-contained tasks that can be executed across different systems, from a local machine to cloud-based clusters. The script accepts parameters such as num_genes, num_samples, and an output_file name, which are passed to the Rust binary. Once the process finishes, Nextflow collects the result—here, a binary file representing the adjacency matrix—and handles logging and orchestration, thus enabling reproducible and portable workflows.

#### 2. Rust
The Rust code generates a synthetic gene expression matrix for a user-defined number of genes and samples, then computes pairwise Pearson correlations in parallel using the rayon crate. It stores data in an ndarray::Array2 for efficient indexing and uses an Arc<Mutex<>> to guard the shared adjacency matrix so that multiple threads can safely write the correlation values. By only calculating the upper triangle of the matrix and mirroring results, the program avoids redundant work. The final adjacency matrix is then written as binary data for space efficiency.

#### Project Structure:

```plaintext
experiment_42/
├── Cargo.toml                     # Rust project configuration and dependencies
└── src/
    ├── main.rs                    # Main Rust script containing program logic
    ├── main.nf                    # Nextflow workflow script
    ├── ouput.txt                  # Output file
    └── partial_adjacency.bin      # Partial adjacency output file in binary format
```

#### How to run:

run main.rs in powershell:

```powershell
cargo run main.rs  --num_genes 1000 --num_samples 50 --output_file partial_adjacency.bin ! tee output.txt
```

run main.nf in powershell:

```powershell
cargo run main.nf  --num_genes 1000 --num_samples 50 --output_file partial_adjacency.bin
```

run main.nf in WSL:

```wsl
nextflow run main.nf --num_genes 1000 --num_samples 50 --output_file partial_adjacency.bin
```

#### [dependencies]

```toml
bio = "2.0.3"
```

#### Explanation of the Output
The output of the Nextflow and Rust program confirms that the execution proceeded correctly and the adjacency matrix was successfully generated and saved. Let's break down the key output lines:
##### 1. Number of genes: 1000
* This line confirms that the --num-genes parameter was correctly parsed and passed to the Rust binary.
* num_genes = 1000 means that the generated synthetic dataset contains 1000 genes (or rows) in the matrix.

##### 2. Number of samples: 50
* This line confirms that the --num-samples parameter was correctly passed.
* num_samples = 50 means that each gene's expression is sampled over 50 samples (or columns).

##### 3. Output file: partial_adjacency.bin
* This line confirms that the --output parameter was received and processed correctly.
* The adjacency matrix was written to partial_adjacency.bin in the current working directory of Nextflow.
* The file path is relative to the execution environment defined by Nextflow.

##### 4. Correlation adjacency matrix written to partial_adjacency.bin
* This line indicates that the matrix computation was completed successfully and saved to the specified binary file.
* The adjacency matrix represents the pairwise Pearson correlation between genes based on their expression profiles.
* The matrix is symmetric because the correlation between gene 𝑖 and gene 𝑗 is the same as between gene  and gene 𝑖
  
##### Structure of partial_adjacency.bin
* The file partial_adjacency.bin contains the adjacency matrix in binary format.
* The matrix is stored as a flat sequence of 64-bit floating-point values (f64) row by row.
* Total size of the file can be calculated as:
  size = 1000 × 1000 × 8 bytes = 8,000,000 bytes = 8 MB

#### Conclusion
* The Nextflow script executed successfully without any missing output or path errors.
* The Rust binary correctly handled the CLI parameters and computed the adjacency matrix using parallel processing with rayon.
* Pearson correlation calculation was handled efficiently in parallel, storing results in a symmetric matrix.
* The binary file was created as expected and contains the full correlation matrix.
