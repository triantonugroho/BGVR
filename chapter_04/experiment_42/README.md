## 4.2. Graph-Based Models for Gene Regulatory Networks (GRNs)

### experiment_42

#### 1. Nextflow 
The following scenario demonstrates how to compute a correlation-based adjacency matrix from a synthetic gene expression dataset using both Rust and Nextflow. Rust (with the rayon and ndarray crates) performs high-performance, parallelized computation of pairwise Pearson correlations, while Nextflow orchestrates the workflow in a reproducible and scalable manner. By combining these technologies, you can easily transition from local development to large-scale HPC or cloud environments without changing your core Rust code.

The main.nf script defines a single process, buildAdjacencyMatrix, which compiles and executes the Rust program. In Nextflow, processes are self-contained tasks that can be executed across different systems, from a local machine to cloud-based clusters. The script accepts parameters such as num_genes, num_samples, and an output_file name, which are passed to the Rust binary. Once the process finishes, Nextflow collects the result—here, a binary file representing the adjacency matrix—and handles logging and orchestration, thus enabling reproducible and portable workflows.

#### 2. Rust
The Rust code generates a synthetic gene expression matrix for a user-defined number of genes and samples, then computes pairwise Pearson correlations in parallel using the rayon crate. It stores data in an ndarray::Array2 for efficient indexing and uses an Arc<Mutex<>> to guard the shared adjacency matrix so that multiple threads can safely write the correlation values. By only calculating the upper triangle of the matrix and mirroring results, the program avoids redundant work. The final adjacency matrix is then written as binary data for space efficiency.

#### Files contents:
* experiment_42/
  * Cargo.toml (Cargo.toml file for dependencies)
* experiment_42/src/
  * main.rs (rust script)
  * main.nf (nextflow script)
  * ouput.txt (output file)
  * partial_adjacency.bin (output file)

#### How to run:

run in powershell:

```powershell
cargo run main.nf --num_genes 1000 --num_samples 50 --output_file partial_adjacency.bin
```

run in WSL:

```wsl
nextflow run main.nf --num_genes 1000 --num_samples 50
```

(run main.nf)

#### [dependencies]

```toml
bio = "2.0.3"
```

#### Explanation of the Output

##### 1. PWM Results (pwm_results.txt)
The Position Weight Matrix (PWM) provides the probabilities of observing each nucleotide (A, C, G, T) at specific positions within the sequence.

* Each row corresponds to a specific position in the sequence.
* The values represent the probability of finding a particular nucleotide at that position.
* For example:
  * Position 0: 𝐴 = 0.977, 𝐶 = 0.017, 𝐺 = 0.004, 𝑇 = 0.002
    → This indicates that at position 0, nucleotide A has a 97.7% chance of appearing.
  * Position 16: 𝐴 = 0.001, 𝐶 = 0.002, 𝐺 = 0.002, 𝑇 = 0.995
    → At position 16, nucleotide T is highly dominant with a 99.5% probability.

###### Interpretation:
* High probabilities for specific nucleotides at certain positions suggest sequence conservation or motif presence.
* Positions with evenly distributed probabilities suggest low specificity or variability at that site.

##### 2. MRF Results (mrf_results.txt)
The Markov Random Field (MRF) results represent the transition probabilities between nucleotides.

* Each row shows the probability of transitioning from one nucleotide to another.
* For example:
  * 𝐴 → 𝐴 = 0.3013 → Probability of A being followed by A is 30.13%
  * 𝐺 → 𝑇 = 0.2411 → Probability of G being followed by T is 24.11%

###### Interpretation:
* Higher transition probabilities between certain pairs suggest sequence patterns or dependencies.
* Lower transition probabilities imply that such transitions are rare in the sequence.

#### Conclusion
* The PWM and MRF outputs provide complementary insights into the sequence:
  * PWM shows the positional preferences of nucleotides, which helps identify conserved motifs.
  * MRF shows the transition patterns between nucleotides, revealing underlying sequence structure or biases.
* These results can be used to:
  * Identify binding motifs in DNA sequences.
  * Predict sequence patterns based on transition probabilities.
  * Improve sequence alignment or motif search algorithms.

