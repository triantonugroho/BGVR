## 2.3. Data Structures and Algorithms in Rust

### experiment_23_2

Real-world genomics pipelines must orchestrate string-based (alignment, assembly) and numerical (factorization, clustering) algorithms at large scale. Rust’s type safety, along with HPC frameworks like Nextflow, helps ensure that tasks can be distributed across multiple nodes without risk of memory corruption or data races (Holmes 2022). Below is a mini-example combining concurrency, the bio crate for reading FASTA, and partial de Bruijn graph assembly.

#### Files contents:
* main.rs (rust script)
* reads.fasta (fasta file)
* Cargo.toml (Cargo.toml file)
* output.txt (output file)

#### How to run:

run in powershell:

```powershell
cargo run | tee output.txt
```

(run main.rs and save the output in output.txt)
  
#### [dependencies]

```toml
rayon = "1.7"
bio = "2.2.0"
anyhow = "1.0"
structopt = "0.3"
log = "0.4"
env_logger = "0.11.6"
```

#### Explanation of the Output
This Rust program constructs a de Bruijn graph from a FASTA file containing DNA reads. The output provides details about the graph construction process, including the number of reads processed and the final number of (k-1)-mer prefixes.

##### Step-by-Step Execution

###### 1. Logging and Setup

* The program initializes a logger to record output both to the console and a log file (output_new.txt).
* It sets k = 21, which means that it will analyze 21-mers (subsequences of length 21).

```rust
Building de Bruijn graph with k-mer size = 21
```

* This confirms that the program is using k = 21 to construct the de Bruijn graph.

###### 2. Reading the FASTA File

* The program opens a FASTA file (reads.fasta) and processes each sequence.
* The length of each DNA sequence is logged:

```sh
Read record of length: 59
Read record of length: 61
Read record of length: 42
Read record of length: 46
Read record of length: 48
```

* This means the FASTA file contains 5 DNA sequences of varying lengths.

###### 3. Processing Reads into a De Bruijn Graph

* The program extracts overlapping k-mers and constructs a de Bruijn graph.
* Reads are processed in parallel (using Rayon) for efficiency.
* Each (k-1)-mer prefix is stored as a node, and suffixes are stored as edges.

```sh
Total reads found: 5
```

* Confirms that 5 reads were successfully processed.

```sh
Built graph with 154 (k-1)-mer prefixes
```

* The final de Bruijn graph contains 154 unique (k-1)-mers, meaning the program found 154 overlapping sequences of length k-1 = 20.
  
#### Conclusion
* The program successfully builds a de Bruijn graph from DNA reads using k = 21.
* Parallel processing (via rayon::par_chunks) improves efficiency when handling large datasets.
* The final graph contains 154 unique prefixes, indicating a dense connectivity in the DNA sequences.
* This method is useful in genome assembly, where overlapping sequences help reconstruct longer DNA fragments.

Overall, the program demonstrates efficient de Bruijn graph construction from FASTA sequences, making it valuable for bioinformatics applications.

