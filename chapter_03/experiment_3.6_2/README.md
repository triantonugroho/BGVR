## 3.6. Putting It All Together—Rust and Nextflow Integration

### experiment_3.6_2

#### 1. Nextflow

Similarly, the alignment crate might accept FASTQ reads and the partial indexes, using a local or global alignment algorithm. HPC tasks parse multiple read files in parallel with rayon or distribute them across ephemeral containers. The summarizer merges alignment outputs by reading JSON files from multiple ephemeral tasks. If scaling to thousands of files or extremely large references, code must handle concurrency overhead and partial merges carefully.

This pipeline “scatters” read files across ephemeral containers in the ALIGN_READS step, then merges partial outputs in SUMMARIZE. Each ephemeral container runs one crate. Users can quickly scale HPC resources on local clusters or in the cloud, launching dozens or hundreds of ephemeral containers as needed. The modularity fosters continuous integration: each Rust crate is tested individually with cargo test, and the containers can be versioned in a Docker registry, guaranteeing reproducible runs.



#### Project Structure:

```plaintext
experiment_3.6_2/
└── Cargo.toml                      # Rust project configuration and dependencies
experiment_36_1/src/                # Note: Directory number mismatch (36_1 vs 36_2)
├── main.rs                         # Main Rust script containing program logic
├── main.nf                         # Nextflow workflow script
├── reference.fa                    # Reference FASTA file
├── ouput.txt                       # Output file (note: typo in filename, missing 't')
└── reads/
    └── sample.fq.rar               # Compressed FASTQ sample file
```

#### Cargo.toml

```toml
[package]
name = "experiment_3.6_2"
version = "0.1.0"
edition = "2021"

[dependencies]
serde = { version = "1", features = ["derive"] }
serde_json = "1"
bincode = "1.3.3"  # Use this version to avoid breaking changes
```

#### How to run:

run in powershell:

```powershell
cargo run main.nf | tee output.txt
```

(run main.nf and save the output in output.txt)
  

#### Explanation of the Output
The program consists of two main components:

##### 1. Nextflow Workflow (main.nf)

* It processes a genomic dataset using three main steps:
  * BUILD_INDEX: Generates a partial FM-index for a reference genome.
  * ALIGN_READS: Aligns read sequences to the reference using the built index.
  * SUMMARIZE: Collects all alignment results and summarizes them into a final output.
* Each process runs in a Docker container.

##### 2. Rust Program (main.rs)

* Runs the Nextflow workflow step by step.
* Handles serialization/deserialization using bincode.
* Calls nextflow using system commands and processes outputs.
* Outputs messages about each step and prints the Global Suffix Array at the end.

##### Understanding the Output
The output in output.txt is:

```rust
Global Suffix Array: [11, 4, 13, 6, 15, 8, 1, 12, 5, 14, 7, 0, 16, 10, 3, 9, 2]
```

This suggests that:

* A suffix array was computed for a given string or genomic sequence.
* The numbers represent positions in the original sequence, sorted lexicographically by suffix.
* The pipeline ran successfully, at least until this part of the Rust program.

#### Conclusion
* The Nextflow pipeline was triggered successfully from Rust.
* The suffix array computation was performed, likely as part of genome indexing.
* The program executed correctly, producing an ordered suffix array output.
* If no errors appeared, the indexing, alignment, and summarization processes completed as expected.
* The pipeline effectively integrates bioinformatics tools with Rust for automation and performance.

