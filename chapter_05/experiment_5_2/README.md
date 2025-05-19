## 5.2. Preprocessing and Quality Control

### experiment_5_2

The Rust program below applies a sliding window trimming algorithm to reads from a FASTQ file, chunk by chunk, and writes partial results to disk. By processing only a limited subset of reads at a time (chunk_size) and using parallel iteration through Rayon’s data parallelism, it addresses memory constraints that can arise with extremely large datasets. In addition, the partial outputs are serialized in JSON, creating a natural checkpointing mechanism that is particularly useful when dealing with large volumes of sequence data in a production or HPC environment.

The sliding window trimming approach itself scans each read’s quality scores from both ends. It uses a configurable window size and average Phred score threshold to determine how far inwards to trim low-quality bases. This design ensures that reads with persistently low-quality regions at the start or end are trimmed back to a region of higher confidence. By decoupling the trimming logic from the parallel iteration, one can easily substitute a more sophisticated algorithm or integrate advanced crates such as linfa for machine learning–based trimming heuristics if desired.

After parsing command-line arguments using the clap crate, the program uses needletail to open and read the FASTQ file, which can be either plain or gzipped. Instead of loading all records at once, it takes in batches of size chunk_size. Each batch is processed via Rayon’s par_iter(), thereby distributing the trimming workload across available CPU cores.

The trimming logic is encapsulated in the sliding_window_trim function, which uses a rolling sum of Phred scores to compute the average quality in a window. If this average drops below the specified threshold, the function updates the boundaries that define the start or end of the trimmed region. The result is a pair of sliced vectors: the retained sequence bases and the corresponding quality scores.

Once the trimming is complete for each batch, the program writes a JSON file containing all trimmed reads for that chunk. This partial output strategy is highly advantageous in large-scale operations. It lets the pipeline resume from the last completed chunk if the process is interrupted, and it mitigates memory spikes by limiting how many records are held in memory simultaneously. Finally, the program reports basic statistics, including how many partial files were written. A subsequent step in the pipeline could merge or further analyze these JSON files, thereby enabling a flexible, HPC-friendly workflow.

#### Project Structure:

```plaintext
experiment_5_2/
├── Cargo.toml              # Rust project configuration and dependencies
└── src/
    ├── main.rs             # Main Rust script containing program logic
    ├── output.txt.rar      # Compressed output.txt file
    └── example.fastq.rar   # Compressed FASTQ example file
```

#### Cargo.toml

```toml
[package]
name = "experiment_5_2"
version = "0.1.0"
edition = "2024"

[dependencies]
rayon = "1.7"
needletail = "0.6.3"
serde = { version = "1.0", features = ["derive"] }
serde_json = "1.0"
anyhow = "1.0"
clap = { version = "4.4", features = ["derive"] }
```

#### How to run:

run in powershell:

```powershell
cargo run -- --input C:\Users\trian\BGVR\chapter_05\experiment_52\src\example.fastq
```

(run main.rs using input data path)
  

#### Explanation of the Output
The output output.txt shows a series of FASTQ reads that have been processed by the sliding window trimming algorithm implemented in main.rs. Each read consists of:

##### 1. Read ID

* Example: SRR11192680.1 1 length=224

* This represents a unique identifier for a sequencing read along with its length.

##### 2. Sequence

* Example:

```rust
ACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGCAGACGGGAAATTAAGTCAGTTGTGAAAGTTTGCGGCTCAACCGTAAAATTGCAGTTGATACTGGTTTCCTTGAGTGCAGTTGAGGCAGGCGGAATTCGTGGTGTAGCGGTGAAATGCTTAGATATCACGAAGAACCCCGATTGCGAAGGCAGCTTGCTAAACTGTATCTGACGCTC
```

* This represents the DNA sequence of the read.

##### 3. Quality Scores (Phred Score)

* Example:

```rust
GGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGGGGFFFGGGGGGGGGGGGGGGGGGDGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGGGGGGGGGGFGGGGGGGGFGGGGGGGGGGEGGGEGEEGGGGGGGGGGFGGGGGGGFFGGGFGFF;:FGGGGC=BEGGGGGGGGGCFGGCGE>@FGDCF?<DBEFFFGCF7<5<A??E?E
```

* These scores indicate the base call confidence from the sequencing machine.

* High-quality bases are represented by higher ASCII characters, while lower-quality bases have lower ASCII values.

Each read is separated by a dashed line (--------------------------------------) to differentiate between multiple reads.

#### Conclusion
* The script successfully parses and processes FASTQ sequences from the input file.

* The sliding window trimming algorithm runs correctly, evaluating the quality scores of each read and trimming low-quality bases.

* The trimmed sequences and corresponding quality scores are saved into JSON files (partial_trim_output_chunk_X.json) for further analysis.

* The script utilizes parallel processing (rayon::par_iter()) to speed up the trimming process.

* The chunk-based approach ensures that large files are processed efficiently without excessive memory usage.


