## 5.3. Alignment and Mapping Algorithms

### experiment_5_3

The following code illustrates how to build partial FM-indexes for a reference genome and then perform parallel alignment of reads in Rust. The reference genome is subdivided into multiple chunks, each of which is indexed separately. Splitting the genome in this way can help manage memory usage when dealing with very large references in HPC environments. The partial indexes are each stored on disk, allowing ephemeral tasks to construct or load only the chunks they need for specific stages of an alignment pipeline.

After building each chunk’s FM-index, the program reads a list of query sequences (one per line in a simple text file) and aligns them against the partial indexes. Rust’s ownership and concurrency model ensures safe parallel processing of all read queries with Rayon, distributing the workload among available CPU cores. The alignment step collects all matching positions from each partial index, adjusts them with an offset, and compiles these into a single list of “global” positions relative to the original reference.

The command-line interface, powered by the clap crate, provides user-configurable parameters for chunk size, FM-index sampling rate, and input/output paths. Once the reference genome is loaded into memory, it is sliced into segments of size defined by --reference-chunk-size. Each segment is used to construct an FM-index, which is then serialized and stored under the output directory.

Alignment is carried out by first loading all the partial index files from disk. Each read is processed with Rayon’s .par_iter(), distributing alignments across threads. For each partial index, the code deserializes the stored FM-index and locates occurrences of the read in that particular segment. The positions reported by the FM-index are relative to the segment being searched, so the code adds a chunk offset to map these hits to the overall reference coordinates.

In HPC or industrial settings, the partial indexes could be built on different nodes, stored on a shared filesystem, and loaded by many ephemeral containers or tasks. One might distribute the read sets in similar chunks, then collect final alignments in a subsequent merge step. Rust’s thread-safe concurrency primitives and data-parallel abstractions simplify scaling out to very large input data, while libraries like ndarray, polars, or linfa can be added to support advanced analytics on alignment results. This approach forms a production-ready foundation for building a high-throughput, distributed alignment pipeline in Rust.

#### Project Structure:

```plaintext
experiment_5_3/
├── Cargo.toml             # Rust project configuration and dependencies
└── src/
    ├── main.rs            # Main Rust script containing program logic
    ├── reads.fq           # FASTQ file containing sequencing reads
    └── reference.fa       # FASTA file containing reference sequences
```

#### Cargo.toml

```toml
[package]
name = "experiment_5_3"
version = "0.1.0"
edition = "2024"

[dependencies]
anyhow = "1.0"
rayon = "1.8"
bio = "2.1.0"
serde = { version = "1.0", features = ["derive"] }
serde_json = "1.0"
bincode = "2.0.1"
clap = { version = "4.0", features = ["derive"] }
```

#### How to run:

run in powershell:

```powershell
cargo run -- --reference reference.fa --reads reads.fq
```

(run main.rs using two input dataset path)
  

#### Explanation of the Output

This Rust program processes a reference genome and a set of reads to perform FM-index-based alignment, meaning it searches for occurrences of each read within the reference genome. The expected output is a JSON file containing the alignment results.

##### 1. Preprocessing Steps

* Reads the reference genome from a file.

* Converts all characters to uppercase and removes invalid characters.

* Appends a sentinel character ($) if missing.

* Constructs the FM-index (Suffix Array + BWT).

##### 2. Processing Reads

* Reads the reads file (each line is treated as a read sequence).

* Uses parallel processing (rayon::par_iter()) to align each read using the FM-index.

* Outputs alignment results as a JSON file.

##### 3. Output
During execution, the program prints various status messages to the console to provide insight into the process:

```rust
✅ Processed reference genome length: 120000
✅ Loaded 50 reads
🔍 Searching for read: "ACTGTC"
✅ Interval found: Interval { lower: 500, upper: 505 }
✅ Found 5 matches for read
🔍 Searching for read: "GGGTTT"
❌ No match found for read.
✅ Wrote alignment results for 50 reads to "alignments.json"
```

##### 4. alignments.json 
The program writes a JSON file containing the read sequences and their match positions in the reference genome. If a read is found in the reference genome, it lists all matching positions; otherwise, it provides an empty list.

```json
[
    {
        "read_seq": "ACTGTC",
        "matches": [500, 501, 502, 503, 504]
    },
    {
        "read_seq": "GGGTTT",
        "matches": []
    }
]
```

#### Conclusion
This program efficiently aligns a set of short DNA sequences (reads) against a reference genome using the FM-index, enabling quick and memory-efficient searching. The parallel processing of reads improves performance significantly. The JSON output provides an easy-to-parse format for downstream analysis, such as variant calling or genome assembly.
