## 5.1. Introduction to Sequence Analysis

### experiment_51_2

This pair of Rust programs demonstrates a practical approach to high-throughput sequence data analysis using concurrency and partial output merging. The first program processes a FASTQ file to compute basic read statistics—specifically the distribution of read lengths, along with totals for the number of reads and cumulative base count. It uses the needletail crate for efficient streaming of large (possibly gzipped) FASTQ data and relies on Rayon to allow for future parallel processing should the data parsing be adapted to multiple threads. The resulting statistics are stored in a JSON file, which can be considered a partial output if the dataset is chunked.

The program reads command-line arguments using clap and opens the specified FASTQ file with the needletail crate. It iterates over each sequence record in a streaming manner, extracting read lengths and aggregating them in a HashMap. The final counts are serialized as a PartialOutput struct, where the read length histogram, total read count, and total base count are stored. This output is written as JSON, and the program prints basic statistics to confirm successful processing.

#### Project Structure:

```plaintext
experiment_51_2/
├── Cargo.toml              # Rust project configuration and dependencies
└── src/
    ├── main.rs             # Main Rust script containing program logic
    ├── example.fastq.rar   # Compressed FASTQ example file
    └── partial_output.json # JSON partial output file for reference
```

#### How to run:

run in powershell:

```powershell
cargo run -- --input C:\Users\trian\BGVR\chapter_05\experiment_51_2\src\example.fastq
```

(run main.rs with input path of fastq file)
  
#### [dependencies]

```toml
rayon = "1.7"
needletail = "0.6.3"
serde = { version = "1.0", features = ["derive"] }
serde_json = "1.0"
anyhow = "1.0"
clap = { version = "4.4", features = ["derive"] }
```

#### Explanation of the Output

Command Execution

```rust
cargo run -- --input C:\Users\trian\BGVR\chapter_05\experiment_51_2\src\example.fastq
```

* This runs the Rust program with an input FASTQ file containing sequencing reads.

* The program processes the file and outputs statistics about read lengths.

* Results are saved in a JSON file (partial_output.json).

Output in partial_output.json

```json
{
  "histogram": {
    "counts": {
      "224": 49141,
      "223": 2133,
      "222": 25,
      "225": 3593
    }
  },
  "total_reads": 54892,
  "total_bases": 12297218
}
```

##### 1. Histogram (Read Length Distribution)

```json
"histogram": {
  "counts": {
    "224": 49141,
    "223": 2133,
    "222": 25,
    "225": 3593
  }
}
```

* The histogram tracks the number of reads for each length:

  * 224 bases → 49,141 reads (most common read length)

  * 223 bases → 2,133 reads

  * 222 bases → 25 reads (very few)

  * 225 bases → 3,593 reads

* The majority of reads are 224 bases long, suggesting a consistent sequencing length.

##### 2. Total Reads Processed

```rust
"total_reads": 54892
```

* The FASTQ file contains 54,892 sequencing reads.

##### 3. Total Bases Counted

```rust
"total_bases": 12297218
```

* The total number of nucleotides across all reads is 12,297,218 bases.

* This is calculated as:

∑(read length×count of reads)

Example calculation:

(224 × 49141) + (223 × 2133) + (222 × 25) + (225 × 3593) = 12,297,218

Printed Output in Terminal

```rust
Wrote partial output to partial_output.json
Reads: 54892, Bases: 12297218
```

* Confirms successful processing and file writing.

#### Conclusion

1. Efficient Read Processing:

   * The program efficiently parses the FASTQ file line by line (streaming).

   * Uses hashmaps to store read length distributions.

2. Uniform Read Lengths:

   * Majority of reads are 224 bases, indicating a standardized sequencing process.

   * Minor variations (223, 222, 225 bases) could be due to sequencing errors or adapter trimming.

3. Scalability:

   * Handles large FASTQ files and outputs summary statistics efficiently.

   * Can be used for quality control in sequencing experiments.
