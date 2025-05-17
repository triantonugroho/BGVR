## 3.5. High-Performance Computing and Parallelization

### experiment_35

The following Rust code demonstrates a robust approach to counting k-mers in large genomic data files (FASTA or FASTQ) using concurrent, local aggregation. It leverages the Needletail crate for reading sequence data, Rayon for parallel processing, and DashMap for merging partial results. This design allows for efficient, multi-core counting of k-mers even in memory-constrained or large-scale scenarios, helping address the complexity of modern genomic datasets.

After parsing the input file through Needletail, the code collects all reads into a vector. It then spawns parallel tasks via Rayon, where each task builds a local HashMap of k-mer frequencies. This local accumulation greatly reduces concurrency overhead compared to updating a global structure for every k-mer. Once local aggregation finishes, each task merges its partial counts into a shared DashMap. Finally, the code extracts the global k-mer counts from the DashMap, sorts them, and serializes the result to JSON, producing a concise record of which k-mers appear and how frequently.

Each crate addresses different concerns. needletail uses efficient buffering for reading, avoiding the overhead of loading entire files at once. dashmap ensures lock-free updates to the shared map, while Arc allows multiple threads to share a single data structure with reference counting. In larger-scale HPC runs, ephemeral containers each run an instance of this program for a subset of input data, merging partial JSON outputs in a final step. Industrial setups often add extensive logging, memory usage checks, and specialized partitioning logic to guarantee consistent load balancing. By embedding such code in Nextflow or Snakemake pipelines, engineers achieve robust HPC workflows that adapt seamlessly from local clusters to large cloud-based HPC environments. We will discuss more about Nextflow in the following section.

#### Project Structure:

```plaintext
experiment_35/
├── Cargo.toml                      # Rust project configuration and dependencies
└── src/
    ├── main.rs                     # Main Rust script containing program logic
    ├── reads.fq.rar                # Compressed FASTQ reads file
    ├── kmer_counts.json.rar        # Compressed JSON file containing k-mer counts
    └── output.txt                  # Text output file
```

#### How to run:

run in powershell:

```powershell
cargo run | tee output.txt
```

(run main.rs and save the output in output.txt)
  
#### [dependencies]

```toml
rayon = "1.7"             # For parallel processing
needletail = "0.6.3"        # For FASTA/FASTQ parsing
dashmap = "6.1.0"           # For concurrent hashmap
serde = { version = "1.0", features = ["derive"] }  # For serialization
serde_json = "1.0"        # For JSON output
```

#### Explanation of the Output
This Rust program processes a FASTQ/FASTA file, extracts k-mers of length 𝑘, counts their occurrences in a parallelized manner, and then saves the counts to kmer_counts.json.

##### Step-by-Step Breakdown of Execution

###### 1. Command-line Arguments Handling
* The program reads the input sequence file (default: "reads.fq") and the k-mer length (default: 31).
* The k-mer length 𝑘 is parsed from the second command-line argument.

###### 2. Parallel K-mer Counting Using DashMap
* Uses DashMap (a thread-safe hashmap) to store k-mer → count mappings.
* Each thread builds a local k-mer count and then merges it into the global DashMap.
* This avoids locking issues that would arise from inserting every k-mer individually.

###### 3. Reading the Sequence File
* The program reads sequences from the FASTQ/FASTA file into a Vec<u8>.
* This approach loads the entire file into memory, which is efficient for moderate-sized datasets but may need streaming for very large datasets.

###### 4. Extracting and Counting k-mers
* Each read is processed in parallel using rayon::par_iter().
* k-mers are extracted from each sequence and stored in a local HashMap.
* After processing a read, the local k-mer counts are merged into the global DashMap.

###### 5. Sorting and Writing to JSON
* After parallel processing, k-mers are collected into a sorted vector.
* The results are serialized into JSON format and saved to kmer_counts.json.

##### Understanding the Output
###### 1. output.txt
```rust
Wrote k-mer counts to kmer_counts.json. Done!
```
* Indicates successful execution of the k-mer counting process.
* The final results were saved to kmer_counts.json.

###### 2. kmer_counts.json
Example snippet:
```json
{
  "k": 31,
  "counts": [
    [
      "AAAAAACCCGTTGCCGAAGGCAGCTTACTGC",
      1
    ],
    [
      "AAAAAACTCCAGTTTACACAAACTAGACTAC",
      1
    ],
    [
      "AAAAAACTCCGATTGCGAAGGCAGCCTGCTA",
      1
    ],
    [
      "AAAAAAGACAACCGTTTGCGAAGGCGCCTTC",
      1
    ],
    [
      "AAAAAAGGTGGAAGGGAAAACCAATGGCTCA",
      1
    ],
    ...
    [
      "TTTTTTTTGAGTAGTGCAGAGGTAGGCGGAA",
      2
    ],
    [
      "TTTTTTTTGATTAGTGCAGAGGTAGGCGGGA",
      1
    ]
  ]
}
```

###### Key Insights from kmer_counts.json
* "k": 31 → The program extracted 31-mers from the sequence data.
* "counts" contains (k-mer, count) pairs:
* Most k-mers appear once (count = 1), indicating they are unique.
* Some k-mers appear multiple times (e.g., "TTTTTTTTGAGTAGTGCAGAGGTAGGCGGAA" appears 2 times).
* The k-mers are sorted lexicographically, making it easier for downstream analysis.

#### Conclusion
##### 1.Efficient Parallel K-mer Counting

* Uses rayon to parallelize k-mer extraction.
* DashMap avoids contention issues while merging k-mer counts.

##### 2. Memory Considerations

* This approach loads all sequences into memory before processing.
* For large datasets, a streaming approach (e.g., processing chunks of reads) would be needed.

##### 3. Accurate Representation of k-mer Distribution

* The output JSON provides a complete frequency distribution of k-mers.
* Lexicographic sorting makes it easier to compare results across datasets.

##### 4. Applications

* This k-mer frequency data can be used for:
* Genome assembly
* Mutation detection
* Species identification
* Data compression in sequencing analysis

