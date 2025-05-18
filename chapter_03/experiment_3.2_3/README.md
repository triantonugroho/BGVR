## 3.2. Sequence Data Structures and Strings Algorithms

### experiment_32_3

#### 1. Nextflow
The Nextflow script below provides two processes: one for building the Rust binary (optional if you prefer a precompiled artifact), and another to run the partial suffix array program on a given FASTA file. Each chunk in the FASTA sequence is processed in parallel using Rayon, yielding multiple partial arrays that are then serialized as JSON. For industrial-scale tasks, ephemeral nodes might each handle a single chunk of the sequence before merging partial JSON outputs into a single suffix array. The Rust code itself exemplifies a simplified approach, illustrating how chunking, parallel sorting, and minimal offset adjustments can be done in a safe, robust manner.

#### 2. Rust
In the Rust code, we first read the FASTA file and remove any header lines that start with “>,” concatenating the remaining lines into a single String. Next, we split this sequence into manageable slices by selecting a chunk size (in this example, one million bases) and calling as_bytes().chunks(chunk_size) on the sequence. For each resulting slice, we construct a naive suffix array by enumerating all possible substring start positions and sorting them; although this method is not practical for extremely large data, it showcases the core idea. By leveraging Rayon’s .par_iter(), each slice is processed on a separate CPU core for improved throughput. Finally, the code serializes each partial suffix array as JSON and writes it to partial_suffix_arrays.json, which can then be merged or further analyzed in subsequent pipeline steps.

To run the code locally, you then either create a new Rust project or clone an existing one containing the partial suffix array code, and place main.nf (the Nextflow script) in that same directory. Finally, launch the pipeline using nextflow run main.nf --fasta /path/to/large_sequence.fa, which triggers the compile process (compiling the Rust program) followed by analysis (executing the compiled binary on your FASTA file). The result is a file named partial_suffix_arrays.json containing serialized partial suffix arrays.

In a large-scale environment, each ephemeral job might generate its own partial JSON file, which can subsequently be merged by a follow-up Nextflow process or separate program—thereby combining all suffix arrays, adjusting offsets for each chunk, and sorting them into a final global index. If your workflow calls for boundary-spanning operations, such as building full BWT or handling k-mer overlaps, you can include additional bases (e.g., k-1 overlap) at each chunk boundary. Because the naive approach in this demonstration sorts all suffix positions in memory, it can be computationally and memory-intensive for extensive data; production-ready solutions often rely on advanced suffix array algorithms (like SA-IS), chunk partitioning, or fully distributed indexing frameworks. Furthermore, one would typically include more robust error handling, file validation, and resource management to handle extremely large inputs. Overall, Nextflow manages the high-level workflow in stand-alone mode—whether on a local PC or a cloud/HPC environment—while the Rust code handles partial suffix array creation in a thread-safe, concurrent manner.

#### Project Structure:

```plaintext
experiment_32_3/
├── Cargo.toml                         # Rust project configuration and dependencies
└── src/
    ├── main.rs                        # Main Rust script containing program logic
    ├── main.nf                        # Nextflow workflow script
    ├── large_sequence.fa              # FASTA file containing sequence data
    ├── partial_suffix_arrays.json     # JSON output file containing partial suffix arrays
    └── output.txt                     # Text output file
```

#### How to run:

run in powershell:

```powershell
cargo run main.nf | tee output.txt
```

(run main.nf and save the output in output.txt)
  
#### [dependencies]

```toml
rayon = "1.7"
serde = { version = "1", features = ["derive"] }
serde_json = "1"
```

#### Explanation of the Output:
The program you are running is a Nextflow pipeline that automates the compilation and execution of a Rust program for building partial suffix arrays on a large DNA sequence stored in a FASTA file. Here’s a breakdown of each part of the output:

##### 1. output.txt:
This file contains a message that reports the number of partial suffix arrays generated and saved to the partial_suffix_arrays.json file. The output you received is:

```rust
Generated 1 partial arrays in partial_suffix_arrays.json
```

* This indicates that the Rust program successfully processed the FASTA sequence and built 1 partial suffix array for the sequence in large_sequence.fa.
* This implies that the FASTA file was processed in 1 chunk, likely because the chunk size is set to 1,000,000 characters, and the input sequence fits within this size.

##### 2. partial_arrays_suffix.json:
This JSON file contains the serialized representation of the partial suffix array that was computed by the Rust program. The array is saved in JSON format for potential merging and further analysis.

Here’s a snippet from the file:

```json
[
  {
    "start_pos": 0,
    "suffix_positions": [
      589,
      979,
      692,
      1126,
      1045,
      365,
      718,
      54,
      104,
      1157,
      633,
      1144,
      1075,
      658,
      596,
      986,
      804,
      ....
    ]
  }
]
```

* start_pos: This is the starting position of the chunk in the original sequence (0 in this case, because the entire sequence is contained in this chunk).
* suffix_positions: This array represents the sorted indices of suffixes for the given chunk of the sequence. These are the starting positions of the suffixes in lexicographical order. The numbers represent positions in the sequence where suffixes begin. For example:
  * 589: This means the suffix starting at position 589 in the sequence is the lexicographically smallest among all suffixes in this chunk.
  * The list is sorted, so each subsequent value represents the next lexicographically smallest suffix.
The full list of suffix positions would contain many entries depending on the length of the sequence and chunk size.

##### 3. main.rs Execution:
* FASTA File Processing: The FASTA file is read, skipping the header lines, and the sequence is extracted.
* Chunking: The sequence is partitioned into chunks of 1,000,000 characters. The program then processes each chunk in parallel to generate a partial suffix array.
* Parallel Suffix Array Generation: Using the rayon library, each chunk is processed in parallel. For each chunk:
* A list of suffix starting positions is created.
* The suffixes are sorted lexicographically based on the substrings starting at each position.

##### 4. Nextflow Pipeline:
* The main.nf Nextflow script automates the process of compiling and running the Rust code. The pipeline first compiles the Rust program (compile process) and then runs it with the specified FASTA file (analysis process).
* The output from the Rust program (partial suffix arrays) is serialized into a JSON file and made available for further processing.

#### Conclusion:

##### 1. Successful Execution:

* The program successfully ran, processed the FASTA file, and generated 1 partial suffix array as expected. This is confirmed by the output message in output.txt and the corresponding partial_suffix_arrays.json file.

##### 2. Efficiency of Parallelization:

* The parallelization via the rayon library allowed for processing the chunks of the sequence in parallel. Even though you only have one chunk in this run, the approach would scale well for larger datasets by splitting the sequence into more chunks and processing them in parallel.

##### 3. Suffix Array Generation:

* The suffix positions stored in partial_arrays_suffix.json are lexicographically sorted indices of suffixes for a specific chunk of the sequence. These will later be used for merging into a full suffix array if needed.

##### 4. Nextflow Workflow:

* The Nextflow pipeline successfully automated the compilation and execution of the Rust code, demonstrating the power of Nextflow for managing complex computational workflows.

##### 5. Scalability:

* The current setup is efficient for handling large FASTA files, especially when the sequence is large enough to be divided into multiple chunks. As the chunk size and sequence size increase, the parallelization will become more valuable.
* To conclude, this output indicates the system successfully implemented a parallel suffix array generation approach using Rust and Nextflow, and it provides a foundation for scaling this method to handle much larger biological sequence datasets.
