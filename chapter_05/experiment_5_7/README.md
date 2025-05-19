## 5.7. RNA-seq and Transcriptomic Analysis

### experiment_57

This Rust program uses rust-htslib to read a BAM file in chunks, tallying up transcript-level expression counts in parallel using Rayon. By chunking the records, you can handle large BAM files in a memory-efficient way. Each chunk’s partial results are written to disk, enabling a strategy where ephemeral jobs (e.g., in an HPC cluster) each process a subset of the data.

When the partial tasks are completed, the code merges all partial results by reading them back and summing the counts for each transcript ID. The final result is then serialized into JSON. This aligns with common HPC workflows, in which data parallelism is used to accelerate computations, and partial outputs mitigate failures by letting pipelines resume from the last successfully processed step.



#### Project Structure:

```plaintext
experiment_57/
├── Cargo.toml                                  # Rust project configuration and dependencies
└── src/
    ├── main.rs                                 # Main Rust script containing program logic
    ├── example.bam                             # BAM input file
    ├── example.gtf                             # GTF input file
    ├── merged_counts.json                      # Merged counts JSON output file
    ├── output.txt                              # Text file output
    └── partial_counts/
        ├── partial_counts_chunk_1.json         # Partial counts in chunk 1 JSON output file
        └── partial_counts_chunk_2.json         # Partial counts in chunk 2 JSON output file
```

#### How to run:

run in wsl:

```wsl
cargo run -- --bam-input example.bam --annotation example.gtf --chunk-size 5000 --partial-outdir partial_counts --merged-output merged_counts.json | tee output.txt
```

(run main.rs with chunk size 5000, input file name example.bam, output directory partial_counts and output file name merged_counts.json and save the output text in output.txt) 
  
#### [dependencies]

```toml
anyhow = "1.0"
rust-htslib = "0.49.0"
rayon = "1.7"
serde = { version = "1.0", features = ["derive"] }
serde_json = "1.0"
clap = { version = "4.3", features = ["derive"] }
```

#### Explanation of the Output and Conclusion
The program processes a BAM file containing aligned sequencing reads and maps them to transcript annotations from a GTF file. The primary output consists of transcript-level read counts, which are then merged into a final JSON file containing the total counts for each transcript.

##### 1. Transcript Mapping and Counting

* The program reads BAM records and determines their chromosomal locations.

* It cross-references these locations with transcript annotations from the GTF file.

* Each read that overlaps a transcript is counted and assigned to the corresponding transcript ID.

* Parallel processing (using Rayon) speeds up the counting process by distributing read-processing tasks across multiple CPU threads.

##### 2. Partial Counts Storage and Merging

* The counts are first stored as intermediate JSON files (partial_counts_chunk_X.json) in a specified directory.

* These partial counts are later merged into a final merged_counts.json file, summing up the counts from all partial files.

* The final JSON output contains a list of transcript IDs and their respective read counts.

##### 3. Printed Debugging Information

* The program prints the regions being processed, displaying chromosome names and genomic positions.

* If no transcripts overlap a region, it prints a message indicating no match was found.

* It also outputs a confirmation message when merging is completed, along with the number of transcripts counted.

#### Conclusion

* The script successfully extracts read counts per transcript from a BAM file using annotation data from a GTF file.

* Parallel processing significantly enhances efficiency, allowing for faster counting of reads.

* The final merged output provides a structured JSON representation of transcript counts, which can be used for downstream analyses such as differential expression analysis or quality control.

* The approach ensures modularity by storing intermediate results, making it easier to debug and verify results at different stages of execution.

