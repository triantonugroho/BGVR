## 5.6. Structural Variant Detection

### experiment_5_6

This Rust code provides a parallelized, chunk-based approach to naive split-read detection for structural variant (SV) analysis. It reads alignment segments—represented by minimal fields such as read ID, chromosome, alignment start, and orientation—from a streaming JSON input. By processing chunks of alignment records, it avoids storing the entire dataset in memory at once, which is crucial for high-performance computing (HPC) environments. The partial results for each chunk are written out so that, in the event of a disruption, the pipeline can resume without reprocessing the entire input file.

After these partial results have been produced, the code merges them into a final list of breakpoints, again leveraging Rust’s ownership and concurrency primitives to ensure safe parallel usage. A naive method is used to detect breakpoints: if a single read has multiple alignments that either span widely on the same chromosome or map to distinct chromosomes, a structural event is presumed. Additional logic, such as orientation checks, read-pair data, and coverage-based filters, can be integrated by extending the breakpoint detection logic.

The program uses the clap crate to handle command-line arguments, allowing the user to specify the path to an alignment file in JSON, a chunk size, an output directory for partial breakpoint results, and a final merged output file. Alignments are loaded in batches of size specified by --chunk-size, grouped by read ID, and then analyzed in parallel with Rayon. For each group of segments tied to the same read, the code sorts the segments by their starting coordinate and calls detect_breakpoints, which checks for large gaps or cross-chromosomal mappings.

Once all chunks are processed, the program scans the partial output directory for any files named partial_breakpoints_*.json, reading each and merging them into a single vector of breakpoints. This final, aggregated list is written to the user-specified output file. In HPC or industrial settings, ephemeral containers can each handle a subset of the data, generating multiple sets of partial outputs that are subsequently merged. Rust’s safety guarantees around shared data structures simplify parallel expansions, and advanced crates like polars, ndarray, or linfa can be introduced to incorporate statistical modeling or machine learning for more precise variant detection.

#### Project Structure:

```plaintext
experiment_5_6/
├── Cargo.toml                                  # Rust project configuration and dependencies
└── src/
    ├── main.rs                                 # Main Rust script containing program logic
    ├── alignment_data.json                     # Alignment data JSON input file
    ├── merged_breakpoints.json                 # Merged breakpoints JSON output file
    ├── output.txt                              # Text file output
    └── partial_breakpoints/
        ├── partial_breakpoints_0.json          # Partial breakpoints in chunk 0 JSON output file
        └── partial_breakpoints_1.json          # Partial breakpoints in chunk 1 JSON output file
```

#### Cargo.toml

```toml
[package]
name = "sv_detector"
version = "0.1.0"
edition = "2024"

[dependencies]
anyhow = "1.0"
rayon = "1.7"
serde = { version = "1.0", features = ["derive"] }
serde_json = "1.0"
clap = { version = "4.3", features = ["derive"] }
```

#### How to run:

run in powershell:

```powershell
cargo run --release -- --alignment-input alignment_data.json --chunk-size 2 --partial-output-dir partial_breakpoints --merged-output merged_breakpoints.json | tee output.txt
```

(run main.rs with chunk size 2, input file name alignment_data.json, output directory partial_breakspoints and output file name merged_breakpoints.json and save the output text in output.txt) 
  

#### Explanation of the Output

##### 1. Chunk Processing

* The program reads alignment_data.json in chunks of size 2 (as specified by --chunk-size 2).

* The dataset contains 3 records, so it is split into two chunks:

  * Chunk 0: First two records (both with read_id="read_1").

  * Chunk 1: Last record (read_id="read_2").

##### 2. Breakpoint Detection (Chunk 0)

* The program groups alignment records by read_id. In Chunk 0, only read_id="read_1" is present.

* The two alignment segments for read_id="read_1":

```json
{"read_id": "read_1", "chrom": "chr1", "start": 100, "cigar": "50M", "orientation": "+"}
{"read_id": "read_1", "chrom": "chr1", "start": 200, "cigar": "50M", "orientation": "-"}
```

* Since both alignments belong to the same chromosome (chr1), a naive intra-chromosomal breakpoint is detected at:

```rust
breakpos = start (100) + parse_cigar_len("50M") = 100 + 50 = 150
```

This breakpoint is saved in partial_breakpoints_0.json:

```rust
{"breakpoints":[{"read_id":"read_1","chrom":"chr1","pos":150,"sv_type":"intra-chr"}]}
```

##### 3. Breakpoint Detection (Chunk 1)

* This chunk contains only read_id="read_2" with a single alignment segment:

```json
{"read_id": "read_2", "chrom": "chr2", "start": 150, "cigar": "30M", "orientation": "+"}
```

* Since breakpoints require at least two segments per read, no breakpoint is detected.

* partial_breakpoints_1.json is empty:

```json
{"breakpoints":[]}
```

##### 4. Merging Partial Results

* The program merges the breakpoints from all partial_breakpoints_*.json files.

* Since partial_breakpoints_1.json contains no breakpoints, the final merged result is the same as partial_breakpoints_0.json:

```json
[{"read_id":"read_1","chrom":"chr1","pos":150,"sv_type":"intra-chr"}]
```

#### Conclusion

* The program correctly detects intra-chromosomal breakpoints based on consecutive alignments from the same read.

* No breakpoints are detected for read_id="read_2" because it has only one alignment segment.

* The chunking mechanism is working as expected, processing the dataset in batches.

* The merging step correctly consolidates detected breakpoints into merged_breakpoints.json.

