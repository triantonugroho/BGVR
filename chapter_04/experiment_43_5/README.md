## 4.3. Motif Discovery and Regulatory Element Identification

### experiment_43_5

#### 1. Nextflow 
Below is an illustrative Nextflow pipeline and Rust code combination that implements motif-discovery workflow. The pipeline demonstrates how to chunk genomic data, scan for motifs in parallel with Rust and the rayon crate, and merge partial results into a final output. In real-world HPC or cloud settings, ephemeral containers can spin up and shut down per chunk, reducing cost while maintaining high throughput (Lee and Park, 2023).

The pipeline first splits the large input FASTA (genome.fa) into smaller chunk files via the splitFasta process, each sized around params.chunk_size base pairs; this enables parallelization by assigning each chunk to a separate HPC task or ephemeral container. Next, the scanMotif process invokes a Rust motif-scanning program on each chunk in parallel, with Nextflow automatically dispatching the chunks to executors such as local machines, HPC clusters, or cloud containers. Finally, the mergeHits process collects all partial JSON outputs into a single, consolidated file. In more advanced implementations, you might parse and merge each JSON file’s data or store the aggregated results in a database.

#### 2. Rust
The code belowshowcased here scans DNA sequences for TATA-like motifs in a robust and scalable way. It defines a customizable TATAPattern structure that can represent both exact and partial (mismatch-tolerant) TATA-box motifs, then uses parallel iteration (via the rayon crate) to rapidly process multiple sequences. This approach is suitable for larger-scale genomic data where performance and flexibility are essential—such as in scanning entire genomes or large collections of promoter regions for putative TATA boxes.

The TATAPattern struct stores the acceptable nucleotides for each position (e.g., ['T'] at position 0, ['A'] at position 1, etc.) and a maximum number of mismatches. The core function, find_tata_boxes, slides a window of the motif’s length across the input sequence, checking if each window contains fewer than or equal to the allowed number of mismatches relative to the pattern. By converting nucleotides to uppercase before matching, the code is case-insensitive. A parallel version, find_tata_boxes_parallel, leverages Rayon’s par_iter to distribute the workload across multiple CPU cores, enabling faster analysis when scanning a large set of sequences.

#### Files contents:
* experiment_43_5/
  * Cargo.toml (Cargo.toml file for dependencies)
  * motif_scanner (output file)
* experiment_43_5/src/
  * main.rs (rust script)
  * main.nf (nextflow script)
  * genome.fa (fasta file)
  * input.fa (fasta file)
  * output.json.rar (compressed output json)
  * tata_scan_merged.json.rar (compressed tata_scan_merged.json)

#### How to run:

run main.rs in powershell:

```powershell
cargo run --release -- "C:\Users\trian\BGVR\chapter_04\experiment_43_5\src\genome.fa" "C:\Users\trian\BGVR\chapter_04\experiment_43_5\src\output.json"
```

run main.nf in WSL:

```wsl
nextflow run main.nf --input_fasta genome.fa --chunk_size 1000000
```

#### [dependencies]

```toml
rayon = "1.10.0"
serde = { version = "1.0", features = ["derive"] }
serde_json = "1.0"
bio = "2.2.0"```

#### Explanation of the Output
My Nextflow pipeline and Rust motif scanner successfully ran and generated the following output:

##### 1. Motif Scanner Executable

* The pipeline compiled and executed the Rust motif-scanning program (experiment_43_5), which scanned the chunked FASTA sequences.
* The motif scanner applied a Position Weight Matrix (PWM)-based scoring algorithm to identify regions in the sequence that matched the given motif.

##### 2. Output File: tata_scan_merged.json
This file contains motif hits in JSON format. Each line represents a detected motif occurrence in the input sequence with:

* "seq_id": The identifier of the sequence from which the motif was found (e.g., "Synthetic").
* "position": The starting position of the motif in the sequence.
* "score": The PWM score calculated for this motif occurrence.

#### Interpretation of the Output
* The motif scanner identified multiple motif occurrences in the input sequence "Synthetic".
* The positions where the motifs were detected vary (e.g., 13, 14, 15, 16, 18, ...).
* The scores represent the strength of the match based on the PWM.
  * A higher score indicates a stronger match to the motif model.
  * Some floating-point precision artifacts (e.g., 2.0999999999999996) are present, which is normal in numerical computations.

#### Conclusion
The pipeline successfully completed the following steps:

* Chunked the input FASTA into smaller files for parallel processing.
* Processed each chunk in parallel using the Rust motif-scanning program.
* Merged the results into a single JSON file (tata_scan_merged.json).

#### Key Takeaways:

* The Nextflow pipeline properly orchestrates multiple processes (splitting, scanning, merging).
* The Rust program efficiently finds motifs using a PWM-based scanning approach.
* The JSON output correctly formats the detected motifs, making it easy to analyze further.
