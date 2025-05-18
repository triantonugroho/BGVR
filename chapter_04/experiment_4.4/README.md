## 4.4. Epigenomic Data Integration and Algorithms

### experiment_4.4

The following Rust code snippet demonstrates a streamlined approach for calling peaks (e.g., regions of high signal intensity) on genomic coverage data, which might arise from experiments like ChIP-seq or ATAC-seq. It addresses real-world tasks where large chromosomal data must be efficiently scanned to identify potential regulatory regions. To ensure scalability in high-performance computing (HPC) or cloud settings, the solution uses the Rayon library for parallel processing and employs a flexible smoothing function to reduce noise in coverage profiles.

After optionally smoothing the coverage data with a rolling average (using a prefix-sum approach), the code slides a user-defined window across each chromosome to compute the mean coverage within that window. Whenever this mean exceeds the specified threshold, it marks the position as a peak. Rayon’s parallel iterator (par_iter) automatically distributes the workload across CPU cores, processing each chromosome (or sub-chromosome data) concurrently. Finally, the identified peaks are gathered into a single collection and written to a file in a simplified BED-like format, ready for further analysis or merging with other datasets.

The resulting partial_peaks.bed file (or multiple such files, one per node) can be merged in a final HPC pipeline step. That step might also incorporate additional heuristics, such as removing spurious peaks or combining nearby calls. This partial merges pattern is standard in large bioinformatics workflows and is readily integrated with container orchestration tools like Nextflow or Snakemake. Because Rust compiles to efficient, static binaries, containers remain lightweight and reproducible, an essential aspect of modern HPC pipelines (Davis et al. (2023)).

#### Project Structure:

```plaintext
experiment_4.4/
├── Cargo.toml                     # Rust project configuration and dependencies
└── src/
    ├── main.rs                    # Main Rust script containing program logic
    └── partial_peaks.bed          # BED format output file containing partial peaks
```

#### Cargo.toml

```toml
[package]
name = "experiment_4.4"
version = "0.1.0"
edition = "2024"

[dependencies]
rayon = "1.10.0"
```

#### How to run:

run main.rs in powershell:

```powershell
cargo run
```
(run main.rs and get the output partial_peaks.bed)


#### Explanation of the Output
My Rust program performs the following steps:

##### 1. Simulates chromosome coverage data

* Two chromosomes (chr1 and chr2) with predefined coverage values.
* Example values for chr2: [0.0, 7.5, 8.0, 6.2, 2.1, 9.4, 10.2, 0.5].

##### 2. Processes the data to detect peaks

* Uses a rolling-window smoothing function (window size = 3).
* Identifies local peaks where the mean coverage in the window exceeds 3.0.

##### 3. Outputs peaks into partial_peaks.bed

* The .bed format (tab-separated) contains:
  * Chromosome (e.g., chr2)
  * Position in the sequence
  * Peak intensity (after smoothing)
    
#### Interpretation of the Output
The partial_peaks.bed file lists detected peaks for chr2, but not for chr1, because:

* chr1's coverage values are lower, and no values exceeded the threshold (3.0).
* chr2 had multiple positions where the smoothed coverage was above the threshold.

##### Example Breakdown:

1. First peak at position 0 on chr2

* Smoothed coverage = 4.458 (above threshold 3.0).

2. Peak at position 2 on chr2
* Smoothed coverage = 5.944.

3. Strongest peak at position 5 on chr2
* Smoothed coverage = 6.611.

#### Conclusion
The program successfully identified peaks based on smoothed coverage data:

* Peak detection works correctly: Peaks appear only where smoothed coverage is above 3.0.
* Parallel processing is effective: The computation is optimized using rayon.
* No peaks in chr1, which confirms that the thresholding logic is working.
