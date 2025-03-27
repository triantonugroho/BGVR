## 5.2. Preprocessing and Quality Control

### experiment_52

The Rust program below applies a sliding window trimming algorithm to reads from a FASTQ file, chunk by chunk, and writes partial results to disk. By processing only a limited subset of reads at a time (chunk_size) and using parallel iteration through Rayon’s data parallelism, it addresses memory constraints that can arise with extremely large datasets. In addition, the partial outputs are serialized in JSON, creating a natural checkpointing mechanism that is particularly useful when dealing with large volumes of sequence data in a production or HPC environment.

The sliding window trimming approach itself scans each read’s quality scores from both ends. It uses a configurable window size and average Phred score threshold to determine how far inwards to trim low-quality bases. This design ensures that reads with persistently low-quality regions at the start or end are trimmed back to a region of higher confidence. By decoupling the trimming logic from the parallel iteration, one can easily substitute a more sophisticated algorithm or integrate advanced crates such as linfa for machine learning–based trimming heuristics if desired.

After parsing command-line arguments using the clap crate, the program uses needletail to open and read the FASTQ file, which can be either plain or gzipped. Instead of loading all records at once, it takes in batches of size chunk_size. Each batch is processed via Rayon’s par_iter(), thereby distributing the trimming workload across available CPU cores.

The trimming logic is encapsulated in the sliding_window_trim function, which uses a rolling sum of Phred scores to compute the average quality in a window. If this average drops below the specified threshold, the function updates the boundaries that define the start or end of the trimmed region. The result is a pair of sliced vectors: the retained sequence bases and the corresponding quality scores.

Once the trimming is complete for each batch, the program writes a JSON file containing all trimmed reads for that chunk. This partial output strategy is highly advantageous in large-scale operations. It lets the pipeline resume from the last completed chunk if the process is interrupted, and it mitigates memory spikes by limiting how many records are held in memory simultaneously. Finally, the program reports basic statistics, including how many partial files were written. A subsequent step in the pipeline could merge or further analyze these JSON files, thereby enabling a flexible, HPC-friendly workflow.

#### Files contents:
* experiment_52/
  * Cargo.toml (Cargo.toml file for dependencies)
*experiment_52/src/
  * main.rs (rust script)
  * output.txt (output file)

#### How to run:

run in powershell:

```powershell
cargo run -- --input C:\Users\trian\BGVR\chapter_05\experiment_52\src\example.fastq
```

(run main.rs using input data path)
  
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

##### 1. Integer Count Estimation

```rust
Actual integer count: 10000
Estimated integer count: 9560.80
```

* A dataset containing 10,000 unique integers (0 to 9,999) is processed using HyperLogLog with precision p=10.
* The estimated count is 9560.80, which is close but slightly under the actual count due to the probabilistic nature of HyperLogLog.
* The estimation error is (10000 - 9560.80) / 10000 = 4.39%, which is expected as HyperLogLog is an approximation algorithm.

##### 2. String Count Estimation

```rust
Actual string unique count: 5
Estimated string count: 13.51
```

* A list of strings is processed: ["apple", "banana", "cherry", "banana", "date", "apple", "elderberry"].

* The actual number of unique strings is 5 ({"apple", "banana", "cherry", "date", "elderberry"}).

* The estimated count is 13.51, which is significantly higher than the actual value.

* The overestimation is likely due to the low precision (p=4), meaning there are only 16 registers (2^4 = 16), leading to more collisions and increased variance.

##### 3. Merging Two HyperLogLogs

```rust
Merging two HyperLogLogs each containing half of the integer range:
Merged estimate of unique integers: 9560.80
Actual unique count (0..10000): 10000
```

* The dataset of 10,000 integers is split into two halves (0..5000 and 5000..9999).

* Two separate HyperLogLogs are created for each half.

* After merging them, the estimated count remains 9560.80, the same as the previous estimate.

* This confirms that merging two HyperLogLogs does not overcount, as they retain only the maximum register values.

#### Conclusion

* HyperLogLog provides an efficient way to estimate cardinality (unique elements) using limited memory.

* The estimation is close but not exact—it depends on the precision (p value) and the nature of the data.

* Higher precision (p) reduces error but requires more memory.

* Smaller datasets (like strings with p=4) suffer from higher variance, leading to more overestimation or underestimation.

* Merging HyperLogLogs maintains accuracy, showing its usefulness for distributed data aggregation.

For better accuracy in real-world applications, choosing an appropriate p value is crucial, balancing memory usage and estimation precision.


