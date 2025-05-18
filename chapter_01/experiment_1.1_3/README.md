## 1.1. Introduction to Rust Programming Language

### experiment_1.1_3

Here is a practical illustration of how Rust’s concurrency model can handle bioinformatics tasks in a parallel environment. For example, consider scanning a large set of DNA reads to tally the total number of occurrences of a particular motif. Using Rayon, you can safely distribute the workload across multiple threads without worrying about data races or improper memory sharing:

The code first opens a FASTA file using the bio crate and collects all DNA reads as Strings, which ensures each sequence is fully owned and easily shared among threads. We define a simple count_occurrences function that searches for a target motif (e.g., "GATTACA") in each sequence, allowing overlapping matches by shifting the search by just one character each time. After loading all sequences, we leverage Rayon’s par_iter() to parallelize the motif-counting step. Each thread independently processes a subset of reads, and .sum() consolidates the results into a final total. Thanks to Rust’s concurrency model and zero-cost iterator abstractions, this high-level code runs efficiently, with no extra runtime overhead compared to manually written loops or lower-level threading approaches.

#### Project Structure:

```plaintext
experiment_1.1_3/
└── Cargo.toml                     # Rust project configuration and dependencies
experiment_1.1_3/src/
├── main.rs                        # Main Rust script containing program logic
├── reads.fasta                    # FASTA file containing sequence reads
└── output.txt                     # Output file
```

#### How to run:

cargo run | tee output.txt

(run main.rs and save the output in output.txt)
  
#### [dependencies]

```toml
rayon = "1.10.0"
bio = "2.0.3"
```

#### Explanation of the Output
The program searches for occurrences of the DNA motif "GATTACA" within sequences read from a FASTA file (reads.fasta).

1. It loads all sequences from the file and stores them in memory.
2. It uses parallel processing (via Rayon) to efficiently search for the motif in each sequence.
3. The function count_occurrences allows overlapping matches, meaning if the motif appears multiple times in a sequence, including overlapping instances, they will all be counted.
4. The total count of "GATTACA" occurrences across all sequences is computed and printed.

##### The output:

Total occurrences of motif 'GATTACA' across all sequences: 621

* indicates that the motif "GATTACA" was found 621 times in the sequences present in reads.fasta.

#### Conclusion
The program successfully scanned sequencing data to detect occurrences of a specific DNA motif ("GATTACA"). By leveraging parallel computation, it efficiently processed large genomic sequences, enabling rapid motif detection. This approach is valuable in genomic analysis, regulatory sequence identification, and mutation studies.










