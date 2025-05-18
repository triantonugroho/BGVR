## 1.3. Data Structures and Algorithms in Rust

### experiment_1.3_1

Classical substring-search algorithms like Knuth–Morris–Pratt (KMP) and Boyer–Moore remain core building blocks in bioinformatics. While tools like BLAST and Bowtie handle complex alignments, simpler tasks—like scanning for motifs or barcodes—often rely on well-known linear or sublinear time solutions (Knuth et al. 1977; Gusfield 1997). Below is a succinct Rust KMP implementation.

The KMP search’s O(n + m) time complexity is vital for high-throughput motif discovery, primer matching, and other foundational string tasks in functional genomics.

#### Project Structure:

```plaintext
experiment_1.3_1/
├── Cargo.toml                     # Rust project configuration and dependencies
└── src/
    ├── main.rs                    # Main Rust script containing program logic
    └── output.txt                 # Output file
```

#### Cargo.toml

```toml
[package]
name = "experiment_1.3_1"
version = "0.1.0"
edition = "2021"

[dependencies]

```

#### How to run:

run in powershell:

```powershell
cargo run | tee output.txt
```

(run main.rs and save the output in output.txt)
  

#### Explanation of the Output
This Rust program implements the Knuth-Morris-Pratt (KMP) algorithm to efficiently search for occurrences of a pattern (ATGC) within a given DNA sequence (ATGCGATATCGATGCGATGCGATGC). The program identifies all the positions in the sequence where the pattern appears.

##### Step-by-Step Execution

###### 1. Building the Prefix Table (build_prefix_table)

* The prefix table (pi array) helps optimize pattern searching by avoiding redundant comparisons.
* It preprocesses the pattern "ATGC" and stores information about repeated prefix-suffix matches.
* The computed prefix table for "ATGC" is [0, 0, 0, 0], meaning there are no repeating prefix-suffix overlaps.

###### 2. Searching the Pattern in the DNA Sequence (kmp_search)

* The function scans the dna_sequence using efficient character comparisons.
* It skips unnecessary comparisons by using the prefix table.
* Whenever a full match is found, it records the starting index of the match.

###### 3. Running the Program (main)

* The DNA sequence: "ATGCGATATCGATGCGATGCGATGC"
* The pattern: "ATGC"
* The KMP search finds the pattern at positions 0, 11, 16, and 21.

##### Output Interpretation

```sh
Pola 'ATGC' ditemukan pada indeks: [0, 11, 16, 21]
```

* "ATGC" is found at index 0 (start of the string).
* "ATGC" is found again at index 11, 16, and 21.
  
#### Conclusion
* The KMP algorithm efficiently finds multiple occurrences of a pattern in a string.
* It avoids unnecessary comparisons by utilizing the prefix table.
* This approach is useful in bioinformatics for DNA sequence analysis, where searching motifs in long genomic sequences is crucial.
* The algorithm runs in O(n + m) time complexity, making it much faster than a brute-force search.
