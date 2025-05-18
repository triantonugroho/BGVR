## 4.3. Motif Discovery and Regulatory Element Identification

### experiment_43_4

The code belowshowcased here scans DNA sequences for TATA-like motifs in a robust and scalable way. It defines a customizable TATAPattern structure that can represent both exact and partial (mismatch-tolerant) TATA-box motifs, then uses parallel iteration (via the rayon crate) to rapidly process multiple sequences. This approach is suitable for larger-scale genomic data where performance and flexibility are essential—such as in scanning entire genomes or large collections of promoter regions for putative TATA boxes.

The TATAPattern struct stores the acceptable nucleotides for each position (e.g., ['T'] at position 0, ['A'] at position 1, etc.) and a maximum number of mismatches. The core function, find_tata_boxes, slides a window of the motif’s length across the input sequence, checking if each window contains fewer than or equal to the allowed number of mismatches relative to the pattern. By converting nucleotides to uppercase before matching, the code is case-insensitive. A parallel version, find_tata_boxes_parallel, leverages Rayon’s par_iter to distribute the workload across multiple CPU cores, enabling faster analysis when scanning a large set of sequences.

#### Project Structure:

```plaintext
experiment_43_4/
├── Cargo.toml                     # Rust project configuration and dependencies
└── src/
    ├── main.rs                    # Main Rust script containing program logic
    └── output.txt                 # Output file
```

#### How to run:

run main.rs in powershell:

```powershell
cargo run | tee output.txt
```
(run main.rs and save the output in output.txt)

#### [dependencies]

```toml
rayon = "1.10.0"
```

#### Explanation of Output
The main.rs code implements a TATA box detection algorithm in DNA sequences using pattern matching with support for mismatches. Here's a detailed explanation of the process and output:

##### Input Data

###### 1. DNA Sequences
The input consists of four DNA sequences:

```rust
let seqs = vec![
    "GGTTTATATAAACTATAATTTTACGT".to_string(), // Sequence 0
    "tttatacccggttttataAa".to_string(),       // Sequence 1
    "AAAAATATA".to_string(),                  // Sequence 2
    "NoTATAhere".to_string(),                 // Sequence 3
];
```

###### 2. Pattern Definitions
Two patterns are defined:

* Default TATA Pattern

  * Fixed consensus "TATA(A|T)A"
  * No mismatches allowed

* Custom Pattern
  * Same core consensus "TATA(A|T)A"
  * Allows one mismatch in the pattern

#### How the Algorithm Works

##### Pattern Matching Process:

1. The algorithm converts the input sequences to uppercase (case-insensitive).
2. It slides a window of size 6 (pattern length) across each sequence.
3. It checks each window against the consensus pattern:
* A window matches if all nucleotides are identical to the consensus OR within the allowed mismatch limit.
4. Matching is performed using matches_window():
* If the number of mismatches exceeds the allowed limit → reject the match.
* If within the allowed limit → accept the match.

##### Parallel Processing:
* The search is executed in parallel using Rayon.
* Each sequence is processed concurrently, improving performance on large datasets.

#### Output Breakdown
Input Sequences:

```rust
Sequences:
[
    "GGTTTATATAAACTATAATTTTACGT",  // Seq 0
    "tttatacccggttttataAa",         // Seq 1
    "AAAAATATA",                    // Seq 2
    "NoTATAhere",                   // Seq 3
]
```

##### Default Pattern Matches (No Mismatch Allowed)
Sequence	Match Positions	Explanation
Seq       0	    [4, 6]	  TATAAA at index 4, TATAAT at index 6
Seq       1	    [14]	    TATAAa at index 14 matches the consensus (case-insensitive)
Seq 2	    []	  No        perfect `TATA(A
Seq 3	    []	  No        perfect match

##### Custom Pattern Matches (Allows 1 Mismatch)
Sequence	Match Positions	      Explanation
Seq       0	    [2, 4, 6, 13]	  Matches even with one mismatch allowed
Seq       1	    [0, 12, 14]	    TATACC, TATA, and TATAAa match
Seq       2	    [3]	            ATATA matches with one mismatch allowed
Seq       3	    []	            No TATA-like sequence close enough

#### Interpretation
* The default pattern only detects strict TATA(A|T)A matches.
* The custom pattern increases sensitivity by allowing mismatches:
  * Catches near-matches and variations of TATA motifs.
  * Useful for identifying biologically relevant TATA-like sequences.

#### Conclusion
* The algorithm successfully:
  * Handles both exact and approximate pattern matching.
  * Supports mismatch tolerance, increasing sensitivity.
  * Uses parallel processing for faster analysis on large datasets.
  * Matches are case-insensitive and robust to sequence variation.
