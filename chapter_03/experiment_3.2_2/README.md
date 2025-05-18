## 3.2. Sequence Data Structures and Strings Algorithms

### experiment_3.2_2

This Rust code demonstrates how to build and merge partial suffix arrays for genomic data in an MPI environment by having each rank receive a segment of the reference (or read) text, construct a local suffix array, and then gather all partial results on rank 0 for global merging, enabling parallelization of suffix array construction—particularly important for large genomes where a single node might encounter time or memory constraints. It further illustrates a simplified HPC scenario for building partial suffix arrays or BWT segments across multiple nodes, merging them into a unified index on rank 0. Each node obtains a portion of the text, constructs a local suffix structure (such as a suffix array or BWT), and sends it to the root node, which combines these partial outputs. Although placeholders are used for demonstration, the example highlights how Rust’s safety guarantees, concurrency features, and MPI-based communication can be harnessed for large-scale string indexing tasks such as suffix array or BWT construction.

Rank 0 determines chunk boundaries by evenly partitioning the text across size ranks, then each rank extracts its assigned substring. The function build_local_suffix_array naively sorts suffix positions based on their lexicographical order, returning a PartialSuffixArray containing those positions plus an offset to indicate where the chunk appears in the global text. A gather operation collects these partial structures at rank 0, which then merges them into a single global suffix array by comparing each suffix substring in a BTreeMap. Because the merging step sorts suffixes across all chunks, the final output is a global ordering of suffix positions that can be used as a baseline for downstream tasks (e.g., LCP array or BWT construction).

#### Project Structure:

```plaintext
experiment_3.2_2/
└── Cargo.toml                     # Rust project configuration and dependencies
experiment_32_3/src/
├── main.rs                        # Main Rust script containing program logic
└── output.txt                     # Text output file
```

#### Cargo.toml

```toml
[package]
name = "experiment_3.2_2"
version = "0.1.0"
edition = "2021"

[dependencies]
bincode = "1.3.3"
mpi = "0.8.0"
serde = { version = "1.0", features = ["derive"] }
```

#### How to run:

run in powershell:

```powershell
cargo run | tee output.txt
```

(run main.rs and save the output in output.txt)
  

#### Explanation of the Output:
The program is designed to create a suffix array for a given text using parallel computing (MPI). Here's a detailed breakdown of the output:

##### Suffix Array:
A suffix array is a sorted array of all possible suffixes of a string. Each element in the array is the index where the corresponding suffix starts. For example, for the string "GATTACAGATTACACAT", the suffixes would be:

* GATTACAGATTACACAT
* ATTACAGATTACACAT
* TTACAGATTACACAT
* ...
The suffix array is simply the sorted list of indices of these suffixes. In the output, each index corresponds to the position in the text where the respective suffix starts.

##### Breakdown of the Output:

```rust
Global Suffix Array: [11, 4, 13, 6, 15, 8, 1, 12, 5, 14, 7, 0, 16, 10, 3, 9, 2]
```

* This array lists the starting positions of the sorted suffixes.
* For instance:
  * The suffix starting at index 11 is "CAGATTACACAT".
  * The suffix starting at index 4 is "ACAT".
  * The suffix starting at index 13 is "AT".
  * And so on.

Each number in the array represents the index of the starting position of a suffix in the string.

##### Steps Performed:
###### 1. Text Partitioning:

* The text "GATTACAGATTACACAT" is divided into chunks and distributed across multiple ranks (MPI processes).
* For example, if there are two ranks, rank 0 might handle the first half of the text, and rank 1 handles the second half.

###### 2. Local Suffix Array Construction:

* Each rank constructs a local suffix array for its chunk of the text. The suffixes within the chunk are sorted lexicographically.

###### 3. Communication between Ranks:

* The partial suffix arrays (constructed by each rank) are serialized and sent to rank 0.

###### 4. Merging:

* Rank 0 receives the partial suffix arrays from all ranks and combines them into a global suffix array.
* The global suffix array is built by merging the local suffix arrays from all ranks. This is done by mapping each local suffix to its global position in the text using the offsets.

###### 5. Final Output:

* After merging the results from all ranks, rank 0 prints the final global suffix array, which represents the sorted order of all suffixes in the entire text.

#### Conclusion:
The output is the sorted order of suffixes of the string "GATTACAGATTACACAT". Each index in the result corresponds to the starting position of a suffix in the original string. The program leverages parallel computing (using MPI) to divide the task of creating the suffix array into chunks, each processed by a different rank, and then collects the results to form the complete suffix array.

* Global Suffix Array: The printed output provides the final sorted order of all suffixes of the text, with each index corresponding to the starting position of the sorted suffix in the original string.
* Scalability: This approach is scalable because it divides the work of sorting the suffixes among multiple ranks, which can be useful for large datasets.
* Parallelism: The program demonstrates how parallel computing can be applied to efficiently process large text data by distributing the workload, merging results, and leveraging collective communication.
