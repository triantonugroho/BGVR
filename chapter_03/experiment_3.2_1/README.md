## 3.2. Sequence Data Structures and Strings Algorithms

### experiment_32_1

This code demonstrates a simple, parallel approach for detecting a known genomic “pattern” (e.g., a short sequence) within large genomic data using MPI and the Knuth-Morris-Pratt (KMP) algorithm. By splitting the genomic reference (or read data) across multiple ranks, each node handles a portion of the text plus any necessary overlap, ensuring that pattern matches spanning chunk boundaries are not missed. This design scales effectively across many ranks, which is common in HPC environments where genomic references can be gigabytes or terabytes in size.

Rank 0 loads the genomic data from a file, computes chunk boundaries, and distributes each chunk (plus overlap, if needed) to the other ranks. Each rank then performs a local KMP search for the target sequence on its assigned chunk and returns partial match indices back to rank 0. Because local match positions are relative to each chunk, rank 0 offsets them properly when assembling the final list of matches. The partial matches are then merged and sorted, yielding a global view of where the target genomic sequence occurs.

#### Project Structure:

```plaintext
experiment_32_1/
├── Cargo.toml                            # Rust project configuration and dependencies
└── src/
    ├── main.rs                           # Main Rust script containing program logic
    ├── big_text_example.txt              # Text file dataset
    ├── generate big_text_example.txt.ipynb  # Python notebook to generate the text dataset
    └── output.txt                        # Text output file
```

#### How to run:

run in powershell:

```powershell
cargo run | tee output.txt
```

(run main.rs and save the output in output.txt)
  
#### [dependencies]

```toml
rand = "0.9.0"
rayon = "1.10.0"
```

#### Explanation of the Output:
The output of the program represents the global positions in the text where the given pattern "ABABABC" was found. These positions are the result of a parallelized string matching algorithm using the Knuth-Morris-Pratt (KMP) pattern matching algorithm, executed across multiple processes using MPI (Message Passing Interface).

Here's a breakdown of the key components in the output:

##### 1. Global Matches for Pattern 'ABABABC':

* The pattern "ABABABC" was found multiple times throughout the entire text.
* The positions listed are the starting indices in the text where this pattern matches.
* The list contains a total of 200 positions where the pattern starts in the text, with each position incrementing by 5000, showing a regular pattern in the text data.
* Each number in the output is an index in the text where the pattern "ABABABC" starts.

##### 2. Text Processing and Parallelization:

* The program is utilizing multiple processes (through MPI) to split the workload of searching for the pattern in a large text file.
* Rank 0 is responsible for distributing the text into chunks for each rank (processor), which then performs the KMP search on their chunk.
* After performing the search, each rank sends its local match positions to rank 0, which merges all the results.

##### 3. Pattern Matching:

* The KMP algorithm is used to efficiently find occurrences of the pattern "ABABABC" in the text. The KMP algorithm builds a prefix table that allows it to skip unnecessary comparisons and directly match parts of the pattern.
* Once a match is found, the starting index (position) of the match in each chunk is stored.
* These positions are eventually combined into a final sorted list of global positions where the pattern occurs in the whole text.

#### Conclusion:
The output confirms that the pattern "ABABABC" is repeated regularly throughout the text, with matches at positions such as 4993, 9993, 14993, and so on. The matches are spaced roughly 5000 characters apart, indicating that the pattern might be embedded in the text at consistent intervals.

#### Key takeaways:

* The program uses parallel computing (via MPI) to efficiently search through large texts, distributing the work across multiple processes.
* The KMP algorithm allows for fast and efficient pattern matching, even in large datasets.
* The output represents the global positions of the pattern matches in the entire text after combining the results from all ranks.
This solution is highly scalable and can be used to search for patterns in very large datasets, making it suitable for text mining and other applications that require fast pattern matching.


