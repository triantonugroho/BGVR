## 1.1. Introduction to Rust Programming Language

## experiment_1.1_1

Below is a Rust code illustrating how to read all sequences from a FASTA file, store them in a vector of strings, and process them in parallel with Rayon. This approach avoids the lifetime and ownership issues that arise when returning references to data held within a local mutex. By reading the file sequentially, converting each record’s bytes to UTF-8, and pushing the resulting strings into a Vec<String>, the code ensures that the data is fully owned, making it safe to be passed around or returned.

Once the sequences are loaded, the snippet uses Rayon’s parallel iterator (par_iter()) to concurrently compute a simple metric—GC content—for demonstration. Each sequence is examined character by character, counting the occurrences of ‘G’ or ‘C’. Because Rust’s ownership model and the Vec<String> both provide clear boundaries on data lifetime, this design not only simplifies concurrency but also allows for scalable processing of large data sets, provided there is sufficient memory to store all sequences. Below is the self-contained code:.

By wrapping the vector of sequences in a Mutex and then converting it into a local copy, Rust ensures that no two threads can mutate the data simultaneously, preventing data races. Any attempt to borrow the sequences vector incorrectly—say, by attempting a second mutable reference while another one is active—would fail at compile time. This guarantees that the code remains memory-safe, even under parallel workloads often found in genomics and functional biology pipelines.

## Project Structure:

```plaintext
experiment_1.1_1/src/
├── main.rs                    # Main Rust script containing program logic
├── example.fasta              # FASTA file
└── output.txt                 # Output file
experiment_1.1_1/               # Empty or duplicate project directory
Cargo.toml                     # Rust project configuration and dependencies file (at root level)
```

## How to run:

cargo run | tee output.txt

(run main.rs and save the output in output.txt)
  
## [dependencies]

```toml
bio = "2.0.3"
rayon = "1.10.0"
```

## Explanation of the Output

Your Rust program reads a FASTA file, processes its sequences in parallel, and writes results to an output file (output.txt). The key metrics computed for each sequence are:
1. Length: The total number of nucleotides in the sequence.
2. GC Content: The count of 'G' (Guanine) and 'C' (Cytosine) bases in the sequence.

### 1. Given the output:

Length: 10000000, GC: 4999007

Successfully processed 1 sequences.

* Length: 10,000,000 → The sequence read from the FASTA file contains 10 million nucleotides.
* GC: 4,999,007 → About 49.99% of the sequence consists of G or C bases.
*Successfully processed 1 sequences. → The FASTA file contained only one sequence, which was fully analyzed.

### 2. How This Happens in the Code:

1. Reading the FASTA File (load_sequences())
   * The program opens example.fasta, reads all sequences, and stores them in a Vec<String>.
    * Each sequence is converted into a String for easy manipulation.
2. Parallel Processing with Rayon (par_iter())
  * Each sequence is processed in parallel using rayon::prelude::*.
  * It calculates the length (seq.len()) and counts the number of G and C bases.
3. Writing to output.txt
  * Each result (Length: ..., GC: ...) is written to the file.
  * The program logs the total number of sequences processed.

## Conclusion

The program successfully analyzed a FASTA-formatted sequencing dataset by computing sequence length and GC content in parallel using Rayon. This efficient processing enables downstream applications such as genome composition analysis, species classification, and quality assessment in bioinformatics workflows.

