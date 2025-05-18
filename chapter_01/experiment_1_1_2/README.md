## 1.1. Introduction to Rust Programming Language

### experiment_1.1_2

Below is a simple example that demonstrates Rust’s zero-cost abstractions in a genomic setting. This example reads a FASTA file and uses Rust’s iterator adapters (such as map, filter, and for_each) to process genomic sequences without incurring additional runtime overhead. Despite the high-level functional style, the compiler optimizes these operations down to efficient machine code, making them comparable to hand-written loops in languages like C or C++.

This code reads a FASTA file using the bio::io::fasta crate, which emits a stream of records (each containing a sequence). For each record, the sequence bytes are converted into a String via String::from_utf8_lossy, then any sequence under 50 nucleotides is discarded through a filter call. The remaining sequences move to a second map step that calculates their GC content by counting the characters ‘G’ or ‘C.’ Finally, the .sum() operation aggregates these individual GC counts into one total. Rust compiles this chain of iterator adapters into efficient loops without constructing intermediary collections, a technique referred to as “zero-cost abstraction.” Consequently, although the code is written in a high-level, functional style, it executes at speeds comparable to hand-tuned loops in lower-level languages.

#### Project Structure:

```plaintext
experiment_1.1_2/
├── Cargo.toml                     # Rust project configuration and dependencies
└── src/
    ├── main.rs                    # Main Rust script containing program logic
    ├── example.fasta              # Fasta file
    └── output.txt                 # Output text file
 ```

#### How to run:

cargo run | tee output.txt

(run main.rs and save the output in output.txt)
  
#### Cargo.toml

```toml
[package]
name = "experiment_1.1_2"
version = "0.1.0"
edition = "2021"

[dependencies]
bio = "2.0.3"
```
#### Explanation of the Output

The program reads a FASTA file and processes each sequence to compute the total GC content (the number of guanine (G) and cytosine (C) bases). It filters sequences that are at least 50 nucleotides long before counting GC bases. The output:


Total GC content in sequences >= 50 nt: 4999007

* Indicates that the total number of GC bases in all sequences meeting the length requirement is 4,999,007.

#### Conclusion
The program successfully analyzed sequencing data by filtering short sequences and computing total GC content. This approach enables applications such as genome composition analysis, species identification, and quality control in bioinformatics workflows.









