## 2.1. Introduction to Rust Programming Language

### experiment_21_2

Below is a simple example that demonstrates Rust’s zero-cost abstractions in a genomic setting. This example reads a FASTA file and uses Rust’s iterator adapters (such as map, filter, and for_each) to process genomic sequences without incurring additional runtime overhead. Despite the high-level functional style, the compiler optimizes these operations down to efficient machine code, making them comparable to hand-written loops in languages like C or C++.

This code reads a FASTA file using the bio::io::fasta crate, which emits a stream of records (each containing a sequence). For each record, the sequence bytes are converted into a String via String::from_utf8_lossy, then any sequence under 50 nucleotides is discarded through a filter call. The remaining sequences move to a second map step that calculates their GC content by counting the characters ‘G’ or ‘C.’ Finally, the .sum() operation aggregates these individual GC counts into one total. Rust compiles this chain of iterator adapters into efficient loops without constructing intermediary collections, a technique referred to as “zero-cost abstraction.” Consequently, although the code is written in a high-level, functional style, it executes at speeds comparable to hand-tuned loops in lower-level languages.

#### Files contents:
* main.rs (rust script)
* reads.fasta (fasta file)
* Cargo.toml (Cargo.toml file)
* output.txt (output file)

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










