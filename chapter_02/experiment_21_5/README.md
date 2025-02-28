## 2.1. Introduction to Rust Programming Language

### experiment_21_5

Below is a simple example that demonstrates Rust’s zero-cost abstractions in a genomic setting. This example reads a FASTA file and uses Rust’s iterator adapters (such as map, filter, and for_each) to process genomic sequences without incurring additional runtime overhead. Despite the high-level functional style, the compiler optimizes these operations down to efficient machine code, making them comparable to hand-written loops in languages like C or C++.

This code reads a FASTA file using the bio::io::fasta crate, which emits a stream of records (each containing a sequence). For each record, the sequence bytes are converted into a String via String::from_utf8_lossy, then any sequence under 50 nucleotides is discarded through a filter call. The remaining sequences move to a second map step that calculates their GC content by counting the characters ‘G’ or ‘C.’ Finally, the .sum() operation aggregates these individual GC counts into one total. Rust compiles this chain of iterator adapters into efficient loops without constructing intermediary collections, a technique referred to as “zero-cost abstraction.” Consequently, although the code is written in a high-level, functional style, it executes at speeds comparable to hand-tuned loops in lower-level languages.

#### Files contents:
* main.rs (rust script)
* main.nf (nextflow script)
* example.fasta (fasta file)
* Cargo.toml (Cargo.toml file)
* output.txt (output file)

#### How to run:

```powershell
cargo run | tee output.txt
```

(run main.rs and save the output in output.txt)
  
#### [dependencies]

```toml
rayon = "1.10.0"
rand = "0.9.0"
```

#### Explanation of the Output

##### 1. "Synthetic graph has 100 total nodes: 50 positions, 50 variants"

* The program generates a synthetic graph consisting of 100 nodes.
* Each node is randomly assigned as either a Position or a Variant, with approximately half falling into each category (50 each).
* Position nodes are represented by a unique number (e.g., 100, 102).
* Variant nodes represent genetic mutations, where one base is substituted for another at a specific position (e.g., A → C at position 201).

##### 2. "Showing first 5 nodes:"

* The program prints the first five nodes in the generated graph:

```sh
Position node at 100
Variant at 201: A -> C
Position node at 102
Variant at 203: A -> T
Variant at 204: G -> C
```

* This confirms that both Position and Variant nodes are correctly generated and stored.

##### 3. "Total occurrences of motif 'AC' in node representations: 3"

* The program searches for the motif "AC" in string representations of the nodes.
* The count_occurrences function looks for "AC" in:
  * Position node numbers (100, 102, etc.).
  * Variant node sequences (A -> C, G -> C, etc.).
* In this particular run, the motif "AC" appears 3 times across all nodes.

#### Conclusion
The program successfully generates a synthetic graph consisting of Position and Variant nodes, demonstrating random mutation events. It utilizes Rayon for parallel computation, efficiently analyzing the structure and searching for sequence motifs. The ability to model genetic variations in this format can support bioinformatics applications, including mutation analysis, sequence pattern detection, and graph-based genomic studies.










