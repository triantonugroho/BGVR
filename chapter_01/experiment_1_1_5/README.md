## 1.1. Introduction to Rust Programming Language

### experiment_1_1_5

Below is a Rust code snippet that demonstrates generating a synthetic variation graph and processing it in a parallelized fashion. Each “node” of the graph can represent a simple genomic position or a base substitution variant. We use a customized rand environment to randomly produce these nodes, and then Rayon to distribute counting tasks—such as tallying how many positions vs. variants exist—across available CPU threads. The snippet also shows how pattern matching and Rust’s algebraic data types allow for explicit, compile-time distinctions between different kinds of genomic nodes, making the code safe, clear, and scalable.

The program defines an enum GraphNode with two variants, Position(usize) and Variant { from_base, to_base, pos }, capturing potential node types in a variation graph. In main, we initialize a random generator and build a vector of GraphNode objects, splitting them roughly half-and-half between position and variant nodes. Parallel aggregation then uses par_iter() from the Rayon crate, where each node is matched in a closure to extract whether it contributes to the position or variant count. This match expression ensures every branch is covered, providing compile-time safety. Finally, we illustrate further pattern matching and an additional parallel motif-counting demonstration. By chaining these approaches—random data generation, typed pattern matching, and parallel execution—the snippet highlights both Rust’s expressiveness and its performance through zero-cost abstractions.

#### Project Structure:

```plaintext
experiment_1_1_5/
├── Cargo.toml                     # Rust project configuration and dependencies
└── src/
    ├── main.rs                    # Main Rust script containing program logic
    ├── main.nf                    # Nextflow workflow script
    ├── example.fasta              # FASTA sequence file
    └── output.txt                 # Output file
```

#### Cargo.toml

```toml
[package]
name = "experiment_1_1_5"
version = "0.1.0"
edition = "2021"

[dependencies]
rayon = "1.10.0"
rand = "0.9.0"
```

#### How to run:

```powershell
cargo run | tee output.txt
```

(run main.rs and save the output in output.txt)
  

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










