## 2.1. Introduction to Rust Programming Language

### experiment_21_8

Demonstrates how Rust’s higher-order functions and iterators can apply a functional programming style to a set of genomic reads. The code below compiles out of the box, does not rely on external crates, and can be executed with either rustc or cargo.

This Rust program demonstrates how to calculate the total GC content across multiple DNA reads. It starts with a small collection of DNA reads, each represented by a string. For each read, the code counts the occurrences of 'G' and 'C' by filtering characters in the string, then sums (or “folds”) all those counts together to produce a single total. Finally, it prints the total GC content, offering a simple yet clear illustration of Rust’s functional-style iteration methods (map and fold) applied to biological sequence data.

At the systems level, Rust still offers low-level features. For example, if you wanted to hand-optimize k-mer counting by directly manipulating bits to compress sequences, Rust’s type system allows you to define bitfields or custom pointer arithmetic. That level of control is reminiscent of classical HPC approaches in C, but Rust’s compile-time checks prevent improper usage of allocated memory.

#### Files contents:
* main.rs (rust script)
* Cargo.toml (Cargo.toml file)
* output.txt (output file)

#### How to run:

```rust
cargo run | tee output.txt
```

(run main.rs and save the output in output.txt)
  
#### [dependencies]

No dependencies.

#### Explanation of the Output

##### Understanding the Code
This Rust program calculates the total GC content (number of G and C bases) in a collection of short DNA reads.

###### 1. DNA Reads Collection

```rust
let reads = vec!["ATCG", "TTGA", "CGTA", "GGTT"];
```

* This vector contains four short DNA sequences: "ATCG", "TTGA", "CGTA", and "GGTT".
* Each sequence is a string representing a segment of DNA.

###### 2. Counting GC Bases in Each Read

```rust
let total_gc = reads
     .iter()
     .map(|read| {
         read.chars()
             .filter(|&c| c == 'G' || c == 'C')
             .count()
     })
```

* ```rust.iter()``` iterates over all reads.
* ```rust.map()``` transforms each read by counting the number of 'G' and 'C' bases.
* ```rust.filter(|&c| c == 'G' || c == 'C')``` keeps only the G and C characters.
* ```rust.count()``` returns the number of G and C bases in each read.

###### 3. Summing Up GC Counts

```rust
.fold(0, |acc, gc_count| acc + gc_count);
```

* ```rust.fold(0, |acc, gc_count| acc + gc_count)``` accumulates the total GC count by summing up individual counts.

##### Breaking Down the Output Calculation

Each read and its GC count:

```sh
Read	  G Content	C Content	GC Count
"ATCG"	1	        1	        2
"TTGA"	1	        0	        1
"CGTA"	1	        1	        2
"GGTT"	2	        0	        2
```

Summing these values:
2 + 1 + 2 + 2 = 7

Thus, the program outputs:

```sh
Total GC content in all reads: 7
```

#### Conclusion
* This program efficiently calculates the total GC content in a given set of DNA reads.
* GC content is an important metric in bioinformatics because it affects DNA stability and genome characteristics.
* The use of functional programming (map, filter, fold) makes the solution concise and efficient.
* This approach can be extended to handle larger datasets or parallelized using Rayon for better performance.


