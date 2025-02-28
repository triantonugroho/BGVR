## 2.1. Introduction to Rust Programming Language

### experiment_21_7

Here is a brief example illustrating how Rust can handle an imperative, pointer-based approach to building a simple suffix array. This style is reminiscent of C++ in that it manipulates raw pointers for performance-critical tasks. Rust’s compiler still enforces strict rules about the pointer’s lifetime, preventing you from using it once the underlying data goes out of scope.

In this snippet, build_suffix_array manually handles pointer arithmetic to compare suffixes, reminiscent of a lower-level C++ style. The reference genome data is stored in a Vec, and seq.as_ptr() provides a raw pointer to the first byte. Because Rust attaches the pointer’s lifetime to the underlying Vec, any misuse—such as accessing memory after seq goes out of scope—would result in a compile-time error. The compare_suffixes function uses unsafe code, but the scope of that unsafety is minimized. Rust checks that the pointer remains valid as long as the Vec has not been dropped, preserving memory safety guarantees even though we resort to pointer-based comparisons for performance or fine-grained control.

## Files contents:
* main.rs (rust script)
* main.nf (nextflow script)
* example.fasta (fasta file)
* Cargo.toml (Cargo.toml file)
* output.txt (output file)

## How to run:

```rust
cargo run | tee output.txt
```

(run main.rs and save the output in output.txt)
  
## [dependencies]

```toml
bio = "2.0.3"
rayon = "1.10.0"
rust-htslib = "0.49.0"
```
## Explanation of the Output

The program reads a FASTA file and processes each sequence to compute the total GC content (the number of guanine (G) and cytosine (C) bases). It filters sequences that are at least 50 nucleotides long before counting GC bases. The output:


Total GC content in sequences >= 50 nt: 4999007

* Indicates that the total number of GC bases in all sequences meeting the length requirement is 4,999,007.

## Conclusion
The program successfully analyzed sequencing data by filtering short sequences and computing total GC content. This approach enables applications such as genome composition analysis, species identification, and quality control in bioinformatics workflows.










