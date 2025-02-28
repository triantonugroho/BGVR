# Contents

## 2.1. Introduction to Rust Programming Language

### experiment_21_1 
Rust code illustrating how to read all sequences from a FASTA file, store them in a vector of strings, and process them in parallel with Rayon.

### experiment_21_2 
A simple example that demonstrates Rust’s zero-cost abstractions in a genomic setting. 

### experiment_21_3 
A practical illustration of how Rust’s scanning a large set of DNA reads to tally the total number of occurrences of a particular motif.

### experiment_21_4 
The following Rust code processes each file format in parallel to find a specific DNA motif (e.g., “GATTACA”) in FASTAQ, BAM and VCF file.

### experiment_21_5 
Rust code snippet that demonstrates generating a synthetic variation graph and tallying how many positions vs. variants exist—across available CPU threads.

### experiment_21_6 
Rust project for Python interoperability, we can specify crate-type = ["cdylib"] to ensure our library is compiled as a dynamic library suitable for Python imports.

### experiment_21_7
A brief example illustrating how Rust can handle an imperative, pointer-based approach to building a simple suffix array.

### experiment_21_8
Demonstrates how Rust’s higher-order functions and iterators can apply a functional programming style to a set of genomic reads which does not rely on external crates, and can be executed with either rustc or cargo.

### experiment_21_9
A simple example demonstrating how a Rust library can use tch (the Rust bindings for PyTorch) and then expose this functionality to Python via PyO3.

### experiment_21_10
Rust code that reads each record from a FASTQ file via the Reader provided by the bio crate, storing the read ID and its length in a HashMap.

