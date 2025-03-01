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

## 2.3. Data Structures and Algorithms in Rust

### experiment_23_1
A succinct Rust Classical substring-search algorithms like Knuth–Morris–Pratt (KMP) implementation.

### experiment_23_2
Rust code of a mini-example combining concurrency, the bio crate for reading FASTA, and partial de Bruijn graph assembly.

## 2.4. AI/ML Implementation with Rust Crates

### experiment_24
A sample Rust code snippet demonstrating how to do K-means clustering with linfa and then, in a separate pipeline, train a neural network with tch-rs.

## 2.5. Acquiring and Cleaning Data

### experiment_25_1
The complete Nextflow script that demonstrates a streamlined pipeline without containers. It downloads .sra files for a list of accessions, performs a simulated checksum verification, and moves the final outputs to a specified directory. 

### experiment_25_2
A single Nextflow script that uses pinned Rust and Samtools versions as if they were installed on our system, references numeric parameters and environment variables, and ensures that future runs use the same conditions by explicitly specifying all necessary details in the pipeline script.

### experiment_25_3
A streamlined Rust tool that demonstrates how to read a FASTQ file using bio and apply a simple filter. It also includes logging for HPC environments, using crates such as env_logger to track progress and potential errors. 

### experiment_25_4
Rust code to parse metadata from a CSV file that describes sample conditions (e.g., diseased vs. healthy) and merge it with alignment statistics. 

## 2.6. Scientific Computation Workflow with Rust and Nextflow

### experiment_26_1
Nextflow DSL2 script demonstrating how to run a Rust read-trimming tool in parallel across multiple FASTQ files. The Rust trimmer itself can be a simple program that removes low-quality bases from both ends of a read.

### experiment_26_2
An example sequence might begin with a Rust-based alignment tool that transforms FASTQ to BAM files, followed by a coverage parser that generates CSV files, and concluding with a machine learning step that trains or applies a model to those coverage statistics.
