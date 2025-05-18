# 3. Data Structures and Algorithms for Bioinformatics

## 3.1. Introduction to Data Structures and Algorithms

### experiment_31_1 
Rust code demonstrates a parallel, thread-safe Bloom Filter implementation in Rust for handling synthetic genomic data.

### experiment_31_2 
Rust code demonstrates a parallel, thread-safe MinHash implementation in Rust that computes approximate set similarities on synthetic genomic data.

### experiment_31_3 
Rust code demonstrates a scenario where different types of genomic data—such as microbial, eukaryotic, pangenome, single-cell transcriptomic, and Hi-C assays—require different indexing approaches.

### experiment_31_4 
Rust code demonstrates how to conditionally offload computations to GPU or FPGA in Rust by relying on feature flags in Cargo. 

### experiment_31_5 
Rust code simulates a distributed genomic indexing workflow. Each MPI rank builds a partial index from its assigned slice of the genome, storing data in a local HashMap.

### experiment_31_6 
Rust code reads each FASTA record, computes basic statistics such as GC content, and then outputs the aggregated results in JSON format.

### experiment_31_7
* Nextflow script declares two processes: one to compile the Rust code (optional if you prefer to compile manually) and another to run the resulting binary on the supplied FASTQ file. 
* Rust code reads FASTQ data using the needletail library, constructs a de Bruijn graph from the resulting k-mers, stores the same k-mers in a Bloom filter (leveraging Rayon for concurrency), and finally writes a minimal adjacency list to graph.json as well as the Bloom filter bits to bloom.json.

## 3.2. Sequence Data Structures and Strings Algorithms

### experiment_32_1
Rust code demonstrates a simple, parallel approach for detecting a known genomic “pattern” (e.g., a short sequence) within large genomic data using MPI and the Knuth-Morris-Pratt (KMP) algorithm.

### experiment_32_2
Rust code demonstrates how to build and merge partial suffix arrays for genomic data in an MPI environment by having each rank receive a segment of the reference (or read) text, construct a local suffix array, and then gather all partial results on rank 0 for global merging, enabling parallelization of suffix array construction—particularly important for large genomes where a single node might encounter time or memory constraints.

### experiment_32_3
* Nextflow script provides two processes: one for building the Rust binary (optional if you prefer a precompiled artifact), and another to run the partial suffix array program on a given FASTA file. 
* The Rust code processes a FASTA file by removing header lines, splitting the sequence into chunks, constructing a naive suffix array for each chunk in parallel, and serializing the results into a JSON file for further analysis.

## 3.3. Graph Data Structures for Genome Assembly and Beyond

### experiment_33
Rust code demonstrates how to build partial de Bruijn graphs from FASTQ data in Rust, leveraging needletail to efficiently parse sequences and Petgraph to represent overlapping k-mers in an undirected graph. 

## 3.4. Indexing and Searching in Large-Scale Biological Datasets

### experiment_34
Rust code illustrates a minimal Bloom filter implementation in Rust for k-mers derived from FASTA or FASTQ data.

## 3.5. High-Performance Computing and Parallelization

### experiment_35
Rust code demonstrates a robust approach to counting k-mers in large genomic data files (FASTA or FASTQ) using concurrent, local aggregation. It leverages the Needletail crate for reading sequence data, Rayon for parallel processing, and DashMap for merging partial results.
 
## 3.6. Putting It All Together—Rust and Nextflow Integration

### experiment_36_1
* Nextflow code illustrates how a Rust-based bioinformatics pipeline might be divided into multiple steps—building partial indexes from genomic data, performing alignments, and merging partial graphs—each handled by a dedicated module (simulating distinct crates).
* Rust code demonstrates a naive “FM-index building” step, in practice. Real HPC solutions for partial FM-index building would replace the placeholder with code that actually computes a partial Burrows–Wheeler Transform. 

### experiment_36_2
Nextflow code distributes read files across ephemeral containers during the ALIGN_READS step, merges outputs in SUMMARIZE, and allows for scalable HPC resource usage, modular testing with cargo test, and reproducibility through Docker container versioning.

### experiment_36_3
* Nextflow script (abyss_pipeline.nf) that demonstrates how you might chain together a Rust-based pre-processing utility and an ABySS assembly step.
* Rust program (rust_preprocess.rs) that demonstrates how you might parse command-line arguments (using the [structopt] or [clap] crate) and perform a rudimentary read-processing step.

