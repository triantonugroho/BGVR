# 1.3. De Bruijn Graphs

Below is a Rust code that demonstrates building a simple De Bruijn graph from a FASTA file by extracting overlapping k-mers and linking them, thereby avoiding pairwise alignments for each read. De Bruijn graphs are commonly employed in genome assembly because they can efficiently represent overlapping k-mers, simplifying reconstruction. By keeping memory usage explicit and parallelizing iteration, the program remains robust and scalable to large data sets containing millions of reads. It leverages Rust’s strong memory safety and concurrency features, along with crates like “rayon” for parallel processing and “nalgebra”/“ndarray” for numeric tasks, ensuring efficient performance in high-throughput sequencing scenarios.

This program imports crates for bioinformatics (bio), linear algebra (nalgebra, ndarray), and parallel computation (rayon). In the build_de_bruijn function, each sequence is processed by extracting overlapping k-mers of length k (referred to as node), along with the subsequent overlapping k-mer (edge). Rather than modifying a single shared hash map in multiple threads, each thread accumulates its own local map of k-mer relationships. Rayon’s .par_iter() enables concurrent iteration over the sequences, and a reduce operation merges these local maps into a single global HashMap<String, Vec<String>>.

After reading a FASTA file named reads.fasta from the src directory, the main function collects all sequences into a vector and invokes build_de_bruijn with a chosen k-mer size (k = 21). Once the De Bruijn graph is built, the code demonstrates HPC-oriented crates by creating example matrices with “nalgebra” and “ndarray.” Finally, it prints out details about the graph’s size and the matrices, confirming that the parallel construction and supporting data structures have been set up correctly.

Several success stories highlight substantial reductions in runtime and cost once legacy Python or Java components are rewritten in Rust, particularly for k-mer counting or parallel motif searches. By combining HPC scheduling with containerized Rust executables, these organizations accelerate the pace of biomarker discovery and gene therapy research while preserving the reproducibility needed for regulatory compliance and collaborative studies.

## Files contents:
* main.rs (rust script)
* main.nf (nextflow script)
* reads.fasta (fasta file)

## How to run:

cargo run main.nf (run the nextflow script that will run the main.rs and save the output in output.txt)
  
## [dependencies]

bio = "2.0.3"

nalgebra = "0.33.2"

ndarray = "0.16.1"

rayon = "1.10.0"

