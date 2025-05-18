## 3.1. Introduction to Data Structures and Algorithms

### experiment_3.1_7

#### 1. Nextflow
An example Nextflow pipeline and accompanying Rust code that demonstrate how to construct both a de Bruijn graph and a Bloom filter for genomic analysis on FASTQ data. This setup is intended to be standalone on a typical PC, using Nextflow to orchestrate the workflow locally, while the Rust code handles reading the FASTQ file, building the data structures, and parallelizing the workload via Rayon. In practice, you would adapt these skeletons for large-scale data and HPC environments, possibly adding chunked data handling, distributed file systems, or MPI-based communication. Nonetheless, the example shows how these pieces fit together in principle.

Create a file named main.nf in a project directory. This script declares two processes: one to compile the Rust code (optional if you prefer to compile manually) and another to run the resulting binary on the supplied FASTQ file. You can run it locally with the Nextflow CLI (nextflow run main.nf --fastq reads.fastq).

#### 2. Rust
Create a Cargo project by running a command like cargo new debruijn_bloom, then replace the src/main.rs file with the sample code provided. In this example, the program reads FASTQ data using the needletail library, constructs a de Bruijn graph from the resulting k-mers, stores the same k-mers in a Bloom filter (leveraging Rayon for concurrency), and finally writes a minimal adjacency list to graph.json as well as the Bloom filter bits to bloom.json. In a real production setting, you would likely expand on these data structures and serialization methods to address larger or more specialized genomic workflows.

To run the pipeline, first install Nextflow on your PC, for example by running curl -s https://get.nextflow.io | bash. Also confirm you have Rust and Cargo installed, typically via rustup. Next, create the Nextflow script main.nf and set up the Rust project with the corresponding Cargo.toml and src/main.rs. Finally, in your terminal, execute nextflow run main.nf --fastq reads.fastq. This command will compile the Rust code (the compile process), then run debruijn_bloom on reads.fastq (the analysis process), writing graph.json and bloom.json to the results directory.

This example is highly simplified. Real genomic use cases often involve splitting large FASTQ files, merging partial data structures across multiple nodes, and handling error correction or advanced graph algorithms. Nevertheless, the above sample shows how Nextflow pipelines can orchestrate Rust HPC tasks—like constructing de Bruijn graphs and Bloom filters—all on a local PC without external cluster dependencies.

#### Project Structure:

```plaintext
experiment_3.1_7/
├── Cargo.toml                     # Rust project configuration and dependencies
└── src/
    ├── main.rs                    # Main Rust script containing program logic
    ├── main.nf                    # Nextflow workflow script
    ├── example.fastq.rar          # Compressed FASTQ example file
    ├── output.txt                 # Text output file
    ├── results/
    │   ├── bloom.json.rar         # Compressed Bloom filter JSON results
    │   └── graph.json.rar         # Compressed graph JSON results
    └── target/
        └── rustc_info.json        # Rust compiler information file
```

#### Cargo.toml
```toml
[package]
name = "experiment_3.1_7"
version = "0.1.0"
edition = "2021"

[dependencies]
rayon = "1.7"
needletail = "0.6"
serde = { version = "1", features = ["derive"] }
serde_json = "1"
```

#### How to run:

run in powershell:

```powershell
cargo run main.nf | tee output.txt
```

(run main.nf  and save the output in output.txt)


#### Explanation of the Output:
This system is a Nextflow pipeline that compiles and runs a Rust program to analyze FASTQ data, building a De Bruijn graph and a Bloom filter from DNA sequences. Let's break down the steps and explain the output:

##### 1. Nextflow Pipeline (main.nf):

* compile process: The compile process is responsible for compiling the Rust program. It ensures that the Rust code is compiled and available as an executable (debruijn_bloom).
* analysis process: Once the binary is compiled, the analysis process takes over. It runs the compiled Rust program, providing it with a FASTQ file and other necessary parameters. The Rust program processes the data, building a De Bruijn graph and a Bloom filter from the k-mers (substrings of length k in the sequences).

##### 2. Rust Program (main.rs):

* The Rust program:
  * Reads the FASTQ file: It uses the needletail crate to read the FASTQ file and stores the sequences.
  * Builds a De Bruijn graph: It constructs a graph where each node is a k-mer, and edges represent connections between k-mers in the sequence.
  * Creates a Bloom filter: It inserts k-mers into a Bloom filter, a space-efficient probabilistic data structure that tells whether an element is present or likely absent.
  * Outputs the results: It saves the De Bruijn graph and Bloom filter as JSON files (graph.json and bloom.json).

##### 3. Parameters:

* --fastq: The input FASTQ file (e.g., example.fastq).
* --kmer: The size of k-mers (default is 31).
* --outdir: The output directory for the results (default is results).

##### Expected Output:
Assuming the pipeline runs successfully, the output will consist of:

###### 1. Console Output: The terminal will show messages indicating the pipeline's progress:

* "Compiling the Rust code..."
* "De Bruijn & Bloom" (indicating the analysis process is running)
* "Reading from FASTQ: [file], k-mer=[size]" (showing the FASTQ file and k-mer size used)
* "Loaded [number of reads] reads. Building de Bruijn graph & Bloom filter..."
* "De Bruijn graph saved to [path]/graph.json"
* "Bloom filter saved to [path]/bloom.json"

###### 2. Files Created: The analysis will produce the following JSON files in the specified output directory (results by default):

* graph.json: Contains the De Bruijn graph, with each k-mer as a node and edges representing the connections between k-mers in the sequences.
* bloom.json: Contains the Bloom filter, representing the set of all distinct k-mers.

##### Explanation of the graph.json and bloom.json files:

###### 1. graph.json:

* This file represents the De Bruijn graph. Each node corresponds to a unique k-mer found in the FASTQ sequences.
* The graph structure can be used to understand how k-mers are connected (which k-mers tend to follow one another in the sequences).
* The format could be a list of nodes, where each node contains the k-mer string and a set of edges representing other k-mers that it connects to.

###### 2. bloom.json:

* This file represents the Bloom filter used to quickly check for the existence of k-mers in the data.
* It’s a probabilistic data structure with a set of bits that are set to true based on the hash values of the inserted k-mers.
* It can tell if a k-mer is "probably" in the data (with some false positives) but guarantees that it’s definitely not present if the result is false.

#### Conclusion:
* Purpose of the Pipeline: The pipeline compiles and runs a Rust program that reads FASTQ data, builds a De Bruijn graph, and generates a Bloom filter.
* Efficiency: The Rust program utilizes parallelism (via the rayon crate) to efficiently build the Bloom filter and process the data, especially useful for large-scale genomic data analysis.
* Output Data: The output files (graph.json and bloom.json) provide the results of the analysis:
* The De Bruijn graph (stored in graph.json) represents the structure of k-mers in the sequences.
* The Bloom filter (stored in bloom.json) is a space-efficient way to check if k-mers exist in the data.
This pipeline is a useful tool for bioinformatics tasks, particularly for DNA sequence assembly and analysis, by leveraging the power of Rust for performance and Nextflow for workflow management.



