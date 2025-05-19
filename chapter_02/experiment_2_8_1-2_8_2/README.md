## 1.8. Tools and Frameworks

### experiment_18_1_2 

#### Rust Code

The companion Rust code illustrates how to Read FASTQ files, parse out the raw sequencing reads, and then build a simple De Bruijn graph.

#### Nextflow Code
The following Nextflow pipeline showcases a minimal yet practical approach to downloading raw FASTQ data from NCBI using SRAToolkit and then invoking a Rust-based program to construct a De Bruijn graph.

The following Nextflow pipeline showcases a minimal yet practical approach to downloading raw FASTQ data from an online source and then invoking a Rust-based program to construct a De Bruijn graph. Written in the DSL2 syntax, this pipeline highlights how processes can be modularized—one process handles data retrieval, while the other executes Rust code to transform reads into a graph structure. By relying on containerization and Nextflow’s dataflow paradigm, the entire workflow remains reproducible and amenable to HPC or cloud environments, ensuring that each step can scale seamlessly.

The pipeline defines two processes: DOWNLOAD_FASTQ fetches the raw sequencing data (e.g., via curl), and BUILD_DEBRUIJN invokes a compiled Rust binary that performs De Bruijn graph construction. Channels orchestrate data movement: the FASTQ URL is transformed into a pipeline artifact that feeds into the download task, whose output is then passed to the Rust process. Each process can be containerized (e.g., using Docker or Singularity images), thereby pinning software versions and system dependencies to fixed states. In large-scale deployments, Nextflow’s native compatibility with HPC schedulers or cloud engines provides dynamic resource allocation, ensuring the pipeline can handle everything from small test sets to multi-terabyte genomic datasets.

The companion Rust code illustrates how to read FASTQ files, parse out the raw sequencing reads, and then build a simple De Bruijn graph. A De Bruijn graph represents k-mer overlaps as edges, forming a powerful framework for tasks like genome assembly or read error correction. This demonstration relies on basic I/O mechanisms in Rust, along with simple data structures (e.g., a hash map) to store adjacency relationships between consecutive k-mers. You might name this file build_de_bruijn.rs, compile it inside your container, and place the resulting executable in a suitable path.

Inside the Rust program, command-line arguments specify a k-mer size, an input FASTQ path, and an output destination. After minimal parsing of FASTQ records—here, implemented by reading lines in groups of four—the code constructs each node by extracting consecutive k-length substrings from the reads. By inserting edges from each k-mer to its successive overlap, the hash map encodes the resultant De Bruijn graph, which is then written to a file. This stand-alone Rust executable is easily integrated into container images and invoked by Nextflow in a larger pipeline, leveraging concurrency and memory safety to handle substantial genomic data at scale.

#### Project Structure:

```plaintext
experiment_18_1_2/
├── Cargo.toml                     # Rust project configuration and dependencies
└── src/
    ├── main.rs                    # Main Rust script containing program logic
    ├── main.nf                    # Nextflow workflow script
    ├── de_bruijn_graph.rar        # Compressed de Bruijn graph text output file
    ├── SRR11192680_1.rar          # Compressed SRR11192680_1.fastq
    ├── SRR11192680_2.rar          # Compressed SRR11192680_2.fastq
    └── reads.rar                  # Compressed reads.fastq (renamed from SRR11192680_1.fastq)
```

#### How to run:

```wsl
nextflow run main.nf
```

(run the nextflow script that will run the workflow and save the output in de_bruijn_graph.txt)
  
#### [dependencies]

only use standard "std"

#### Explanation of the Output: de_bruijn_graph.txt

The output file represents a De Bruijn graph constructed from the sequencing data obtained from an SRA (Sequence Read Archive) file. The De Bruijn graph is a commonly used data structure in bioinformatics, particularly for genome assembly. Below is an explanation of the key components of the output:

##### 1. Header Section

The first lines of the file provide metadata about the constructed graph:

De Bruijn Graph (k=21)

Number of nodes: 290694

* De Bruijn Graph (k=21): The k-mer size used in the construction of the graph is 21 nucleotides. This means that each node in the graph represents a 21-mer (a sequence of 21 bases).
* Number of nodes: 290694: The total number of unique k-mers (nodes) present in the dataset.

##### 2. Node-Edge Representation

Each line in the file represents a node (a k-mer) and its outgoing edges (connections to other k-mers):

  GTCCGGTGTGAAAGTCTATCG => ["TCCGGTGTGAAAGTCTATCGC"]

  TTAAGCCAGTGGGGAAAGTTT => ["TAAGCCAGTGGGGAAAGTTTG"]

* The left-hand side of the => represents a node (a 21-mer).
* The right-hand side (inside ["..."]) represents the successor k-mer(s), which is obtained by shifting one nucleotide to the right.

  For example:

  GTCCGGTGTGAAAGTCTATCG => ["TCCGGTGTGAAAGTCTATCGC"]
  
* The node "GTCCGGTGTGAAAGTCTATCG" is a 21-mer.
* Its successor "TCCGGTGTGAAAGTCTATCGC" is formed by shifting one base forward.
  
##### 3. Multiple Edges

Some nodes have multiple outgoing edges, indicating alternative paths in the sequencing data:

  CAGGGGCTCAACCCCGGTACT => ["AGGGGCTCAACCCCGGTACTG", "AGGGGCTCAACCCCGGTACTG", "AGGGGCTCAACCCCGGTACTG",     "AGGGGCTCAACCCCGGTACTG", "AGGGGCTCAACCCCGGTACTG"]
  
* The k-mer "CAGGGGCTCAACCCCGGTACT" connects to multiple instances of "AGGGGCTCAACCCCGGTACTG", suggesting repeated occurrences in the sequencing data.
  
##### 4. Highly Repetitive K-mers

Some nodes are connected to a large number of identical successor k-mers:


  AATTCCCGGTGTAGCGGTGGA => ["ATTCCCGGTGTAGCGGTGGAA", "ATTCCCGGTGTAGCGGTGGAA", ..., "ATTCCCGGTGTAGCGGTGGAA"]

* The high frequency of this k-mer suggests that it may belong to a repetitive or conserved region in the genome.
  
##### 5. Terminal Nodes

Nodes without further connections represent the end of certain sequence paths:

  ACTGCTTTTGAAACTGCCAGA => ["CTGCTTTTGAAACTGCCAGAC"]

* These nodes are likely found at the end of sequencing reads or fragmented regions.

#### Conclusion

The program successfully constructed a De Bruijn graph from sequencing data by extracting overlapping k-mers (k=21) and mapping their connections. This graph representation captures the local structure of the genome, enabling downstream applications such as genome assembly, variant detection, and error correction in sequencing data analysis.

