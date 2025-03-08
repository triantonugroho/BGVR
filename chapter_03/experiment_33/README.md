## 3.3. Graph Data Structures for Genome Assembly and Beyond

### experiment_33

The following code demonstrates how to build partial de Bruijn graphs from FASTQ data in Rust, leveraging needletail to efficiently parse sequences and Petgraph to represent overlapping k-mers in an undirected graph. By splitting the reads into chunks and processing them in parallel with Rayon, each segment of data contributes to a partial de Bruijn graph, which is then serialized for further merging or analysis. In practice, this design can handle large genomic datasets by distributing work across multiple cores and managing memory usage more effectively than a single, monolithic graph construction.

The program reads all sequences from a FASTQ file, divides them into equal-sized chunks, and maps them to individual threads. Each chunk is processed by extracting every k-mer from each read, creating or reusing a node in the graph, and adding edges between consecutive k-mers. After building the local subgraph, the code collects the node labels and edges in a data structure, serializes it, and writes it out as a JSON file. This facilitates modularity by allowing intermediate or partial results to be saved and later combined into a final, global de Bruijn graph if desired.

#### Files contents:
* experiment_33/
  * Cargo.toml (Cargo.toml file for dependencies)
* experiment_33/src/
  * main.rs (rust script)
  * reads.fq.rar (compressed reads.fq)
  * partial_debrujin_graphs.json.rar (compressed partial_debrujin_graphs.json output file)
  * output.txt (output file)

#### How to run:

run in powershell:

```powershell
cargo run | tee output.txt
```

(run main.rs and save the output in output.txt)
  
#### [dependencies]

```toml
rayon = "1.7"
needletail = "0.6"
serde = { version = "1", features = ["derive"] }
serde_json = "1"
petgraph = "0.7.1"
```

#### Explanation of the Output
The Rust program implements a parallelized de Bruijn graph construction for sequencing reads using k-mers (substrings of length  𝑘. The program reads a FASTA/FASTQ file, extracts sequences, splits them into chunks, and constructs partial de Bruijn graphs in parallel using the rayon library.

##### 1. output.txt
```rust
Wrote partial de Bruijn graphs to partial_debruijn_graphs.json
```

* This message confirms that the de Bruijn graph construction completed successfully and the results were saved to partial_debruijn_graphs.json.

##### 2. partial_debruijn_graphs.json
This JSON file contains the partial de Bruijn graphs generated from the sequencing reads. The structure includes:

* nodes: A list of unique k-mers (31-mers in this case).
* edges: A list of tuples (i, j) indicating an edge between k-mer nodes[i] and nodes[j].
* k: The length of each k-mer (set to 31 in this case).
* Example JSON Snippet
```json
{
  "nodes": [
    "ACGGAGGATGCGAGCGTTATCCGGATTTATT",
    "CGGAGGATGCGAGCGTTATCCGGATTTATTG",
    "GGAGGATGCGAGCGTTATCCGGATTTATTGG",
    "GAGGATGCGAGCGTTATCCGGATTTATTGGG",
    "AGGATGCGAGCGTTATCCGGATTTATTGGGT",
    "GGATGCGAGCGTTATCCGGATTTATTGGGTT"
  ],
  "edges": [
    [0, 1],
    [1, 2],
    [2, 3],
    [3, 4],
    [4, 5]
  ],
  "k": 31
}
```

Understanding the Data
* Each node represents a unique k-mer (substring of length 31).
* Each edge represents a connection between consecutive k-mers (i.e., overlapping by 𝑘−1 characters).
* The graph structure enables genome assembly by linking overlapping k-mers.

#### 3. Explanation of the Algorithm
* Step 1: Read Input FASTQ/FASTA File
  * The program reads sequencing reads from reads.fq.
  * Uses needletail library to parse the FASTA/FASTQ format efficiently.
* Step 2: Partition Reads into Chunks
  * The sequencing reads are split into smaller chunks (100,000 reads per chunk).
  * This enables parallel graph construction.
* Step 3: Construct Partial de Bruijn Graphs
  * For each chunk, the program:
  * Extracts k-mers (substrings of length k = 31).
  * Adds k-mers as nodes to the graph.
  * Creates edges between consecutive k-mers.
  * Stores the resulting graph.
* Step 4: Parallel Processing
  * Uses rayon to process multiple chunks in parallel.
  * Each chunk produces a partial de Bruijn graph.
* Step 5: Output Partial Graphs
  * Each partial graph is serialized into JSON (partial_debruijn_graphs.json).
  * This allows further processing (e.g., merging into a full de Bruijn graph).

#### Conclusion
##### 1. Successful Construction of Partial de Bruijn Graphs

* The output JSON file confirms that the program successfully generated partial de Bruijn graphs.

##### 2. Efficient Parallelization

* The use of rayon allows for parallel processing, improving performance for large datasets.

##### 3. Scalability

* The chunk-based approach makes it possible to handle massive sequencing datasets without running out of memory.

##### 4. Graph Representation is Correct

* The edges properly reflect overlapping k-mers, ensuring the graph structure correctly represents sequence connectivity.

##### 5. Next Steps

* The partial graphs can be merged to form a global de Bruijn graph.
* Can be used for genome assembly, error correction, or variant detection.

##### Final Thoughts
* This pipeline efficiently constructs partial de Bruijn graphs from sequencing reads.
* The parallel execution improves performance and scales well for large sequencing datasets.
* The graph-based approach is useful for genome assembly and sequence analysis.
