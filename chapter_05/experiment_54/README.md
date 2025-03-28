## 5.4. De Novo Assembly Approaches

### experiment_54

The following Rust program demonstrates a chunk-based strategy for building a k-mer count table and constructing a minimal de Bruijn graph from potentially large FASTQ inputs. By splitting the FASTQ reading process into chunks, it avoids loading all reads into memory at once. Each chunk’s records are then processed in parallel with Rayon, and the partial de Bruijn graphs are serialized to disk. In high-performance computing (HPC) or cloud-based environments, ephemeral tasks could run this chunking step in parallel for different portions of the input, writing separate partial graphs for downstream merging.

After parsing command-line arguments with clap, the program creates an output directory for storing partial de Bruijn graphs. The parse_fastx_file function from needletail is used to read the FASTQ records in a streaming manner, reading only a fixed number of records (chunk_size) per iteration. This buffer is then processed using Rayon’s .par_iter(), distributing the k-mer counting workload across available CPU cores.

During each chunk’s processing, the code constructs a local FnvHashMap of k-mers to their counts, then converts that into a small de Bruijn graph by treating each k-mer’s prefix as a node and the last base of the k-mer as an edge. This partial graph is serialized to disk using bincode. Once all chunks are exhausted, the code scans the partial outputs, reads each of them in turn, and merges them into a single de Bruijn graph. The final thresholding pass removes edges below a certain coverage value, creating a cleaner graph. This pattern of partial outputs followed by a merge step is a staple of HPC pipelines, allowing easy scaling to massive datasets while maintaining safe concurrency.

#### Files contents:
* experiment_54/
  * Cargo.toml (Cargo.toml file for dependencies)
*experiment_54/src/
  * main.rs (rust script)
  * reads.fq (fastq file)
  * reference.fa (fasta file)

#### How to run:

run in powershell:

```powershell
cargo run -- --input C:\Users\trian\BGVR\chapter_05\experiment_54\src\example.fastq
```

(run main.rs using two input dataset path)
  
#### [dependencies]

```toml
anyhow = "1.0"
rayon = "1.8"
needletail = "0.6.3"
fnv = "1.0"
serde = { version = "1.0", features = ["derive"] }
serde_json = "1.0"
bincode = "2.0.1"
clap = { version = "4.4", features = ["derive"] }
```

#### Explanation of the Output

##### Step 1: Reading the FASTQ file
* The program reads the FASTQ file in chunks (default: 10,000 sequences per chunk).

* Each sequence is processed to extract k-mers (default: k=31).

##### Step 2: Counting k-mers
* Each chunk is processed in parallel using rayon, where:

  * Every sequence is scanned for overlapping k-mers.

  * A hashmap stores the count of each k-mer.

* Only k-mers that meet a minimum threshold (default: 2 occurrences) are included in the graph.

##### Step 3: Constructing Partial de Bruijn Graphs
* The program builds a partial de Bruijn graph for each chunk and saves it to disk.

* Partial graphs are stored in binary format using bincode inside the directory partial_kmer_maps/.

##### Step 4: Merging Partial Graphs
* The program reads all partial graphs and merges them into a final graph.

* Any edges (k-1-mers → next base) with a count below the threshold are removed.

##### Step 5: Writing Final de Bruijn Graph
* The final de Bruijn graph is stored as a binary file (final_debruijn.bin).

##### Output
Assume we run the following command:

```powershell
cargo run -- --input example.fastq --k 5 --threshold 2 --chunk_size 5000
```

##### Output 

```rust
Processed chunk 0 with 5000 records, wrote partial de Bruijn to "partial_kmer_maps/partial_debruijn_0.bin"
Processed chunk 1 with 5000 records, wrote partial de Bruijn to "partial_kmer_maps/partial_debruijn_1.bin"
Processed chunk 2 with 4200 records, wrote partial de Bruijn to "partial_kmer_maps/partial_debruijn_2.bin"
Merging partial de Bruijn graphs...
Final de Bruijn graph has 1,234 prefix nodes. Written to "final_debruijn.bin".
```

* Generated Files
  * partial_kmer_maps/partial_debruijn_0.bin

  * partial_kmer_maps/partial_debruijn_1.bin

  * partial_kmer_maps/partial_debruijn_2.bin

  * final_debruijn.bin

#### Interpretation of the Output

* Processed chunk X with Y records
  → Each FASTQ chunk is processed, and a partial graph is saved.

* Final de Bruijn graph has N prefix nodes
  → The final graph contains N unique (k-1)-mer prefixes, meaning N unique k-mers were found.

* The binary file (final_debruijn.bin) contains the merged and filtered de Bruijn graph, which can be later deserialized for further analysis.

#### Conclusion
* The program successfully constructs a de Bruijn graph in a memory-efficient and parallelized manner.
* The chunking mechanism allows processing large FASTQ files without high memory usage.
* The final graph can be used for genome assembly, sequence error correction, or bioinformatics analysis.

