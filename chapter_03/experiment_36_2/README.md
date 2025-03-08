## 3.6. Putting It All Together—Rust and Nextflow Integration

### experiment_36_2

#### 1. Nextflow

Similarly, the alignment crate might accept FASTQ reads and the partial indexes, using a local or global alignment algorithm. HPC tasks parse multiple read files in parallel with rayon or distribute them across ephemeral containers. The summarizer merges alignment outputs by reading JSON files from multiple ephemeral tasks. If scaling to thousands of files or extremely large references, code must handle concurrency overhead and partial merges carefully.

This pipeline “scatters” read files across ephemeral containers in the ALIGN_READS step, then merges partial outputs in SUMMARIZE. Each ephemeral container runs one crate. Users can quickly scale HPC resources on local clusters or in the cloud, launching dozens or hundreds of ephemeral containers as needed. The modularity fosters continuous integration: each Rust crate is tested individually with cargo test, and the containers can be versioned in a Docker registry, guaranteeing reproducible runs.



#### Files contents:
* experiment_36_2/
  * Cargo.toml (Cargo.toml file for dependencies)
* experiment_36_1//src/
  * main.rs (rust script)
  * main.nf (nextflow script)
  * reference.fa (reference fasta file)

#### How to run:

run in powershell:

```powershell
cargo run main.nf | tee output.txt
```

(run main.nf and save the output in output.txt)
  
#### [dependencies]

```toml
serde = { version = "1", features = ["derive"] }
serde_json = "1"
bincode = "2.0.0"
```

#### Explanation of the Output

##### 1. Output from main.rs
The main.rs script processes a reference genome file (reference.fa), chunks it into segments of a specified size (default: 1,000,000 characters), and creates partial FM-index placeholders for each chunk. The output indicates:

* Message in output.txt:
```rust
Created partial_fm_indexes.json with 1 chunks.
```

This means that only one chunk was created, which suggests that the input reference file was smaller than the defined chunk size.

* partial_fm_indexes.json Content:

```json
[
  {
    "reference_chunk": "chunk_0",
    "bwt_data": [47, 47, 32, 84, 104, 105, 115, 32, 101, 120, 97, 109, ...]
  }
]
```

This JSON file contains:

* "reference_chunk": "chunk_0" → Indicates the name of the chunk (first chunk).
* "bwt_data": [...] → Stores the raw bytes of that chunk. This data would typically be used for suffix arrays or BWT transformations in indexing.

##### 2. Output from main.nf
The main.nf script simulates an HPC pipeline where multiple processing steps happen in parallel, including:

* Indexing (rust_index_builder) → Creates a partial index from reference regions.
* Alignment (rust_aligner) → Simulates read alignment.
* Graph Assembly (rust_graph_assembler) → Builds graphs from indexes and alignment scores.

The output includes:

```rust
Building partial graph for region: chr1_partA
Found alignment result: query_id=readA, score=10
Found alignment result: query_id=readB, score=10
Building partial graph for region: chr1_partB
Found alignment result: query_id=readC, score=10
Found alignment result: query_id=readD, score=10
Global Graph: PartialGraph { node_count: 72 }
```

This means:

* Two partial graphs were built (for chr1_partA and chr1_partB).
* Each read gets a score based on length.
* The final merged graph has 72 nodes.

#### Conclusion

##### 1. main.rs (Partial FM-Index Generation)

* Successfully creates partial FM-index placeholders for reference genome segments.
* The reference file was small, so only one chunk was generated.
* This method could be used in large-scale genomic indexing for fast searching.

##### 2. main.nf (HPC-Style Bioinformatics Pipeline)

* Simulates an HPC workflow, where:
  * Indexing builds local indexes.
  * Alignment assigns scores to reads.
  * Graph assembly merges results.
* Successfully merges data from multiple regions into a single final graph.

#### Key Takeaways
* The pipeline demonstrates how bioinformatics workflows scale using modular Rust components.
* The FM-indexing approach allows efficient reference genome compression and searching.
* The pipeline could be expanded with real bioinformatics algorithms for actual sequencing data analysis.

