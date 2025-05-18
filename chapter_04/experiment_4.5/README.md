## 4.5. Transcriptomics and Alternative Splicing Algorithms

### experiment_4.5

This Rust code snippet demonstrates a straightforward approach to building and merging partial splicing graphs, which represent gene transcripts by connecting exons and splice junctions. In typical RNA-seq or transcriptomics workflows, large numbers of read alignments must be processed to infer splicing patterns across the genome. By chunking these alignments, processing them in parallel, and merging partial results, this solution efficiently scales to large datasets commonly encountered in high-throughput sequencing experiments.

The code defines three primary data structures. First, ExonSegment represents a genomic segment—often used to store exon boundaries. In this example, it is kept minimal, with _start and _end fields underscored to indicate they are currently unused. Second, Alignment holds information about individual read alignments, such as the chromosome name and start/end positions. Third, SplicingGraph is a map-based structure that stores connections (or “junctions”) between exonic segments and tracks coverage or counts associated with those connections.

To handle large datasets efficiently, the input alignments are divided into chunks of a fixed size (chunk_size). Using Rust’s .chunks(chunk_size) method, the code produces multiple smaller vectors of alignments from the overall set. These chunks are processed in parallel using Rayon’s .into_par_iter(), which automatically distributes the workload across available CPU cores. Each chunk is passed to process_alignment_chunk, where a local SplicingGraph is built by iterating over each alignment and adding the corresponding exon-junction data.

Instead of locking a single global graph structure during processing, each chunk produces its own local splicing graph, capturing exons and junctions present in that subset of data. After all chunks are processed, Rayon’s .reduce(...) method merges each local graph into one final SplicingGraph. The merge function in SplicingGraph unifies adjacency data by extending any overlapping exon connections, effectively combining partial coverage information from all chunks into a single, coherent splicing graph.

Once all chunks have been processed and merged, the final splicing graph is written to disk. In a production scenario, this step would typically involve binary serialization (e.g., using bincode) or a more advanced format (e.g., JSON) to store adjacency maps and coverage information. In this demonstration, however, a placeholder text output is used to illustrate that the consolidated graph is readily available for downstream analysis, visualization, or further data integration.

#### Project Structure:

```plaintext
experiment_4.5/
├── Cargo.toml                        # Rust project configuration and dependencies
└── src/
    ├── main.rs                       # Main Rust script containing program logic
    ├── partial_splicing_graph.bin    # Binary format output file containing partial splicing graph
    └── output.txt                    # Text output file
```

#### Cargo.toml

```toml
[package]
name = "experiment_4.5"
version = "0.1.0"
edition = "2024"

[dependencies]
rayon = "1.10.0"
```

#### How to run:

run main.rs in powershell:

```powershell
cargo run | tee output.txt
```
(run main.rs and get the partial_splicing_graph.bin output and output.txt)


#### Explanation of the Output
My Rust program constructs a splicing graph based on a set of alignment data. The output consists of:

##### 1. Console Message (output.txt)

```rust
Splicing graph has been successfully written.
```

This confirms that the program executed successfully and generated the splicing graph.

##### 2. Binary File (partial_splicing_graph.bin)

```rust
Serialized splicing graph example
SplicingGraph {
    adjacency: {
        (100, 200): [(ExonSegment { _start: 150, _end: 300 }, 1)],
        (150, 300): [(ExonSegment { _start: 200, _end: 400 }, 1)],
        (500, 700): [(ExonSegment { _start: 550, _end: 800 }, 1)],
    },
}
```

This is a serialized representation of the splicing graph, showing how exons are connected.

#### Interpretation of the Splicing Graph
* Nodes: Each exon is represented by (start, end).
* Edges: A connection exists between exons if there is a splicing event.
* Example Connections:
  * (100, 200) → (150, 300)
  * (150, 300) → (200, 400)
  * (500, 700) → (550, 800)

This structure represents alternative splicing, where exons can be skipped, included, or linked in different ways.

#### Conclusion
The program successfully builds a splicing graph from alignment data:

* Uses parallel processing (rayon) to improve performance.
* Correctly links exon segments based on alignment data.
* The splicing graph is stored and serialized for further analysis.
