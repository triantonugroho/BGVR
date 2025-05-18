## 4.6. Single-Cell Functional Genomics

### experiment_46_1

This Rust code demonstrates a parallelized approach to constructing a k-nearest neighbor (k-NN) graph for a set of single-cell data, a common task in bioinformatics pipelines where large datasets must be efficiently processed. After building the k-NN graph, a simple breadth-first search (BFS) is used to derive a “pseudotime” measure, which can provide a rough chronological ordering of cells based on their similarity.

The main data structure is a Cell, representing an individual sample’s gene expression profile as a vector of floating-point values. Each Cell is assigned an unused _id field (prefixed with an underscore to avoid compiler warnings) for potential referencing or labeling. We also define a KnnGraph struct to hold edges, where edges[i] gives the indices of the k nearest neighbors to cell i. This adjacency-based graph representation allows us to quickly traverse connections among cells once the graph is constructed. Before building the graph, the code creates a small synthetic dataset of cells. In a real-world scenario, one could load thousands or millions of cells from a file or a streaming interface, which makes the parallel design critical for scalability.

The core function, build_knn_graph, leverages Rayon to run in parallel. First, the total number of cells n is determined. An Arc<Mutex<_>> protects a shared vector of adjacency lists. Each index i corresponds to one cell, and in a parallel loop, the code calculates the distance from cell i to every other cell j. These distances are stored in a local vector, sorted, and truncated to the top k nearest neighbors. Writing those neighbors back into the shared structure requires briefly locking a mutex to avoid data races. For smaller datasets, the straightforward O(n^2) distance computation is acceptable. For larger datasets, approximate nearest neighbor methods or chunk-based merging can drastically improve performance. Regardless of the method used to find neighbors, the resulting data structure is a clear and efficient representation of cell-to-cell adjacency.

Once the k-NN graph is available, the compute_pseudotime function uses a breadth-first search to measure how many “steps” each cell is from a chosen root cell. This BFS distance—stored in the pseudotime vector—serves as a basic proxy for progression, such as in single-cell developmental trajectories. The BFS continues until all reachable nodes are processed, resulting in a vector that can be used for downstream tasks like ordering or clustering cells based on their distances from the root.

#### Project Structure:

```plaintext
experiment_4.6_1/
├── Cargo.toml                     # Rust project configuration and dependencies
└── src/
    ├── main.rs                    # Main Rust script containing program logic
    └── output.txt                 # Text output file
```

#### Cargo.toml

```toml
[package]
name = "experiment_4.6_1"
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
(run main.rs and save the output in output.txt)


#### Explanation of the Output
My Rust program constructs a k-Nearest Neighbor (k-NN) graph based on the expression profiles of single-cell data and then computes a pseudotime trajectory using Breadth-First Search (BFS).

##### Output Breakdown

###### 1. k-NN Graph Construction
The program computes pairwise Euclidean distances between cells and selects the top-2 nearest neighbors (k = 2) for each cell. The results:

```rust
k-NN edges:
Cell 0 neighbors = [4, 1]
Cell 1 neighbors = [4, 0]
Cell 2 neighbors = [3, 1]
Cell 3 neighbors = [2, 1]
Cell 4 neighbors = [0, 1]
```

* Each cell connects to two nearest neighbors based on expression similarity.
* Example:
  * Cell 0 is closest to Cell 4 and Cell 1.
  * Cell 2 is closest to Cell 3 and Cell 1.

###### 2. Pseudotime Computation
The program computes pseudotime from Cell 0 (used as the root node). Pseudotime is calculated as BFS distance from this root.

```rust
Pseudotime from root=0: [0.0, 1.0, inf, inf, 1.0]
```

* Cell 0 is the root, so its pseudotime is 0.0.
* Cell 1 and Cell 4 are direct neighbors of Cell 0, so they have a pseudotime of 1.0.
* Cell 2 and Cell 3 are not connected to Cell 0, so they have infinite (inf) pseudotime (unreachable from root).

#### Conclusion
The program successfully constructs a k-NN graph and computes pseudotime using BFS.

* The k-NN structure correctly identifies local neighborhoods based on expression similarity.
* The BFS-based pseudotime metric captures relative progression from the root cell.
* Disconnected cells (Cell 2, Cell 3) result in infinite pseudotime, indicating that k=2 may be too small to fully connect the dataset.

