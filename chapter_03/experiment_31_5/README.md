## 3.1. Introduction to Data Structures and Algorithms

### experiment_31_5

This code demonstrates a parallel, thread-safe MinHash implementation in Rust that computes approximate set similarities on synthetic genomic data. The program generates random DNA sequences of a specified length and uses them to create two datasets. Each dataset is processed by a MinHasher that applies multiple hash seeds in parallel, determining a “minimum hash value” across all items for each seed. This approach is particularly relevant for high-throughput genomic or large-scale text data, where efficiently estimating set overlap becomes necessary.

Internally, the MinHasher maintains a specified number of hash seeds. For each item in a dataset, it computes a distinct hash based on the item and seed, retains only the smallest hash value for that seed, and repeats the process for all seeds. The resulting collection of smallest hash values, one per seed, becomes the MinHash signature. Two datasets can then be compared by examining their signatures: the fraction of identical positions in the signatures approximates the Jaccard similarity of the original sets. The code uses Rayon’s parallel iterators to distribute these computations across CPU cores, thereby increasing throughput in scenarios where large numbers of seeds or large datasets are involved.

#### Files contents:
* experiment_31_2/Cargo.toml (Cargo.toml file for dependencies)
* src/
  * main.rs (rust script)
  * output.txt (output file)

#### How to run:

run in powershell:

```powershell
cargo run | tee output.txt
```

(run main.rs and save the output in output.txt)
  
#### [dependencies]

```toml
mpi = "0.8"
serde = { version = "1.0", features = ["derive"] }
bincode = "2.0.0"
```

#### Explanation of the Output:
The program is a simulation of parallel computation using MPI (Message Passing Interface) to distribute the task of building and merging genomic indexes across multiple nodes. Let's walk through the key sections of the code and the resulting output:

##### Key Components:

###### 1. MPI Initialization:

* The MPI environment is initialized, and a world communicator is created.
* Each rank (or node) in the MPI world performs some computation and communicates with other ranks.

###### 2. Partial Index Creation:

* Each rank creates a local partial index for a slice of the genome. This index is represented by the PartialIndex struct, which contains:
  * rank_id: the rank of the node in the MPI world.
  * index_data: a HashMap<String, usize> where each rank generates a key-value pair with a unique string (key_rank_<rank_id>) and its corresponding rank_id as the value.
* Example: On rank 0, the local index data would be {"key_rank_0": 0}.

###### 3. Serialization and Data Gathering:

* The PartialIndex is serialized into a byte vector using bincode, and each rank sends its serialized chunk to the other ranks.
* The lengths of the serialized chunks are gathered using all_gather_into to ensure all ranks know the size of data from each other.
* The program calculates the total received buffer size and displacements (offsets where each rank’s data will be placed in the receive buffer).

###### 4. MPI All-Gather:

* Using all_gather_varcount_into, the serialized data from all ranks is collected into a single buffer on every rank.
* This function ensures that each rank receives the correct chunk from each other rank, based on their respective data lengths and displacements.

###### 5. Deserialization and Merging:

* After gathering the data, rank 0 processes the collected data:
* The serialized chunks are deserialized back into PartialIndex structs.
* These partial indexes are then merged into a global index, which is essentially a combination of all the index_data from each rank.

###### 6. Output:

* Finally, rank 0 prints the merged global index:
  * The merged index is expected to be a combination of all the partial index data from each rank.
* In this specific case, with only rank 0 running (since the output shows key_rank_0), the merged index is just:

```json
{"key_rank_0": 0}
```

* Reason for Output:
  * The output {"key_rank_0": 0} indicates that there is only one rank (rank 0) in this specific MPI run, or the other ranks are not contributing any data to the global index.
  * Rank 0 creates its partial index with the key key_rank_0 and value 0. Since this is the only rank involved, the global index that is merged only contains this single entry.

#### Conclusion:
* The program is correctly simulating parallel computation where multiple ranks would typically build partial indexes and merge them into a final global index. However, because the program seems to be running with a single rank (rank 0), the final global index contains only data from that rank.
* In a multi-rank run, if more ranks were involved, the merged global index would combine the data from all the ranks, leading to a more complex global index with multiple entries. The code is set up to handle this scenario, but in this case, only rank 0 is involved, resulting in the output showing just the partial index for that rank.

#### Key Takeaways:
* MPI communication is correctly set up for parallel computation, with serialization and deserialization to exchange data across nodes.
The output depends on the number of ranks used during the execution. In a real multi-rank environment, the merged index would contain more data from different ranks.



