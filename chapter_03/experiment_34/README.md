## 3.4. Indexing and Searching in Large-Scale Biological Datasets

### experiment_34

The code illustrates a minimal Bloom filter implementation in Rust for k-mers derived from FASTA or FASTQ data. It reads a specified input file, splits the reads into k-mers in parallel using Rayon, and then sequentially inserts those k-mers into a Bloom filter, which is ultimately serialized to JSON. By relying on parallel iteration for k-mer extraction, the code effectively distributes the workload across CPU cores, while avoiding data races or the need for locks by performing Bloom filter insertions in a single-threaded pass.

After reading in all sequences, each record is broken down into overlapping k-mers, which are gathered in a thread-safe manner using a parallel iterator. These k-mers are then inserted into the Bloom filter, which allocates a bit vector large enough for a given false-positive rate (determined by the number of bits and hash functions). The filter uses a simple FNV-like hashing approach, with multiple seeds ensuring that each k-mer sets several bits in the array. Finally, the completed Bloom filter can be queried for membership, and its bits are written to a JSON file for potential merging with partial filters from other nodes or processes.

In HPC practice, ephemeral containers can use this code to build partial Bloom filters from distinct input segments, which are then bitwise-ORed to yield a unified Bloom filter. This final structure can help detect whether a given k-mer appears in the entire dataset, thus enabling quick membership checks in large-scale read classification or reference indexing workflows. Carefully tuning the number of bits and number of hash functions is important: too few bits lead to high false-positive rates, while too many bits lead to excessive memory usage. In some industrial-scale applications, engineers split k-mers by certain prefixes before hashing, distributing the load across HPC nodes.

#### Files contents:
* experiment_34/
  * Cargo.toml (Cargo.toml file for dependencies)
* experiment_34/src/
  * main.rs (rust script)
  * reads.fq.rar (compressed reads.fq)
  * bloom.json
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
bitvec = "1.0.1"
```

#### Explanation of the Output
The Rust program implements a Bloom filter to store and check for the presence of k-mers (substrings of length 𝑘) in sequencing reads. The Bloom filter is a probabilistic data structure that efficiently checks set membership while allowing for false positives but no false negatives.

##### 1. output.txt
```rust
Bloom contains first k-mer 'ACGGAGGATGCGAGCGTTATCCGGATTTATT': true
Constructed Bloom filter (10650458 k-mers processed), result written to bloom.json
```
Breakdown of the Output
* The program successfully inserted 10,650,458 k-mers into the Bloom filter.
* The first extracted k-mer (ACGGAGGATGCGAGCGTTATCCGGATTTATT) is confirmed to exist in the Bloom filter.
* The Bloom filter is serialized and saved to bloom.json.

##### 2. bloom.json
This file contains the serialized Bloom filter, including:

* bits: The bit array storing hashed k-mer positions.
* num_bits: The total number of bits in the Bloom filter (10,000,000).
* num_hashes: The number of hash functions used (3).
* k: The k-mer length (31).

Example JSON Snippet
```json
{
  "bits": [4, 0, 48, 128, 10, 0, 0, 66, 0, 0, 16, 0, 8, 0, 1, ...],
  "num_bits": 10000000,
  "num_hashes": 3,
  "k": 31
}
```

##### Understanding the Data
* Step 1: Read FASTQ/FASTA Sequences
  * The program reads sequencing data from reads.fq.
  * It extracts all possible k-mers (31-mers) from the reads.
* Step 2: Initialize the Bloom Filter
  * The Bloom filter is initialized with:
  * 10,000,000 bits (approx. 1.25 MB of memory).
  * 3 hash functions (to reduce false positives).
* Step 3: Extract k-mers in Parallel
  * Uses rayon for parallel processing.
  * 10,650,458 k-mers were extracted.
* Step 4: Insert k-mers into the Bloom Filter
  * Each k-mer is hashed 3 times and inserted into the Bloom filter.
* Step 5: Check Membership
  * The program checks whether the first extracted k-mer exists in the Bloom filter.
  * Since it was just inserted, the result is true.
* Step 6: Serialize and Save Bloom Filter
  * The Bloom filter is saved to bloom.json for future use.

#### Conclusion
##### 1. Efficient k-mer Storage

* The Bloom filter stores 10+ million k-mers using only 10MB of space.
* This is much more memory-efficient than using a hash table.

##### 2. Probabilistic Nature

* The Bloom filter allows false positives (a k-mer might be reported as present when it isn't).
* However, false negatives are impossible (if a k-mer was inserted, it will always be found).

##### 3. Parallel Processing Boost

* The use of rayon significantly speeds up k-mer extraction.

##### 4. Scalability

* The Bloom filter can be used to efficiently filter and search large sequencing datasets.

##### 5. Next Steps

* The Bloom filter can be used in genome assembly to quickly check for sequencing errors or shared k-mers between datasets.

#### Final Thoughts
* The implementation successfully constructs a Bloom filter for k-mers.
* The parallelized k-mer extraction ensures high performance.
* The compact representation enables handling large genomic datasets efficiently.
