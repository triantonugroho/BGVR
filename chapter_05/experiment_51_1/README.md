## 5.1. Introduction to Sequence Analysis

### experiment_51_1

The following code showcases a HyperLogLog (HLL) implementation in Rust, a probabilistic data structure designed for estimating the cardinality (the count of unique elements) of a data set. Unlike traditional methods that store every distinct element, HyperLogLog employs a fixed-size array of small counters, making it extremely memory efficient, particularly for very large data streams. The structure is defined by a parameter p, which determines the number of registers m=2^p. Each incoming item is hashed to a 64-bit value, and the top p bits index into one of these registers. The remaining bits are then examined to find the longest run of leading zeros, and the register is updated with the maximum such run observed over time. By only storing these small counters, HyperLogLog achieves a compact representation of a potentially massive set.

A key feature of HyperLogLog is its mergeability. Two HLL instances of the same precision can be combined by taking the maximum value of each corresponding register, reflecting the union of elements from both data streams. This allows distributed or parallel systems to maintain separate sketches for subsets of data and aggregate them into a single unified estimate at any time. Another strength is the ability to handle any data type implementing Rust’s Hash trait, from simple strings to complex objects. Although HyperLogLog provides an approximate count rather than an exact value, its relative error typically remains low, and it can be tuned by adjusting the precision parameter p to trade off memory usage against accuracy.

The HyperLogLog structure is defined with a precision field p, a vector of 8-bit registers, and a constant alpha used in the cardinality estimation formula. During construction in new, the registers are initialized to zero, and alpha is set based on the number of registers. The add method takes an element, hashes it to 64 bits, and uses the top p bits to select a register. It shifts away these top bits and counts the leading zeros in the remainder to obtain a run-of-zeros value. If this value exceeds the register’s current value, the register is updated to reflect the newly observed maximum.

The estimation of unique elements happens in the estimate function, which applies the canonical HyperLogLog formula. The formula multiplies a bias-correcting term alpha by the square of the number of registers m, then divides by the sum of 2^(-register_value) across all registers. Because each register only holds the largest run of zeros encountered, this set of counters succinctly summarizes the distribution of hash values in the data. The from_iter constructor provides a convenient way to build a HyperLogLog by passing in an iterator of elements, while the merge method allows two sketches of the same precision to be combined, an operation that is particularly useful in large-scale or distributed applications. Finally, the main function demonstrates how to apply HyperLogLog to integer values, string data with duplicates, and how merging two separate sketches approximates the count of their combined sets.

#### Files contents:
* experiment_51_1/
  * Cargo.toml (Cargo.toml file for dependencies)
*experiment_51_1/src/
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

