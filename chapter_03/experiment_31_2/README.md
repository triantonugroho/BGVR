## 3.1. Introduction to Data Structures and Algorithms

### experiment_31_2

This code demonstrates a parallel, thread-safe MinHash implementation in Rust that computes approximate set similarities on synthetic genomic data. The program generates random DNA sequences of a specified length and uses them to create two datasets. Each dataset is processed by a MinHasher that applies multiple hash seeds in parallel, determining a “minimum hash value” across all items for each seed. This approach is particularly relevant for high-throughput genomic or large-scale text data, where efficiently estimating set overlap becomes necessary.

Internally, the MinHasher maintains a specified number of hash seeds. For each item in a dataset, it computes a distinct hash based on the item and seed, retains only the smallest hash value for that seed, and repeats the process for all seeds. The resulting collection of smallest hash values, one per seed, becomes the MinHash signature. Two datasets can then be compared by examining their signatures: the fraction of identical positions in the signatures approximates the Jaccard similarity of the original sets. The code uses Rayon’s parallel iterators to distribute these computations across CPU cores, thereby increasing throughput in scenarios where large numbers of seeds or large datasets are involved.

#### Files contents:
* experiment_31_2/Cargo.toml ((Cargo.toml file for dependencies)
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
rand = "0.9.0"
rayon = "1.10.0"
```
Explanation of the Output:
The program computes MinHash signatures for two sets of synthetic DNA sequences and calculates the approximate Jaccard similarity between them.

MinHash Signature for set1:

The program generates the MinHash signature for the first set of DNA sequences. The signature consists of a list of hash values, one for each of the 100 hash functions (or "seeds").
The first few values are:
csharp
Copy
Edit
[2942580149134852, 6899882479743125, 21997229414397169, 3675138294618832, 674304587691108]
These values represent the minimum hash values computed for each hash function applied to the sequences in the first set.
MinHash Signature for set2:

Similarly, the program computes the MinHash signature for the second set of DNA sequences.
The first few values for set2 are:
csharp
Copy
Edit
[5375613036628250, 1580167804507962, 10283091463069476, 1055940623035697, 8502197830624736]
Approximate Jaccard Similarity:

The Jaccard similarity measures the similarity between two sets by comparing the number of common elements to the total number of unique elements.
In this case, the MinHash signatures are used to estimate the Jaccard similarity in a computationally efficient way.
The similarity between sig1 (the signature of set1) and sig2 (the signature of set2) is calculated.
The result is:
java
Copy
Edit
Approximate Jaccard similarity = 0.000
This means that, according to the MinHash signatures, the two sets are very dissimilar. The similarity score of 0.000 indicates that there are no matching hash values across all the hash functions between the two sets.
Conclusion:
MinHash Signatures are a useful technique for estimating the similarity between two large sets in a way that is computationally efficient. In this case, the program demonstrates this by comparing synthetic genomic data (DNA sequences).
The Jaccard similarity between the two sets in this example is very low (0.000), indicating that the two sets are very different according to their MinHash signatures. This is expected, as the two sets are generated independently and consist of random sequences.
The program demonstrates how parallelism (rayon::prelude::*) can speed up the computation of hash values and the MinHash signature, making it suitable for large datasets.
In real-world applications, MinHash is often used for tasks like duplicate detection, clustering, and recommendation systems where sets may be too large to compare directly.


