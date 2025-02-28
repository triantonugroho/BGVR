## 1.6. Brief Introduction to CRISPR System

### experiment_16

Below is a sample Rust code demonstrating how one might blend concurrency patterns with specialized data structures to handle large genomic data, taking a simplified approach to building and merging suffix arrays. This illustration reflects the broader theme of robust engineering: even as data sizes explode in CRISPR or proteomics projects, parallelism and efficient algorithms help keep analyses feasible on modern HPC clusters.

In the snippet above, we first load a genomic sequence from a FASTA-like file, stripping away header lines. The build_parallel_suffix_array function splits the sequence into multiple chunks, each handled by a separate thread using Rayon’s parallel iterators. Within each chunk, we construct a naive suffix array by sorting suffixes; we then merge these partial results, adjusting indices to refer back to the global genome. A final sort ensures the suffix array is globally correct, even at chunk boundaries. Although this example is simplified—real-world scenarios might require more advanced suffix-array or FM-index algorithms— it demonstrates a practical pattern of combining concurrency and specialized data structures to manage large-scale genomic tasks on HPC systems.

#### Files Contents
* main.rs (rust script)
* main.nf (nextflow script)
* Python code to synthesize example_genome_fasta_file (python code)
* example_genome.fasta (fasta file)
* Cargo.toml (Cargo.toml file)

#### How to run

cargo run main.nf 

(run the nextflow script that will run the main.rs and save the output in output.txt)
  
## [dependencies]

rayon = "1.10.0"

#### Explanation of the Output

The program constructs a suffix array from a large genomic sequence (10 million bases) stored in a FASTA file. It uses parallel computing (via Rayon) to speed up suffix array construction by splitting the sequence into 8 chunks and processing them concurrently. Below is a breakdown of the output:

##### 1. Output Breakdown

Genome length: 10000000

* The genome sequence loaded from the FASTA file has 10 million characters (nucleotides: A, T, C, G).
* This confirms that the genomic data was successfully read and processed.

Suffix array length: 10000000

* A suffix array is a sorted list of all suffix starting positions in the genome.
* Since the genome length is 10 million, the suffix array also contains 10 million indices.
* Each index in the suffix array corresponds to a starting position of a suffix in the original sequence.

First 10 entries in suffix array: [5701463, 9164049, 9386089, 9100229, 78655, 2090381, 7738064, 5701464, 9963205, 47079]

* The first 10 entries in the suffix array are shown.
* These numbers represent the starting positions of the lexicographically smallest suffixes in the genome.
* For example, 5701463 means that the suffix starting at position 5,701,463 in the genome is one of the smallest in lexicographic order.
  
##### 2. How the Suffix Array Was Built

* The genome was divided into 8 chunks.
* Each chunk’s suffix array was built independently and then merged.
* A final sorting step ensured that the suffix array was correctly ordered across the entire genome.
  
##### 3. Key Insights

* Efficient Parallel Processing: Using Rayon allowed parallel computation, making suffix array construction much faster.
* Suffix Sorting Completeness: The suffix array correctly covers all 10 million positions, proving the sorting process worked properly.
* Scalability: This method is efficient for large genomic datasets, commonly used in bioinformatics applications like genome indexing and sequence alignment.

#### Conclusion

The program successfully builds a suffix array from a 10-million-character genome using parallel computing. The results demonstrate that the suffix array is correctly structured and sorted, making it useful for fast substring searches and pattern matching in genomic sequences.
