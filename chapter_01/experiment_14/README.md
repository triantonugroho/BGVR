# 1.4. Pangenome Graph Theorems

In the scenario below, we imagine having multiple genomic haplotypes, each stored in a separate FASTA file. We wish to merge these sequences into a single pangenome structure that captures shared subpaths and preserves haplotype-specific variants. By leveraging Rust’s concurrency features (via “rayon”), we can scale to large sets of FASTA files while ensuring memory safety and efficiency. In addition, we employ succinct indexing concepts by labeling subpaths (k-mers) and merging them into equivalence classes that unify homologous regions of different haplotypes.

The following three synthetic FASTA files exemplify distinct haplotypes. Each file should be stored under src/

>haplotype1_chr1

ACGTACGTACGTACGTACGTACGTACGTTTGGGCCCACGTACGTACGTAAAAC

>haplotype1_chr2

TTGGGCCCACGTACGTACGTACGTAAAACACGTACGTACGTACGTACGTGGGC

>haplotype2_chr1

ACGTACGTACGTACGTACGTACGTACGTTTTTGAAAACCCTGGGACGTACGTA

>haplotype2_chr2

TTTTGAAAACCCTGGGACGTACGTACGTACGTACGTACGTACGTTTTTGAAAA

>haplotype3_chr1

ACGTACGTACGTACGTACGTACGTACGTTTTGGGCCCCCCCACGTACGTACGT

>haplotype3_chr2

TTTGGGCCCCCCCACGTACGTACGTCACGTTTTGAAAACGTACGTACGTACGA

Once these three FASTA files have been saved to the src/ folder, the next step is to configure the Rust project dependencies and provide the main program logic. The first snippet below illustrates the necessary entries for your Cargo.toml. This is then followed by the main Rust code, which reads the FASTA files, constructs overlapping k-mers, and merges them into a single pangenome graph.

## Files contents:
* main.rs (rust script)
* main.nf (nextflow script)
* haplotype1.fasta (1st fasta file)
* haplotype2.fasta (2nd fasta file)
* haplotype3.fasta (3rd fasta file)

## How to run:

cargo run main.nf (run the nextflow script that will run the main.rs and save the output in output.txt)
  
## [dependencies]

bio = "2.0.3"

rayon = "1.10.0"

## Explanation of the Output

The output represents the result of constructing a pangenome graph from multiple haplotype sequences stored in FASTA files. Below is a breakdown of the key parts of the output:

### 1. Pangenome Graph Construction

Constructed a pangenome graph with 72 nodes

* The program reads three haplotype FASTA files (haplotype1.fasta, haplotype2.fasta, and haplotype3.fasta).
* It constructs a pangenome graph, where:
  * Each node represents a k-mer (a substring of length k = 21).
  * Each edge represents an overlap between consecutive k-mers in the sequences.
* A total of 72 nodes (unique k-mers) were created in the graph.
  
### 2. Sample Nodes and Their Connections

  Node: CTAGCTAGCTAGCTAAGCTAG -> ["TAGCTAGCTAGCTAAGCTAGC"]
  Node: GCTAGCTAGCTAGCTAGCTAC -> ["CTAGCTAGCTAGCTAGCTACT"]
  Node: GCTAGCTAGCTAGCTAAGCTA -> ["CTAGCTAGCTAGCTAAGCTAG"]
  Node: GCTAGCTACTAGCTAGCTAGC -> ["CTAGCTACTAGCTAGCTAGCT"]
  Node: TAGCTAGCTACTAGCTAGCTA -> ["AGCTAGCTACTAGCTAGCTAG"]

* Each line represents a node (k-mer) and its connected k-mers (edges).
  
* Example:

  Node: CTAGCTAGCTAGCTAAGCTAG -> ["TAGCTAGCTAGCTAAGCTAGC"]

  * The node CTAGCTAGCTAGCTAAGCTAG is a 21-base k-mer.
  * It connects to TAGCTAGCTAGCTAAGCTAGC, which is the next overlapping k-mer in the sequence.
* This pattern repeats for other k-mers, forming a graph structure where edges represent sequence continuity.

### 3. Conclusion

* The program successfully builds a pangenome graph by extracting k-mers from multiple haplotypes.
* The graph has 72 unique k-mer nodes.
* Each node connects to one or more k-mers, forming a directed graph.
* The output file (output.txt) contains the full graph structure, showing how different k-mers are linked.

This approach allows for comparative genomics by capturing shared and unique sequences across different haplotypes.

