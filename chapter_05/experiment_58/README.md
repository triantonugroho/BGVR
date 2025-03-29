## 5.8. Integrating Short and Long Reads

### experiment_58

This Rust program demonstrates a “partial polishing” approach for hybrid assemblies, processing short-read alignments in chunks to build correction patches. It leverages rust-htslib to read BAM files, rayon for parallel iteration, and anyhow for streamlined error handling. By grouping alignment records by their contig (TID), it accumulates potential mismatches relative to a simplistic “reference base” assumption, packaging each mismatch into a CorrectionPatch.

In an HPC environment, ephemeral tasks might each receive a subset of the BAM or a partial chunk of reads, generate correction patches for their assigned data, and then store these patches in JSON files. A merge step at the end of the pipeline consolidates all partial patches into a final global correction set. The code can be further extended to incorporate advanced reference lookups, coverage-based weighting, or more precise heuristics for base calls. Rust’s concurrency guarantees and chunk-based strategy keep the pipeline memory-safe and scalable.

The program begins by parsing command-line arguments with the clap crate, allowing the user to specify the input BAM path, the chunk size for reading, the output directory for partial patch files, and the final output file for merged patches. It then creates the partial output directory if needed.

Records from the BAM file are read in segments using the helper function read_chunk, which retrieves up to a user-defined number of bam::Records per iteration. After reading a chunk, the code groups its records by TID (the numeric contig identifier) and uses Rayon’s data parallelism to examine all records associated with each TID. A naive mismatch detection step collects positions where the read base disagrees with a hypothetical reference base (assumed to be ‘A’ for demonstration) and increments a coverage count for each discrepancy. Each mismatch is stored as a CorrectionPatch.

The partial patches from each chunk are immediately serialized to JSON and written to disk. When all chunks are processed, the program scans the partial output directory, merges the patch data by concatenation, and writes the final combined results as JSON. In production HPC scenarios, additional steps might unify patches across overlapping chunk boundaries or incorporate a more nuanced reference base retrieval. Rust’s strong concurrency model, with no global data races, supports distributed or ephemeral container tasks without risking memory corruption, and the approach naturally adapts to large-scale data by distributing chunk processing across multiple nodes.


#### Files contents:
* experiment_58/
  * Cargo.toml (Cargo.toml file for dependencies)
* experiment_58/src/
  * main.rs (rust script)
  * input.bam (bam input file)
  * annotations.gtf (gtf input file)
  * merged_corrections.json (merged counts json output file)
  * output.txt (text file output)
* experiment_58/src/partial_patches
  * partial_patches_chunk_0.json (partial patches in chunk 0 json output file)
  * partial_patches_chunk_1.json (partial patches in chunk 0 json output file)
  * partial_patches_chunk_2.json (partial patches in chunk 0 json output file)
  * partial_patches_chunk_3.json (partial patches in chunk 0 json output file)

#### How to run:

run in wsl:

```wsl
cargo run -- --bam-input input.bam --gtf annotations.gtf --output merged_correction.json.json | tee output.txt
```

(run main.rs with input filename input.bam and output file name merged_correction.json output text in output.txt) 
  
#### [dependencies]

```toml
anyhow = "1.0"
bam = "0.1.2"              
rust-htslib = "0.49.0"      
rayon = "1.8.0"             
serde = { version = "1.0", features = ["derive"] } 
serde_json = "1.0"           
log = "0.4.20"               
env_logger = "0.11.7"       
clap = { version = "4.4", features = ["derive"] }  
bio = "2.2.0"   
```

#### Explanation of the Output

The output consists of several key components that provide insights into the execution of the polish_patches program:

##### 1. Partial Patch Files (partial_patches_chunk_*.json)
The program processes a BAM file (input.bam) in chunks of 10,000 records at a time. Each chunk is analyzed, and correction patches (mismatches found in sequencing reads) are stored as JSON files inside the partial_patches directory.
From output.txt, we see that:

* Chunk 0 contains 10,000 records and is saved as partial_patches_chunk_0.json.

* Chunk 1 contains 10,000 records and is saved as partial_patches_chunk_1.json.

* Chunk 2 contains 10,000 records and is saved as partial_patches_chunk_2.json.

* Chunk 3 contains 4,298 remaining records and is saved as partial_patches_chunk_3.json.

This indicates that the BAM file had 34,298 total records.

##### 2. Merged Correction File (merged_corrections.json)

After processing all chunks, the program merges all correction patches into a single JSON file (merged_corrections.json). The final output contains 8 correction patches, each represented as follows:

```json
{
  "contig": "chr1",
  "position": 12345,
  "ref_base": 65,
  "suggested_base": 67,
  "coverage": 50
}
```

Each correction patch contains:

* contig: Chromosome or sequence identifier (e.g., chr1, chr2).

* position: The genomic position where a mismatch was found.

* ref_base: The expected nucleotide at this position (ASCII representation of bases).

* suggested_base: The corrected nucleotide suggested by the program.

* coverage: The number of sequencing reads supporting the correction.

##### 3. Coverage of Mismatches

The coverage value represents how many reads support the mismatch at a given position. Higher coverage suggests a stronger confidence in the suggested correction.

#### Conclusion

1. Efficient Chunk-Based Processing

  * The program successfully processed the BAM file in chunks of 10,000 records each, ensuring memory efficiency.

  * The final chunk had 4,298 records, meaning the total BAM file contained 34,298 reads.

2. Detection of Nucleotide Mismatches

  * The algorithm detected 8 significant correction patches where the reference base differs from observed sequencing data.

  * These corrections are useful for refining genome assemblies or detecting sequencing errors.

3. Parallel Processing with Rayon

  * The use of rayon allows for parallel processing of reads, improving efficiency when handling large datasets.

4. Integration with Annotations (GTF Support)

* While the --gtf argument is present, it is not explicitly used in this version of the code.

* Future improvements could integrate GTF annotations to provide biological context (e.g., if mismatches occur in exons or regulatory regions).
