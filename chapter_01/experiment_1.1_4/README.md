## 2.1. Introduction to Rust Programming Language

### experiment_21_4

In many bioinformatics workflows, one might receive a set of files representing different stages of an analysis pipeline: a FASTQ file with raw reads, a BAM file showing how those reads align to a reference genome, and a VCF file listing detected variants. To illustrate a unified approach, the following Rust code processes each file format in parallel to find a specific DNA motif (e.g., “GATTACA”). While real-world projects often require more complex logic (like checking CIGAR alignments in BAM or genotype fields in VCF), this demo showcases how Rust can handle these formats consistently and efficiently, leveraging multi-threaded execution through Rayon.

The code first defines a helper function, count_occurrences, to locate overlapping instances of a target motif in a string. It then provides three specialized functions—process_fastq, process_bam, and process_vcf—each handling a different file format using community crates like bio (for FASTQ) and rust-htslib (for SAM/BAM and VCF). In each function, Rust reads the input records, extracts the relevant sequences or alleles into a Vec<String>, and then calls par_iter() from Rayon to distribute the motif-counting workload across available CPU cores. Finally, main orchestrates these steps by specifying a motif (e.g., “GATTACA”), applying each function to its corresponding file, and printing out how many total motif matches are found. This design cleanly separates reading logic from parallel processing, remains memory-safe, and can be scaled or adapted to handle larger inputs or more complex analyses without sacrificing performance or maintainability.

#### Project Structure:

```plaintext
experiment_21_4/
├── Cargo.toml                     # Rust project configuration and dependencies
└── src/
    ├── main.rs                    # Main Rust script containing program logic
    ├── main.nf                    # Nextflow workflow script
    ├── example.fastq.rara         # Compressed FASTQ example file (note: unusual .rara extension)
    ├── example.bam                # BAM alignment file
    ├── example.vcf                # VCF variant call file
    └── output.txt                 # Output file
```

#### How to run:

```nextflow
nextflow run main.nf | tee output.txt
```

(run main.nf in WSL terminal and save the output in output.txt)

#### [dependencies]

```toml
bio = "2.0.3"
rayon = "1.10.0"
rust-htslib = "0.49.0"
```

#### Explanation of the Output

##### 1. FASTQ: Found motif 'GATTACA' 0 times.

* The program read sequences from example.fastq and searched for occurrences of "GATTACA".
* Since the count is 0, it means the motif was not found in any of the sequences in the FASTQ file.
* This could be due to the motif being absent in the dataset, the sequences being too short, or having mutations that prevent exact matches.

##### 2. BAM: Found motif 'GATTACA' 230 times.

* The program processed example.bam (a binary alignment file containing mapped sequencing reads).
* "GATTACA" was detected 230 times across all aligned reads.
* This suggests that the motif appears frequently in the aligned sequencing data.
* The presence of this motif could indicate a biologically significant region, a conserved sequence, or a sequencing artifact.

##### 3. VCF: Found motif 'GATTACA' 0 times.

* The program analyzed example.vcf, which contains variant data (reference and alternative alleles).
* Since "GATTACA" was found 0 times, it means the motif is not present in any reference or alternative alleles in the VCF file.
* This could be because the dataset does not contain variations that introduce this motif, or the motif does not overlap with the recorded variant positions.

#### Conclusion
The program successfully analyzed sequencing and variant data from different file formats (FASTQ, BAM, and VCF) to detect occurrences of the "GATTACA" motif. While the motif was absent in the raw sequencing reads (FASTQ) and variant calls (VCF), it was detected 230 times in the aligned reads (BAM). This indicates that the motif appears in mapped reads, potentially suggesting its biological relevance in the aligned genome region. Such motif searches can be applied in genome annotation, regulatory element detection, and sequence motif enrichment analysis.


















