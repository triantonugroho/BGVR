## 2.1. Introduction to Rust Programming Language

### experiment_21_4

Below is a simple example that demonstrates Rust’s zero-cost abstractions in a genomic setting. This example reads a FASTA file and uses Rust’s iterator adapters (such as map, filter, and for_each) to process genomic sequences without incurring additional runtime overhead. Despite the high-level functional style, the compiler optimizes these operations down to efficient machine code, making them comparable to hand-written loops in languages like C or C++.

This code reads a FASTA file using the bio::io::fasta crate, which emits a stream of records (each containing a sequence). For each record, the sequence bytes are converted into a String via String::from_utf8_lossy, then any sequence under 50 nucleotides is discarded through a filter call. The remaining sequences move to a second map step that calculates their GC content by counting the characters ‘G’ or ‘C.’ Finally, the .sum() operation aggregates these individual GC counts into one total. Rust compiles this chain of iterator adapters into efficient loops without constructing intermediary collections, a technique referred to as “zero-cost abstraction.” Consequently, although the code is written in a high-level, functional style, it executes at speeds comparable to hand-tuned loops in lower-level languages.

#### Files contents:
* main.rs (rust script)
* main.nf (nextflow script)
* example.fastq.rara (compressed example.fastq)
* example.bam (bam file)
* example.vcf (vcf file)
* Cargo.toml (Cargo.toml file)
* output.txt (output file)

#### How to run:

nextflow run main.nf | tee output.txt

(run main.nf in WSL terminal and save the output in output.txt)
  
#### [dependencies]

```sh
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


















