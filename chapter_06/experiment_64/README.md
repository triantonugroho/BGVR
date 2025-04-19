## 6.4. Parallel and Distributed Processing of HTS Data

### experiment_64

The following Rust program is demonstrating parallelized read counting in a BAM file across multiple genomic regions. The code is designed for both HPC cluster schedulers and cloud platforms, making it suitable for ephemeral compute instances that rapidly process a subset of data before shutting down. The snippet uses the rust-htslib crate for reading BAM files, rayon for concurrency, and includes enhancements for clearer error handling and logging. Additional strategies, such as streaming data through channels or adding advanced numeric libraries (e.g., ndarray, linfa) or deep learning frameworks (e.g., tch-rs), can be easily built on top of this foundation.

In this improved version, the process_bam_chunk function returns a Result type from the anyhow crate, consolidating error messages into a single error object that can be tracked or retried if necessary. Each region is processed in parallel using rayon’s parallel iterators, leveraging multiple CPU cores or nodes in an HPC cluster. Failures in specific regions do not crash the entire program; rather, they are logged, and the pipeline can continue processing other intervals.

For industrial-scale usage, teams often add even more sophisticated features. These can include structured logging (e.g., JSON logging) for ingestion by observability stacks such as ELK or OpenTelemetry, or the use of channels to stream partial results into a dedicated writer thread. In HPC or cloud setups, containerizing this Rust tool ensures a reproducible runtime environment with minimal overhead and simplifies deployment to popular orchestration systems like Nextflow or Kubernetes. If the analysis calls for more intensive numeric or machine-learning tasks, libraries like ndarray or tch-rs can be seamlessly integrated, all while retaining Rust’s memory safety and performance benefits.

The Nextflow script below orchestrates multiple stages: fetching BAM files (or alignment), running a Rust-based variant-calling or read-counting tool, and finally merging the resulting outputs. Each stage runs in an ephemeral container or environment, leveraging HPC or cloud resources for rapid parallelization. The Rust binary, compiled from the code you developed earlier, can be substituted in place of rust_caller_tool. Because Rust binaries are typically statically linked, they introduce minimal overhead when spun up in containerized Nextflow processes.

This pipeline demonstrates how Nextflow simplifies the orchestration of complex tasks. The alignmentOrFetch stage pulls or generates BAM files for each sample. The variantCalling stage then applies a Rust-based tool, which can execute parallel read-counting or variant-calling logic using the rayon crate, ensuring multi-core efficiency. Finally, mergeVcfs aggregates all VCF outputs into a single file, suitable for downstream steps like annotation or population-level analyses.

Many AI engineers and bioinformaticians have adopted these Rust and Nextflow solutions to power scalable genomics pipelines, enabling them to run cost-effectively on public clouds or local HPC clusters. One success story involves a large hospital system that processes thousands of clinical exomes per week, using ephemeral containers to align reads, call variants, and combine results seamlessly. Rust’s concurrency guarantees helped them avoid race conditions that might otherwise have caused data corruption in long, multi-threaded analyses. They also appreciated the minimal overhead of static binaries, which improved performance and reduced container image sizes.

#### Files contents:
* experiment_64/
  * Cargo.toml (Cargo.toml file for dependencies)
* experiment_64/src/
  * main.rs (rust script)
  * main.nf (nextflow script)
  * ref.fasta (reference fasta file as input file)
  * ref.fasta.fai (indexed ref.fasta.fai)
  * sample1.bam (sample 1 bam file as input file)
  * sample1.bam.bai (indexed sample1.bam file)
  * sample2.bam (sample 2 bam file as input file)
  * sample2.bam.bai (indexed sample2.bam file)
  * samples.txt (samples text file as input file)
  * output.txt (text file output)
* experiment_64/target/debug/
  * bam_read_counter.rar (compressed bam_read_counter execution file output from running main.rs)
* experiment_64/src/results/
  * merged.vcf (merged vcf output file)

#### How to run:

run main.rs in wsl:

```wsl
cargo run -- --bam sample1.bam --region 'chr1:1-10' 'chr1:1-32' 2>&1 | tee output.txt
```

(run main.rs with sample1.bam input, region : chr1:1-10, chr:1-32 parameters and save the output in output.txt)

run main.nf in wsl:

```wsl
nextflow run main.nf
```

run main.nf with this parameters:
params.sample_list = "samples.txt"
params.output_dir = "results"
params.region = "chr1:1-32"
params.mock = true  // Set to true to use mock commands instead of actual tools

#### [dependencies]

```toml
clap = { version = "4.0", features = ["derive"] }
anyhow = "1.0"
log = "0.4"
rayon = "1.7"
rust-htslib = "0.49.0"
```

#### Explanation of the Output
We ran a complete Nextflow pipeline (main.nf) to simulate a variant calling workflow using mock data and a custom Rust tool (bam_read_counter) for read processing. Here's a breakdown of what each part of the pipeline did and what the final output means:

##### 🧪 Input

* samples.txt:

```text
sample1
sample2
```

This file lists the sample IDs to process.

* BAM files:

 * sample1.bam: Contains alignments from sample1.sam, already tested using your Rust tool.

 * sample2.bam: Assumed to exist with similar structure.

* ref.fasta: Not used in mock mode, but assumed as reference in real variant calling.

##### 🛠️ Workflow Steps Summary

###### 1. Channel Setup

* samples_ch: Reads each line of samples.txt into the workflow as individual sample IDs.

###### 2. alignmentOrFetch Process

* Copies BAM files (sample1.bam, sample2.bam) from a fixed directory into the working directory for processing.

###### 3. variantCalling Process
* Mock mode is ON: Instead of real variant calling, it generates dummy VCF files per sample with 3 hardcoded variants.

* For both sample1 and sample2, a VCF file is created with:

```rust
chr1	14653	.	A	G	100	PASS	.	GT	0/1
chr1	14907	.	A	G	100	PASS	.	GT	0/1
chr1	15211	.	G	A	100	PASS	.	GT	1/1
```

####### 4. mergeVcfs Process

* Mock merging: Simply concatenates all generated VCF files.

* Keeps one header and appends all variant records without deduplication or sample separation.

##### 📄 Output: results/merged.vcf

```text
##fileformat=VCFv4.2
##source=MockVariantCaller
##contig=<ID=chr1,length=248956422>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1
chr1	14653	.	A	G	100	PASS	.	GT	0/1
chr1	14907	.	A	G	100	PASS	.	GT	0/1
chr1	15211	.	G	A	100	PASS	.	GT	1/1
chr1	14653	.	A	G	100	PASS	.	GT	0/1
chr1	14907	.	A	G	100	PASS	.	GT	0/1
chr1	15211	.	G	A	100	PASS	.	GT	1/1
```

* Contains 6 variants total (3 from sample1, 3 from sample2).

* Header only shows one sample column (sample1) due to mock merge simplicity.

* Duplicate lines reflect repeated variant positions across samples.

#### 🧾 Conclusion

✅ Your Rust tool (bam_read_counter) successfully processed regions and read counts as intended. It's ready to be incorporated if needed.

✅ The Nextflow pipeline ran successfully, handling multiple samples and generating merged results.

🧪 Because you used params.mock = true, no real alignment or variant calling occurred. The variants are just placeholders.

🔍 The output merged.vcf confirms both samples were processed and their dummy variants were combined.

