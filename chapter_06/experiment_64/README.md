## 6.4. Parallel and Distributed Processing of HTS Data

### experiment_64

The following Rust program is demonstrating parallelized read counting in a BAM file across multiple genomic regions. The code is designed for both HPC cluster schedulers and cloud platforms, making it suitable for ephemeral compute instances that rapidly process a subset of data before shutting down. The snippet uses the rust-htslib crate for reading BAM files, rayon for concurrency, and includes enhancements for clearer error handling and logging. Additional strategies, such as streaming data through channels or adding advanced numeric libraries (e.g., ndarray, linfa) or deep learning frameworks (e.g., tch-rs), can be easily built on top of this foundation.

In this improved version, the process_bam_chunk function returns a Result type from the anyhow crate, consolidating error messages into a single error object that can be tracked or retried if necessary. Each region is processed in parallel using rayon’s parallel iterators, leveraging multiple CPU cores or nodes in an HPC cluster. Failures in specific regions do not crash the entire program; rather, they are logged, and the pipeline can continue processing other intervals.

For industrial-scale usage, teams often add even more sophisticated features. These can include structured logging (e.g., JSON logging) for ingestion by observability stacks such as ELK or OpenTelemetry, or the use of channels to stream partial results into a dedicated writer thread. In HPC or cloud setups, containerizing this Rust tool ensures a reproducible runtime environment with minimal overhead and simplifies deployment to popular orchestration systems like Nextflow or Kubernetes. If the analysis calls for more intensive numeric or machine-learning tasks, libraries like ndarray or tch-rs can be seamlessly integrated, all while retaining Rust’s memory safety and performance benefits.

The Nextflow script below orchestrates multiple stages: fetching BAM files (or alignment), running a Rust-based variant-calling or read-counting tool, and finally merging the resulting outputs. Each stage runs in an ephemeral container or environment, leveraging HPC or cloud resources for rapid parallelization. The Rust binary, compiled from the code you developed earlier, can be substituted in place of rust_caller_tool. Because Rust binaries are typically statically linked, they introduce minimal overhead when spun up in containerized Nextflow processes.

This pipeline demonstrates how Nextflow simplifies the orchestration of complex tasks. The alignmentOrFetch stage pulls or generates BAM files for each sample. The variantCalling stage then applies a Rust-based tool, which can execute parallel read-counting or variant-calling logic using the rayon crate, ensuring multi-core efficiency. Finally, mergeVcfs aggregates all VCF outputs into a single file, suitable for downstream steps like annotation or population-level analyses.

Many AI engineers and bioinformaticians have adopted these Rust and Nextflow solutions to power scalable genomics pipelines, enabling them to run cost-effectively on public clouds or local HPC clusters. One success story involves a large hospital system that processes thousands of clinical exomes per week, using ephemeral containers to align reads, call variants, and combine results seamlessly. Rust’s concurrency guarantees helped them avoid race conditions that might otherwise have caused data corruption in long, multi-threaded analyses. They also appreciated the minimal overhead of static binaries, which improved performance and reduced container image sizes.

#### Files contents:
* experiment_63/
  * Cargo.toml (Cargo.toml file for dependencies)
* experiment_63/src/
  * main.rs (rust script)
  * main.nf (nextflow script)
  * chunks_list.txt (chunk list input text file)
  * wgs_cohort.bcf (bcf input file)
  * wgs_cohort.bcf.csi (indexed bcf input file)
  * filtered_output.bcf (bcf output file from running main.rs)
* experiment_63/target/debug/
  * bcf_filter_tool.rar (compressed bcf_filter_tool execution file output from running main.rs)
* experiment_63/src/work/69/fd2fa2cc9f7401386d609c509b0544/
  * filtered_chr1_1000-1200.bcf (bcf output file from running main.nf)

#### How to run:

run main.rs in wsl:

```wsl
cargo run -- --input wgs_cohort.bcf --output filtered_output.bcf --min-qual 30.0 --min-depth 10 --chunk-size 50000
```

(run main.rs with input wgs_cohort.bcf and output file name filtered_ouput.bcf with params.min_qual = 30, params.min_depth = 10, and params.chunk_size  = 50000)

run main.nf in wsl:

```wsl
nextflow run main.nf
```

(run main.nf with input wgs_cohort.bcf and output file name filtered_ouput.bcf with params.min_qual = 30, params.min_depth = 10, and params.chunk_size  = 50000)

#### [dependencies]

```toml
anyhow = "1.0"
clap = { version = "4.4", features = ["derive"] }
env_logger = "0.11.7"
log = "0.4"
rayon = "1.8"
rust-htslib = "0.49.0"
```

#### Explanation of the Output

▶️ Running main.rs via Cargo

Command:

```wsl
cargo run -- --input wgs_cohort.bcf --output filtered_output.bcf --min-qual 30.0 --min-depth 10 --chunk-size 50000
```

🔍 Output:

* A single output BCF file: filtered_output.bcf

* Contains all records from wgs_cohort.bcf that:

  * Have QUAL >= 30

  * Have average DP (depth) across samples >= 10

* Processing is done in parallel batches of 50,000 records using Rayon.

* Log output (if using RUST_LOG=info) will show:

  * Total records processed

  * Total records retained after filtering

✅ Interpretation:

This execution processes the entire BCF file at once (no genomic region splitting). It’s great for:

* Benchmarking speed of the filter

* Validating correctness of filtering logic

* Producing one consolidated result

▶️ Running main.nf via Nextflow

Command:

```wsl
nextflow run main.nf
```

🔍 Output:

* One BCF files, named:

  * filtered_chr1_1000-1200.bcf

✅ Interpretation:

* Instead of processing the full file at once, this splits the input by genomic regions (chunks).

* Each chunk is processed in parallel via the filterVariants process.

* All filtering is still done using your Rust binary, with the same filtering thresholds.

* If a chunk has no variants or no variants passing the filters, its corresponding output .bcf may be empty or missing.

🧠 Conclusion and Insights

| Execution           | Input Size                | Output Files                                | Parallelism         | Use Case                                      |
|---------------------|----------------------------|----------------------------------------------|----------------------|-----------------------------------------------|
| `main.rs` (cargo run) | Full BCF                   | `filtered_output.bcf`                        | Parallel per batch   | Full file filtering, standalone run           |
| `main.nf` (Nextflow)  | Full BCF split by region   | One output per region: `filtered_<chunk>.bcf` | Parallel per chunk   | Scalable filtering over genomic regions       |

🧾 Key Takeaways:

* Tool is Modular and Efficient: It's reusable in both standalone and pipeline contexts, and benefits from multithreading.

* Nextflow Pipeline Adds Scalability: Makes it easier to split tasks and parallelize across a cluster or local cores.

* You’re Setup for Bigger Datasets: With chunks_file, you can scale to thousands of regions (like full chromosomes or exomes).

* Flexibility: You can update main.nf later to generate regions dynamically, or use actual genomic intervals from .bed or .fai files.

