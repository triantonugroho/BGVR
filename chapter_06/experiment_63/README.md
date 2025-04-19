## 6.3. Variant Call Format (VCF/BCF) Handling

### experiment_63

Below is a production-level example that demonstrates how to filter a BCF file in parallel using rust-htslib, rayon, and additional Rust crates for robust error handling and logging. It processes records in chunks to avoid excessive memory use, which is critical for large population-level datasets. This setup is well-suited for HPC environments and cloud-based workflows alike, where reliability, speed, and memory efficiency are of paramount importance.

In this updated code, the CLI parsing is handled by the clap crate, making it easy to specify input files, output files, and filter parameters on the command line. The anyhow and env_logger crates provide more robust error handling and logging, ensuring that issues such as missing files or corrupt data are reported clearly. By processing the file in configurable chunks, the memory footprint remains more predictable, which is crucial when operating on large datasets often found in population genomics projects.

Each batch of records is filtered in parallel using rayon’s parallel iterator, which leverages multiple CPU cores. This design scales efficiently in HPC contexts and cloud environments with autoscaling capabilities, allowing you to handle massive variant datasets without rewriting core logic. The thread-safety and ownership guarantees provided by Rust eliminate data races and memory corruption issues that can arise under high concurrency, preserving both correctness and performance.

The following is a Nextflow workflow snippet aligned with the previously discussed Rust code. This example uses bcftools to slice out a specific genomic interval (the “chunk”) from a larger BCF file and then pipes the resulting records into our Rust-based bcf_filter_tool. Each chunk is processed by an ephemeral Nextflow task, which allows for highly parallelized operations on an HPC cluster or cloud environment. The ephemeral containers used by Nextflow ensure that each task spins up with just enough resources to process a chunk, then gracefully shuts down.

In this workflow, bcftools view isolates the records corresponding to the specified chunk, which might be a chromosome interval like chr1:1-1000000. That subset of variant data is then piped into the bcf_filter_tool, which corresponds to the Rust code discussed earlier. By specifying --input /dev/stdin and --output /dev/stdout, the tool reads directly from the pipeline’s stream and writes filtered records back to standard output, which Nextflow redirects into a file named filtered_${chunk}.bcf.

For large genomic datasets spanning many intervals, Nextflow automatically distributes these tasks across multiple compute nodes, each running an ephemeral container or environment. Once a chunk is processed, the task finishes, freeing resources for other intervals or other pipelines. This approach is highly scalable, minimizing total runtime for population-level analyses. Subsequent steps in Nextflow might merge the per-chunk BCF files back together using bcftools concat or a custom merging process, maintaining data integrity and allowing further population-wide computations.

By combining the memory-safe, high-performance Rust filtering tool with Nextflow’s orchestration, teams can tackle large workflows while ensuring that each step is both resilient (capable of retries on transient errors) and efficient in resource consumption. More advanced use cases can incorporate deep learning packages via the tch-rs crate or numerical libraries like ndarray and linfa for sophisticated analyses, all within the same Rust-based codebase.

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

Execution | Input Size | Output Files | Parallelism | Use Case
main.rs (cargo run) | Full BCF | filtered_output.bcf | Parallel per batch | Full file filtering, standalone run
main.nf (Nextflow) | Full BCF split by region | One output per region: filtered_<chunk>.bcf | Parallel per chunk | Scalable filtering over genomic regions

🧾 Key Takeaways:

* Tool is Modular and Efficient: It's reusable in both standalone and pipeline contexts, and benefits from multithreading.

* Nextflow Pipeline Adds Scalability: Makes it easier to split tasks and parallelize across a cluster or local cores.

* You’re Setup for Bigger Datasets: With chunks_file, you can scale to thousands of regions (like full chromosomes or exomes).

* Flexibility: You can update main.nf later to generate regions dynamically, or use actual genomic intervals from .bed or .fai files.
