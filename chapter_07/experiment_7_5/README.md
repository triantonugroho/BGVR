## 7.5. Integrating Rust Noodles into Nextflow Pipelines

### experiment_7_5

This Rust code example demonstrates how to open and index a BAM file if necessary, calculate base‐by‐base coverage over a specific genomic region using the noodles-bam and noodles-core crates, and finally serialize the results to JSON via serde. It showcases how to integrate multiple crates to perform essential bioinformatics tasks in a performant and maintainable manner, laying the groundwork for more advanced genomic operations.

In a single pass, main checks for an existing BAM index and creates one if it does not exist, then calls compute_coverage to iterate over all reads that overlap the region. That function calculates per‐base coverage by determining the overlap between each read and the query range, accumulating results into a vector of coverage counts that are serialized via serde_json for easy output or further downstream processing.

This code snippet ensures that the given BAM file has an index, then reads the specified genomic region to count how many reads fall into that range. The coverage count is then written to a JSON file, which is suitable for ingestion by downstream analytics or reporting steps. To further enhance industrial-scale robustness, logs can be collected via crates like tracing for structured logging, and the code could be encapsulated in Docker or Singularity containers with pinned versions to guarantee reproducibility.

Below is a minimal Nextflow script that orchestrates multiple coverage tasks in ephemeral containers. Each container runs the above Rust tool on different inputs, illustrating how Nextflow’s DAG scheduling synergizes with Rust’s concurrency. As with the Rust code, everything is enclosed in a panel for direct usage, and each ephemeral container can run the coverage computation in parallel, merging the JSON outputs in a final step.

Below is a minimal Nextflow script demonstrating how to run the previously described Rust coverage tool inside ephemeral containers on multiple BAM inputs, with each coverage task executed in parallel and a final step merging the JSON outputs. This setup showcases how Nextflow’s DAG scheduling seamlessly integrates with Rust’s concurrency, allowing each container‐based job to process a specific input and produce JSON results for subsequent merging.

The RUN_COVERAGE process receives tuples that include a sample identifier along with its BAM and BAI files, then spins up ephemeral containers (one per sample) running the Rust tool to compute coverage. Each container outputs a JSON file, which is collectively consumed by the MERGE_COVERAGE process. In the final step, a simple JSON merge combines these separate coverage results into a single JSON, illustrating how Nextflow’s dataflow model and Rust’s concurrency align to perform scalable, reproducible genomics analysis in containerized environments.

In this Nextflow script, coverageCalc is invoked once for each combination of BAM file and genomic region. Nextflow takes care of scheduling these tasks on available compute resources, whether those are local cores, HPC cluster nodes, or cloud instances. The ephemeral containers each run the Rust coverage tool, generating partial JSON files that are then merged in the final stage of the pipeline. This pattern, known as scatter-gather, exemplifies a standard approach for handling large genomic datasets by splitting them into manageable subsets and recombining the results.

#### Project Structure:

```plaintext
experiment_7_5/
├── Cargo.toml                  # Rust dependencies
├── src/
│   ├── main.rs                 # Rust implementation
│   ├── main.nf                 # Nextflow workflow
│   ├── test.bam                # Test BAM file
│   ├── test.bam.bai            # Indexed test.bam file
│   ├── test.sam                # Test SAM file to create test.bam
│   ├── output.txt              # Text file output
│   └── work/                   # Nextflow work directory
│       ├── 65/54b3be71f96f81cbb7987c87cd42f1/
│       │   └── test.coverage.json        # Test coverage JSON file
│       ├── f0/c029dbc59728cf12ca4ea10d38edb5/
│       │   └── test.coverage.json        # Test coverage JSON file
│       └── fb/a19c2d0203adcbff8bc1a8c54dc6c6/
│           └── merged_coverage.json      # Merged coverage JSON file
└── target/
    └── debug/
        └── rust_coverage_tool.rar  # Compressed Rust coverage tool executable
```

#### Cargo.toml

```toml
[package]
name = "rust_coverage_tool"
version = "0.1.0"
edition = "2024"

[dependencies]
noodles-bam = "0.5"
noodles-core = "0.5"
serde = { version = "1.0", features = ["derive"] }
serde_json = "1.0"
```

#### How to run:

run main.rs in wsl:

```wsl
cargo run | tee output.txt
```

(run main.rs and save the text output in output.txt)

run main.nf in wsl:

```wsl
nextflow run main.nf
```

run main.nf and create test_coverage.json and merged_coverage.json


#### Explanation of the Output
##### 🛠 main.rs (Rust Program)

What happens in main.rs:

1. Input:
   * test.bam (BAM file) and test.bam.bai (its BAI index).
2. Settings:
   * Region: chr1:10000-10100.
3. Processing:
   * Opens and reads the BAM file manually record by record (not using random access).
   * Filters only records that:
     * Map to chr1.
     * Overlap with 10000..=10100.
   * Calculates coverage:
     * For each position between 10000 and 10100, counts how many reads cover that position.
     * Coverage is an array of 101 elements (because 101 positions from 10000 to 10100 inclusive).
4. Output:
   * A JSON object like this:

```text
{
  "reference_name": "chr1",
  "start": 10000,
  "end": 10100,
  "coverage": [0, 0, 0, 0, ..., 1, 1, 0, 0, ..., 0]
}
```

* This output is exactly what you saved in test_coverage.json.

###### ✅ Main point:
The Rust program produces the coverage profile over chr1:10000-10100, and prints it to JSON.

##### 🛠 main.nf (Nextflow Workflow)

What happens in main.nf:

1. Input:
   * Takes BAM + BAI file pairs — in your example, only test.bam and test.bam.bai.

2. Step 1: RUN_COVERAGE:
   * For each BAM:
     * Runs your Rust tool.
     * Saves output JSON to ${sample_id}.coverage.json (here it becomes test_coverage.json).
3. Step 2: MERGE_COVERAGE:
   * Merges all individual JSON outputs into one single file called merged_coverage.json.
   * It does this by:
     * Opening a bracket [
     * Appending each JSON
     * Adding commas between them
     * Closing with ]
  * Since you only have one sample (test.bam), the result is a JSON array with one object inside:

```text
[
  {
    "reference_name": "chr1",
    "start": 10000,
    "end": 10100,
    "coverage": [...]
  }
]
```

###### ✅ Main point:
The Nextflow workflow automates running the Rust tool on multiple samples and merging results into a clean JSON array.

#### 📋 Summary Table

| Program         | Input Files         | Output Files                  | Content Description                        |
|-----------------|----------------------|--------------------------------|--------------------------------------------|
| `main.rs` (Rust) | `test.bam`, `test.bam.bai` | Console output, `test_coverage.json` | Coverage over `chr1:10000–10100`, 1 BAM file |
| `main.nf` (Nextflow) | `test.bam`, `test.bam.bai` | `test_coverage.json`, `merged_coverage.json` | Same coverage result, but merged into a JSON array |

#### 📢 Conclusion

* Both main.rs and main.nf produce exactly the same per-sample coverage JSON (test_coverage.json).
* main.nf just automates running main.rs for multiple BAM files, and merges results into a single merged_coverage.json.
* Since you only have one BAM (test.bam), merged_coverage.json just contains a single object inside an array.
* ✅ The content of the outputs matches — your Rust program and Nextflow workflow are working correctly and consistently.
