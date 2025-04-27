# 8.1. Fundamentals of Genetic Variation

## 8.1.3 Practical Perspectives and HPC Concurrency

### experiment_8_1_3

Below is an illustrative Rust code snippet that computes genotype frequencies, performs a chi-square–based HW p-value calculation, and exploits concurrency to handle separate chunks of the genome. This example is adapted for industrial-scale usage by incorporating crates that enhance numerical stability, data handling, and concurrency.

In this code, rust-htslib is employed for reading and parsing VCF/BCF files. The rayon crate can be used for parallel iterators if multiple shards of the genome need processing concurrently. ndarray assists with numerical operations for genotype frequency calculations, while statrs ensures robust statistical distributions. polars offers efficient DataFrame manipulation to keep results organized in memory. With thoughtful adjustments to memory usage, chunk sizes, and concurrency levels, this approach scales up to industrial-size datasets typically encountered in large consortia or AI-driven pharmaceutical pipelines.

#### Files contents:
* experiment_8_1_3/
  * Cargo.toml (Cargo.toml file for dependencies)
* experiment_8_1_3/src/
  * main.rs (rust script)
  * synthetic.vcf (synthetic vcf file for input file)
  * synthetic.vcf.hw_results.csv (synthetic.vcf result csv file)
  * output.txt (text file output)
* experiment_8_1_3/src/work/65/54b3be71f96f81cbb7987c87cd42f1/
  * test.coverage.json (test coverage json file)
* experiment_8_1_3/src/work/f0/c029dbc59728cf12ca4ea10d38edb5/
  * test.coverage.json (test coverage json file)
* experiment_8_1_3/src/work/fb/a19c2d0203adcbff8bc1a8c54dc6c6/
  * merged_coverage.json (merged coverage json file)

#### How to run:

run main.rs in wsl:

```wsl
cargo run -- synthetic.vcf 0 1000000
```

(run main.rs and create synthetic.vcf.hw_results.csv output)

#### [dependencies]

```toml
rust-htslib = "0.49.0"
rayon = "1.5.1"
ndarray = "0.16.1"
statrs = "0.18.0"
polars = { version = "0.46", features = ["lazy"] }
```

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

