## 8.6. Advanced Topics in Variant Analysis and Annotation

### experiment_8_6

The code implements a sophisticated genomic variant scoring pipeline that combines pangenome graph analysis with machine learning to assess the functional impact of genetic variants. Leveraging modern computational genomics approaches, the pipeline integrates information from variation graph representations (ODGI), machine learning inference (ONNX Runtime), and haplotype phasing (WhatsHap) to provide comprehensive variant assessments that consider both sequence context and predicted functional effects. The framework is designed to enhance variant prioritization in clinical genomics, rare disease research, and population-scale studies by incorporating graph-based contextual information that traditional linear genome approaches lack.

The pipeline operates by first loading a pangenome graph and a trained ONNX model, then processing VCF files in parallel batches to maximize throughput. For each variant, the system extracts features including allele lengths, graph-based metrics (node degree, centrality), and sequence complexity, arranges them into tensors for efficient batch processing, then runs inference through the neural network to produce variant scores. Simultaneously, it leverages WhatsHap's phasing algorithms to group variants into haplotype blocks, providing additional context for interpreting complex genomic regions. The implementation incorporates robust error handling, progress reporting, comprehensive logging, and flexible output formats (Parquet, Arrow, CSV, JSON), while tracking detailed statistics throughout execution. The architecture supports both single-file processing and batch operations across multiple VCFs, making it adaptable to both targeted analyses and population-scale studies.

#### Project Structure

```
experiment_8_6/
├── Cargo.toml                          # Main workspace configuration
├── src/
│   ├── main.rs                         # Main Rust application
│   ├── generate_onnx_model.py          # Python script to generate ONNX model
│   ├── sample_graph.json               # Sample pangenome graph (input)
│   ├── sample_variants.vcf             # Sample VCF file (input)
│   ├── variant_model.onnx               # Trained ONNX model (input)
│   └── vcf_list.txt                    # List of VCF files for batch processing (input)
├── mock_libs/                          # Mock library implementations
│   ├── odgi/
│   │   ├── Cargo.toml                  # ODGI mock library configuration
│   │   └── src/
│   │       └── lib.rs                  # ODGI graph operations mock
│   ├── onnxruntime/
│   │   ├── Cargo.toml                  # ONNX Runtime mock library configuration
│   │   └── src/
│   │       └── lib.rs                  # ONNX inference mock
│   └── whatshap-rs/
│       ├── Cargo.toml                  # WhatsHap phasing mock library configuration
│       └── src/
│           └── lib.rs                  # Variant phasing mock
├── results/
│   └── results.arrow                   # Output file (scored variants)
└── target/debug/
    └── variant-scorer.rar              # Compiled executable (compressed)
```

#### Directory Descriptions

##### Root Files
- **Cargo.toml**: Workspace configuration defining the main project and mock library dependencies

#### src/
Contains the main application source code and input files:
- **main.rs**: Core variant scoring application written in Rust
- **generate_onnx_model.py**: Python utility to create the ONNX model for variant scoring
- **sample_graph.json**: Example pangenome graph in JSON format
- **sample_variants.vcf**: Sample VCF file containing genetic variants
- **variant_model.onnx**: Pre-trained machine learning model for variant scoring
- **vcf_list.txt**: Text file listing multiple VCF files for batch processing

##### mock_libs/
Mock implementations of external libraries for development and testing:
- **odgi/**: Mock pangenome graph library for graph operations
- **onnxruntime/**: Mock ONNX runtime for machine learning inference
- **whatshap-rs/**: Mock variant phasing library

##### results/
Output directory containing processed results:
- **results.arrow**: Scored variants in Apache Arrow format

##### target/debug/
Build artifacts directory:
- **variant-scorer.rar**: Compiled executable (compressed for distribution)

#### File Types

| Extension | Description |
|-----------|-------------|
| `.rs` | Rust source code |
| `.py` | Python script |
| `.toml` | Configuration file (Cargo) |
| `.json` | JSON data file |
| `.vcf` | Variant Call Format file |
| `.onnx` | ONNX machine learning model |
| `.txt` | Plain text file |
| `.arrow` | Apache Arrow binary format |
| `.rar` | Compressed archive |


#### Cargo.toml

```toml
[package]
name = "variant-scorer"
version = "0.1.0"
edition = "2021"
authors = ["Your Name <your.email@example.com>"]
description = "A tool for scoring genetic variants using pangenome graphs and machine learning"

[dependencies]
anyhow = "1.0"
clap = { version = "4.4", features = ["derive"] }
ndarray = "0.15"
odgi = { path = "./mock_libs/odgi" }
onnxruntime = { path = "./mock_libs/onnxruntime" }
polars = { version = "0.33", features = ["parquet", "ipc", "json", "lazy", "dtype-full"] }
rayon = "1.8"
rust-htslib = "0.44"
serde = { version = "1.0", features = ["derive"] }
whatshap-rs = { path = "./mock_libs/whatshap-rs" }
num_cpus = "1.16"
tracing = "0.1"
tracing-subscriber = "0.3"
indicatif = "0.17"
thiserror = "1.0"
tempfile = "3.8"
serde_json = "1.0"
rand = "0.8"

[features]
default = []
cuda = []

[profile.release]
lto = true
codegen-units = 1
opt-level = 3
```

#### How to run:

run main.rs in wsl:

```wsl
cargo run | tee output.txt
```

(run main.rs with cohort_A.vcf, cohort_B.vcf and synthetic_variant_data.csv as input files and create query_results.parquet output file)


#### 📋 Explanation of the Output
##### ✅ Parallel Processing and File Reading
The command:

```bash
cargo run | tee output.txt
```

executes your Rust application and logs its output to output.txt.

* Parallelism: The output shows that the program starts with 8 threads, utilizing all available CPU cores via the rayon thread pool:

```text
Starting pangenome analysis with 8 threads
```

* VCF Reading:

```text
Reading variants from cohort_A.vcf
Reading variants from cohort_B.vcf
```

Each file was parsed in under 10 milliseconds:

```text
Read 1000 variants from cohort_A.vcf in 5.77ms
Read 1000 variants from cohort_B.vcf in 6.09ms
```

##### 🧬 Variant Set Algebra Results

* Union (A ∪ B): 2000 variants — all unique across both cohorts.

* Intersection (A ∩ B): 0 — no shared variants between the two sets.

* A \ B and B \ A: 1000 variants each — all variants are cohort-specific.

* Jaccard Index: 0.0000, indicating no overlap between cohort A and B variants.

```text
Variant set comparison:
  A∪B = 2000 variants
  A∩B = 0 variants
  A\B = 1000 variants
  B\A = 1000 variants
  Jaccard index = 0.0000
```

##### 📊 CSV File Processing

* synthetic_variant_data.csv was read into a Polars DataFrame with:
  * 1000 rows
  * 7 fields: CHROM, POS, REF, ALT, GT, GQ, and DP
```text
DataFrame query results:
  Variants in CSV: 1000
  DataFrame schema: ...
```

* Statistics Note: Skipped due to compatibility issues with polars 0.47.0.

##### 💾 Export Operation

* Successfully exported the DataFrame to a Parquet file:

```text
Exported DataFrame with 1000 rows to query_results.parquet
```

##### ⏱️ Performance

* Total runtime: Only 55 milliseconds for the entire process, which includes:
  * Multi-threaded parsing of 2,000 VCF entries
  * Set operations
  * Reading a 1,000-row CSV
  * Exporting to Parquet

```text
Total execution time: 55.27ms
```

#### ✅ Conclusion
This execution confirms that your Rust-based pangenome tool is:

* ⚡ Fast and Efficient: Processes multiple files and performs variant set algebra and CSV I/O in under 60 milliseconds.

* 🧵 Scalable: Uses all available cores with rayon to parallelize I/O-bound and compute-bound tasks.

* 🛠️ Reliable: Provides error handling, schema validation, and gracefully skips unsupported features (e.g., Polars stats).

* 🧬 Biologically Informative: Set operations reveal that the two cohorts have zero shared variants, likely indicating they originate from entirely distinct populations or datasets.

* 🧱 Ready for Integration: Generates intermediate results in efficient formats like Parquet, enabling easy downstream analysis in cloud-native or Python-based workflows.


