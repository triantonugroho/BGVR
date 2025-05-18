## 3.6. Putting It All Together—Rust and Nextflow Integration

### experiment_3.6_3

#### 1. Nextflow

Rust can be used alongside ABySS to enhance and customize the assembly workflow in multiple ways. First, researchers can implement specialized preprocessing routines in Rust, such as advanced quality filtering, read classification, or read sub-sampling, before passing the data to ABySS. Second, Rust’s highly optimized data structures make it suitable for post-assembly analysis, enabling efficient detection of coverage patterns or other genomic features in the resulting contigs. Finally, lightweight Rust-based command-line utilities can be developed to automate tasks like merging FASTQ files, checking file integrity, or converting between various bioinformatics file formats, thereby streamlining the entire pipeline.

By combining Rust-based utilities with Nextflow, each stage of the pipeline can be treated as a separate process—one step for Rust-based data transformation, another step for ABySS assembly, and further steps for QA/QC. Since Nextflow natively supports calling any command-line tool, integrating Rust binaries is seamless.

Below is a simplified example of a Nextflow script (abyss_pipeline.nf) that demonstrates how you might chain together a Rust-based pre-processing utility and an ABySS assembly step. This example assumes you have a Rust program called rust_preprocess for filtering reads, and that you have ABySS installed and accessible in your path. Adapt it as needed for your environment.

In this pipeline, the first stage, preprocessReads, applies a hypothetical Rust tool called rust_preprocess to filter or otherwise transform the raw FASTQ data. Next, the runAbyss process carries out the actual ABySS assembly, which can be configured to use either the single-processor or the MPI-enabled version, depending on how params.abyssExecutable is set. Finally, the assembleStats step uses a utility (in this example, abyss-fac) to generate assembly statistics for the resulting contigs, providing insights into metrics such as total assembled length and the distribution of contig sizes.

#### 2. Rust
Below is a minimal Rust program (rust_preprocess.rs) that demonstrates how you might parse command-line arguments (using the [structopt] or [clap] crate) and perform a rudimentary read-processing step. In practice, you would add sophisticated logic for quality trimming, adapter removal, or custom filtering.

This Rust tool can be compiled into a binary (cargo build --release), placed in your PATH, and then called by your Nextflow pipeline. You could extend this basic pattern to implement more advanced filtering—e.g., removing low-quality reads, trimming adapters, or performing k-mer-based classification—before passing the processed reads to ABySS.

When these technologies are combined, the result is a powerful, modular, and efficient workflow. Rust programs handle computation-heavy tasks such as read filtering with memory safety and high performance, while ABySS assembles the genome from short reads using either a single node or a multi-node MPI environment. Nextflow then orchestrates the entire pipeline, managing data ingestion, parallel execution, and the generation of reproducible reports on final assembly statistics.

This setup is well-suited to real-world scenarios where large volumes of sequencing data need to be processed reliably and quickly. By leveraging containers (e.g., Docker or Singularity images that include Rust, ABySS, and Nextflow), scientists and bioinformaticians can share their pipelines and deploy them consistently across different compute infrastructures—local servers, HPC clusters, or the cloud. The end result is a reproducible and scalable framework that can handle ever-increasing volumes of genomic data.

#### Project Structure:

```plaintext
experiment_3.6_3/
├── Cargo.toml                       # Rust project configuration and dependencies
└── src/
    ├── main.rs                      # Main Rust script containing program logic
    ├── main.nf                      # Nextflow workflow script
    ├── sample.fastq.rar             # Compressed FASTQ sample file
    ├── filtered_sample.fastq.rar    # Compressed filtered FASTQ sample file
    ├── output_nf.txt                # Nextflow output file
    ├── output_rs.txt                # Rust output file
    ├── nextflow.log.9               # Nextflow log file
    ├── assembly-1.dot               # Assembly graph in DOT format
    └── assembly-1.path              # Assembly path file
```

  
#### Cargo.toml

```toml
[package]
name = "rust_preprocess"
version = "0.1.0"
edition = "2021"

[dependencies]
clap = { version = "4.4.7", features = ["derive"] }
```

#### How to run:

##### 1) run main.rs in powershell:

```powershell
cargo run | tee output.txt
```

(run main.rs and save the output in output_rs.txt)

##### 2) run main.nf in powershell:

```powershell
cargo run main.nf | tee output_nf.txt
```

(run main.nf and save the output in output_rs.txt)


#### Explanation of the Output

##### 1. Rust Preprocessing Output (output_rs.txt)
* Preprocessing complete. Output saved to "filtered_sample.fastq"
* This confirms that the Rust script successfully processed the input FASTQ file.
* However, it does not indicate whether any reads were actually written to the output file.

#### 2. Nextflow Execution Output (output_nf.txt)
* The pipeline starts and executes preprocessReads, runAbyss, and assembleStats sequentially.
* preprocessReads completes successfully:

```rust
executor >  local (1)
[75/e55c40] preprocessReads (1) [100%] 1 of 1 ✔
```

* runAbyss encounters an error:
```rust
ERROR ~ Error executing process > 'runAbyss'
...
Process `runAbyss` terminated with an error exit status (2)
```

* The detailed error messages indicate:
```rust
assembly-1.fa:0: warning: file is empty

* The input file (filtered_sample.fastq) for ABySS wasn't empty and contain valid reads, the problem is from runAbyss.

* The process ultimately crashes with a segmentation fault:
```rust
make: *** [/home/trian/miniconda3/bin/abyss-pe.Makefile:568: assembly-2.dot] Segmentation fault (core dumped)
```
* This occurs because ABySS tries to process an empty file, causing memory access errors.

#### Conclusion

##### 1. The preprocessing step produced an empty output file.
* Possible reasons:
  * The input FASTQ file was already empty.
  * The Rust filtering logic removed all reads.
  * A file handling error occurred in the Rust script.
##### 2. ABySS failed due to an empty input file.
  * ABySS expects valid sequences to assemble. Without data, it cannot proceed.
  * This led to warnings about empty files and eventually a segmentation fault.
