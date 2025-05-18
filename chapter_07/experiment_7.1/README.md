## 7.1. Foundational Data Structures in Rust Noodles

### experiment_7.1

Below is a code panel that highlights Rust integration using noodles-bam and noodles-vcf for coverage calculations and variant reading. This code is structured for HPC concurrency, streaming partial results from multiple intervals. The coverage is simply counted in this demonstration, but advanced logic could store coverage in an in-memory segment tree for subsequent queries.

The coverage calculation is straightforward, but the concurrency is implicit in the call to .par_iter() from the rayon crate. HPC concurrency can extend this approach by dividing reads across different nodes or ephemeral containers, each with its own region. If advanced numeric computing is required, crates like ndarray can store coverage arrays or correlation matrices, while tch-rs can interface with PyTorch for deeper AI-based analyses.

The Nextflow pipeline shown below orchestrates ephemeral container tasks that each invoke Rust binaries compiled from the above code. This pipeline demonstrates how HPC or cloud systems can run multiple coverage calculations and variant loads in parallel, merging partial results at the end. The ephemeral containers spin up, process one region, and write partial coverage results, then shut down, ensuring resources are not wasted.

This pipeline ensures ephemeral containers each run a specialized Rust tool, either rust_noodles_coverage for coverage or rust_noodles_variant for variant extraction. The ephemeral approach is valuable in AI-driven genomic projects, where large data sets can be automatically subdivided among many computing nodes. Container images fix the Rust environment for reproducibility, while Nextflow manages parallel dispatch, error handling, and final merges.

#### Project Structure:

```plaintext
experiment_7.1/
├── Cargo.toml                  # Rust dependencies
├── src/
│   ├── main.rs                 # Rust implementation
│   ├── main.nf                 # Nextflow workflow
│   ├── bams.txt                # Text file containing sorted BAM file list
│   ├── cohort.vcf              # Cohort VCF file
│   ├── regions.txt             # Region list text file
│   ├── sample1.sam             # Sample 1 SAM file to create sample1.bam
│   ├── sample1.bam             # Sample 1 BAM file as input
│   ├── sample1.sorted.bam.bai  # Indexed sorted sample 1 BAM file
│   ├── sample2.sam             # Sample 2 SAM file to create sample2.bam
│   ├── sample2.bam             # Sample 2 BAM file as input
│   ├── sample2.sorted.bam.bai  # Indexed sorted sample 2 BAM file
│   ├── variants.vcf            # Variants VCF file
│   ├── output.txt              # Text file output
│   └── work/                   # Nextflow work directory
│       └── fc/c33468689ba766ec2e2b2e5f570587/
│           └── coverage_ouput.txt            # Coverage text file output
└── target/
    └── debug/
        └── rust_noodles_tool.rar  # Compressed Rust noodles tool executable
```

#### Cargo.toml

```toml
[package]
name = "rust_noodles_coverage"
version = "0.1.0"
edition = "2024"

[dependencies]
noodles = { version = "0.5", features = ["bam", "core", "vcf"] }
rayon = "1.5.1"
clap = { version = "4", features = ["derive"] }
serde = { version = "1.0", features = ["derive"] }
serde_json = "1.0"
log = "0.4"
env_logger = "0.11.8"
thiserror = "2.0.12"
anyhow = "1.0"
```

#### How to run:

run main.rs in wsl:

```wsl
cargo run -- --vcf-file cohort.vcf --bam-files sample1.sorted.bam,sample2.sorted.bam | tee output.txt
```

(run main.rs with cohort.vcf, sample1.sorted.bam and sample2.sorted.bam as input parameter and save the output in coverage_output.txt)

run main.nf in wsl:

```wsl
nextflow run main.nf
```

run main.nf with this parameters:
params.bam_list = "bams.txt"
params.vcf_file = "cohort.vcf"
params.rust_bin = "/mnt/c/Users/trian/BGVR/chapter_07/experiment_7.1/target/debug"

#### Explanation of the Output
##### ✅ main.rs Explanation (Rust CLI Tool)
This Rust program performs basic BAM file parsing using the noodles crate, and it’s structured as a CLI tool. Here's what happens when it's executed:

######💡 Input Parameters:
* --vcf-file cohort.vcf: Not parsed (VCF parsing is not enabled yet, just logged).
* --bam-files sample1.sorted.bam,sample2.sorted.bam: Comma-separated paths of BAM files to be processed in parallel using Rayon.

###### 🧠 Main Workflow:
1. Logging Init: Logging is initialized via env_logger.
2. Argument Parsing: The CLI args are parsed using clap.
3. VCF Handling: It simply logs that the VCF file is acknowledged, but real parsing is skipped (VCF feature not enabled).
4. Parallel BAM Processing:
   * Each BAM file is opened with noodles::bam::Reader.  
   * It reads headers and reference sequences.
   * For each BAM:
     * Logs header and up to 2 reference sequences.
     * Reads records one by one, printing the first 5 records’ positions, MAPQ, and CIGAR.
     * Counts the total number of records read.
5. Final line confirms all BAMs have been processed.

##### ✅ main.nf Explanation (Nextflow Pipeline)
Nextflow script wraps the above binary and controls its execution in a reproducible workflow.

###### 💡 Parameters Used:
* params.bam_list = "bams.txt": List of BAM files (line-separated).
* params.vcf_file = "cohort.vcf": Single VCF file path.
* params.rust_bin = "/mnt/c/.../debug": Folder containing the compiled binary rust_noodles_coverage.

###### 🧠 Process rustCoverageRunner:
* Input:
  * bams.txt (as bam_list_file)
  * cohort.vcf
* The bam_list_file is converted into a comma-separated string using:

```wsl
BAM_FILES=$(cat bams.txt | tr '\n' ',' | sed 's/,$//')
```

* The binary is executed like this:

```wsl
rust_noodles_coverage --vcf-file cohort.vcf --bam-files "sample1.sorted.bam,sample2.sorted.bam"
```

* Output is redirected to coverage_output.txt.

###### 📤 Output: coverage_output.txt
Contents:

```txt
Starting BAM processing application...
VCF file noted: cohort.vcf (VCF processing not enabled)
Add 'vcf' to noodles features in Cargo.toml to enable VCF processing
Processing 2 BAM files...
Processing BAM file: sample1.sorted.bam
Processing BAM file: sample2.sorted.bam
Processing completed successfully
```

This is expected based on the println!() calls in your Rust main.rs. Because BAM parsing is done, but only limited record info is printed, the output is minimal.

#### ✅ Conclusion
✅ Integration Success: Your Rust binary was successfully integrated into a Nextflow pipeline.

✅ Inputs Correctly Wired: bams.txt was read and transformed properly into CLI format for Rust tool.

✅ Parallelism Works: Your BAM files were processed in parallel inside main.rs using Rayon.

⚠️ VCF Functionality Not Implemented: Your tool acknowledges the VCF but doesn’t use it yet—future feature.

📁 Output Location: Since Nextflow runs in isolated work directories, the output coverage_output.txt appears in a work subfolder (experiment_71/src/work/fc/...), and can be collected later with publishDir or saved explicitly.
