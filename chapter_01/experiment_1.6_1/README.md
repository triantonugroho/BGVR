## 1.6. Scientific Computation Workflow with Rust and Nextflow

### experiment_1.6_1

#### Nextflow Script:
Below is an example Nextflow DSL2 script demonstrating how to run a Rust read-trimming tool in parallel across multiple FASTQ files. Each file is processed by a separate job, and the pipeline remains agnostic to whether it executes on-premises under SLURM or in the cloud through other Nextflow executors.

In this example, Nextflow discovers all FASTQ files in the reads directory, sends each file to a process named RUST_TRIM, and manages parallel scheduling on the specified HPC system. By relying on container images for this process, the read-trimming tool and its dependencies remain consistent across nodes.

#### Rust Script:
The Rust trimmer itself can be a simple program that removes low-quality bases from both ends of a read. A short demonstration is shown here, but in a real application, it would include robust error handling, concurrency features, and thorough logging.

#### Project Structure:

```plaintext
experiment_1.6_1/
└── Cargo.toml                      # Rust project configuration and dependencies
experiment_1.6_1/src/
├── main.rs                         # Main Rust script containing program logic
├── main.nf                         # Nextflow workflow script
├── nextflow.log.9                  # Nextflow log file
├── output_nf.txt                   # Nextflow output file
├── output_rs.txt                   # Rust output file
└── reads/
    ├── sample.fastq.rar            # Compressed FASTQ sample file
    └── trimmed_sample.fastq.rar    # Compressed trimmed FASTQ sample file
```

#### How to run:

run in WSL:

```wsl
nextflow run main.nf | tee output.txt
```

(run main.nf and save the output in output.txt)
  
#### [dependencies]

```toml
[package]
name = "my_rust_trim"
version = "0.1.0"
edition = "2021"

[dependencies]

```

#### Explanation of the Output
The provided outputs (output_nf.txt and output_rs.txt) describe the execution results of the Rust-based FASTQ trimmer within a Nextflow pipeline. Below is an analysis of each output:

##### 1. output_rs.txt
This file contains a message indicating the successful completion of the Rust trimmer program:

```rust
Trimming complete: reads/sample.fastq -> reads/trimmed_sample.fastq
```

This message confirms that:

* The Rust script (main.rs) successfully processed the input FASTQ file (reads/sample.fastq).
* It trimmed the sequences and quality scores by removing 5 bases from both ends.
* The trimmed output was saved as reads/trimmed_sample.fastq.

##### 2. output_nf.txt
This file logs the execution of the Nextflow workflow, which runs the Rust-based trimming tool on multiple FASTQ files.

##### Workflow Execution Details

```nextflow
N E X T F L O W   ~  version 24.10.4
```

* Launching `main.nf` [nauseous_feynman] DSL2 - revision: 50463510c4
* This confirms that Nextflow version 24.10.4 was used.
* The workflow script main.nf was successfully launched.

Process Execution
```nextflow
executor >  local (2)
[7e/a51a5f] RUST_TRIM (2) [  0%] 0 of 2

executor >  local (2)
[6d/ef8b4f] RUST_TRIM (1) [100%] 2 of 2 ✔
```

* The process RUST_TRIM was executed locally instead of using SLURM.
* Two FASTQ files were processed.
* The task completed successfully (100% 2 of 2 ✔).

Output File Paths

```sh
/mnt/c/Users/trian/BGVR/chapter_02/experiment_26_1/src/work/7e/a51a5f8c64944739abd4e89d067c9d/trimmed_trimmed_sample.fastq
/mnt/c/Users/trian/BGVR/chapter_02/experiment_26_1/src/work/6d/ef8b4f966fad114a8651e7d27a67f1/trimmed_sample.fastq
```
* The trimmed FASTQ files were saved at specific locations.
* A possible issue is that one file is named trimmed_trimmed_sample.fastq, indicating that the trimming process may have been applied twice unintentionally.

Completion Summary

```nextflow
Completed at: 28-Feb-2025 10:54:13
Duration    : 1m 34s
CPU hours   : 0.1
Succeeded   : 2
```

* The workflow finished successfully in 1 minute and 34 seconds.
* 2 FASTQ files were processed without errors.

####Conclusion

##### 1. Successful Execution

* The Rust-based trimming script executed correctly and produced the expected outputs.
* The Nextflow pipeline successfully processed two FASTQ files.
* The process was executed locally without issues.

##### 2. Potential Issue: Duplicate Trimming

* The presence of trimmed_trimmed_sample.fastq suggests that the script might be running twice on the same file.
* This could be due to how files are passed to the Nextflow process.
* It is recommended to check main.nf and ensure that files are only trimmed once.

##### 3. Efficiency and Performance

* The pipeline completed in 1m 34s with 0.1 CPU hours, showing good performance for this dataset.
* The workflow is scalable and can handle more FASTQ files efficiently.

##### 4. Next Steps

* Debug why trimmed_trimmed_sample.fastq was created.
* Ensure input files are only passed once to the Rust script.
* If more files need processing, Nextflow can be extended to handle additional parallelization.

This setup demonstrates an efficient, automated, and reproducible FASTQ preprocessing pipeline using Rust for performance and Nextflow for workflow management. 










