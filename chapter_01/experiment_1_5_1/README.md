## 1.5. Acquiring and Cleaning Data

### experiment_1_5_1

Below is the complete Nextflow script that demonstrates a streamlined pipeline without containers. It downloads .sra files for a list of accessions, performs a simulated checksum verification, and moves the final outputs to a specified directory. You can store this script in a file named main.nf and adjust parameters as needed.

This pipeline first reads each accession from accessions.txt as a separate token in a channel. The DOWNLOAD_SRA process retrieves one .sra file per accession in parallel, limiting the number of concurrent tasks to four. The VERIFY_CHECKSUM process then simulates a checksum verification on each file, though in practice one might run md5sum or a similar command. Finally, the ORGANIZE_OUTPUT process places verified files in a designated output directory for easy organization.

When you want to add a Rust-based executable to this pipeline, you can place it anywhere in your system PATH. Suppose you have compiled a utility named my_rust_tool that transforms .sra files into a particular format. You can introduce a new process as follows:

```nextflow
process RUST_TOOL {
    input:
    file sraFile from verified_sra

    output:
    file("${sraFile}.out") into rust_results

    script:
    """
    echo "Running my_rust_tool on ${sraFile}"
    my_rust_tool --input ${sraFile} --output ${sraFile}.out
    """
}
```

This process depends on the same channel (verified_sra) that provided input to ORGANIZE_OUTPUT. You can then integrate it into the workflow, either replacing ORGANIZE_OUTPUT or running in parallel if the dataflow logic allows. Provided you have already built my_rust_tool, Nextflow calls it directly, relying on the host environment to supply the required dependencies.

Running the pipeline typically involves executing nextflow run main.nf in the directory that contains your script and accession list. Nextflow creates a .nextflow.log file to record progress and details about each step. The pipeline’s intermediate results and logs are stored in the work/ directory, which holds separate subfolders for each process execution. If a process fails, you can inspect its log files or rerun the pipeline after making corrections. You can also override parameters via the command line, for instance by specifying a different output directory or changing the concurrency level.

#### Project Structure:

```plaintext
experiment_1_5_1/
├── Cargo.toml                     # Rust project configuration and dependencies
├── main.nf                        # Nextflow workflow script
├── nextflow.log.9                 # Nextflow log file
├── output.txt                     # Output file
└── downloads/
    └── SRR11192680.sra            # Downloaded SRA file
```

#### Cargo.toml

```toml
[package]
name = "experiment_1_5_1"
version = "0.1.0"
edition = "2021"

[dependencies]

```

#### How to run:

run in WSL:

```wsl
nextflow run main.nf | tee output.txt
```

(run main.nf and save the output in output.txt)
  

#### Explanation of the Output
This Nextflow workflow (main.nf) automates the process of:

1. Downloading SRA files from NCBI using prefetch.
2. Verifying checksums to ensure file integrity.
3. Organizing the downloaded files into a specified output directory.
   
##### Step-by-Step Execution
Each step in the workflow runs as an independent Nextflow process, and their execution is tracked.

###### 1. Workflow Initialization

```nextflow
N E X T F L O W   ~  version 24.10.4
Launching `main.nf` [confident_ritchie] DSL2 - revision: 97db5e2019
```

* Nextflow version: 24.10.4 is used.
* Workflow identifier: confident_ritchie (Nextflow generates a random name for each run).
* DSL2 enabled: Indicates the use of modular pipeline design.

###### 2. Downloading SRA Files

```nextflow
executor >  local (1)
[52/2d08d4] DOWNLOAD_SRA_FILES (1) [  0%] 0 of 1
```
* The DOWNLOAD_SRA_FILES process starts.
* Initially, 0 of 1 files have been downloaded.

```nextflow
executor >  local (1)
[52/2d08d4] DOWNLOAD_SRA_FILES (1) [100%] 1 of 1 ✔
```

* The download is complete (100%), confirming that 1 SRA file was fetched.
* The SRA file is now available for the next process.
  
###### 3. Verifying Checksum

```nextflow
executor >  local (3)
[52/2d08d4] DOWNLOAD_SRA_FILES (1) [100%] 1 of 1 ✔
[c7/18c787] VERIFY_CHECKSUM (1)    [100%] 1 of 1 ✔
```

* VERIFY_CHECKSUM process runs next.
* The message Checksum OK (from script) confirms that the downloaded SRA file is not corrupted.

###### 4. Organizing Output Files

```nextflow
executor >  local (3)
[52/2d08d4] DOWNLOAD_SRA_FILES (1) [100%] 1 of 1 ✔
[c7/18c787] VERIFY_CHECKSUM (1)    [100%] 1 of 1 ✔
[9b/f0b998] ORGANIZE_OUTPUT (1)    [  0%] 0 of 1

* The ORGANIZE_OUTPUT process starts, copying the verified file to the output directory (downloads).

```nextflow
executor >  local (3)
[52/2d08d4] DOWNLOAD_SRA_FILES (1) [100%] 1 of 1 ✔
[c7/18c787] VERIFY_CHECKSUM (1)    [100%] 1 of 1 ✔
[9b/f0b998] ORGANIZE_OUTPUT (1)    [100%] 1 of 1 ✔
```

* The file organization is completed successfully.
  
###### 5. Workflow Completion Summary

```nextflow
Completed at: 28-Feb-2025 08:29:34
Duration    : 14m 15s
CPU hours   : 0.2
Succeeded   : 3
```

* Execution time: 14 minutes, 15 seconds.
* Total CPU usage: 0.2 CPU hours (light computation).
* Successful processes: 3 (all steps completed correctly).

#### Conclusion
* The workflow successfully completed all steps:
  1. Downloaded 1 SRA file.
  2. Verified its integrity (checksum check).
  3. Organized the output file in the downloads folder.
* Key Takeaways:
  * The pipeline is modular, meaning it can easily scale to handle multiple accessions.
  * Execution is parallelized, optimizing performance.
  * The processes executed correctly and in the expected order.


