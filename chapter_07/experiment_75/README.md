## 7.5. Integrating Rust Noodles into Nextflow Pipelines

### experiment_75

This Rust code example demonstrates how to open and index a BAM file if necessary, calculate base‐by‐base coverage over a specific genomic region using the noodles-bam and noodles-core crates, and finally serialize the results to JSON via serde. It showcases how to integrate multiple crates to perform essential bioinformatics tasks in a performant and maintainable manner, laying the groundwork for more advanced genomic operations.

In a single pass, main checks for an existing BAM index and creates one if it does not exist, then calls compute_coverage to iterate over all reads that overlap the region. That function calculates per‐base coverage by determining the overlap between each read and the query range, accumulating results into a vector of coverage counts that are serialized via serde_json for easy output or further downstream processing.

This code snippet ensures that the given BAM file has an index, then reads the specified genomic region to count how many reads fall into that range. The coverage count is then written to a JSON file, which is suitable for ingestion by downstream analytics or reporting steps. To further enhance industrial-scale robustness, logs can be collected via crates like tracing for structured logging, and the code could be encapsulated in Docker or Singularity containers with pinned versions to guarantee reproducibility.

Below is a minimal Nextflow script that orchestrates multiple coverage tasks in ephemeral containers. Each container runs the above Rust tool on different inputs, illustrating how Nextflow’s DAG scheduling synergizes with Rust’s concurrency. As with the Rust code, everything is enclosed in a panel for direct usage, and each ephemeral container can run the coverage computation in parallel, merging the JSON outputs in a final step.

Below is a minimal Nextflow script demonstrating how to run the previously described Rust coverage tool inside ephemeral containers on multiple BAM inputs, with each coverage task executed in parallel and a final step merging the JSON outputs. This setup showcases how Nextflow’s DAG scheduling seamlessly integrates with Rust’s concurrency, allowing each container‐based job to process a specific input and produce JSON results for subsequent merging.

The RUN_COVERAGE process receives tuples that include a sample identifier along with its BAM and BAI files, then spins up ephemeral containers (one per sample) running the Rust tool to compute coverage. Each container outputs a JSON file, which is collectively consumed by the MERGE_COVERAGE process. In the final step, a simple JSON merge combines these separate coverage results into a single JSON, illustrating how Nextflow’s dataflow model and Rust’s concurrency align to perform scalable, reproducible genomics analysis in containerized environments.

In this Nextflow script, coverageCalc is invoked once for each combination of BAM file and genomic region. Nextflow takes care of scheduling these tasks on available compute resources, whether those are local cores, HPC cluster nodes, or cloud instances. The ephemeral containers each run the Rust coverage tool, generating partial JSON files that are then merged in the final stage of the pipeline. This pattern, known as scatter-gather, exemplifies a standard approach for handling large genomic datasets by splitting them into manageable subsets and recombining the results.

#### Files contents:
* experiment_75/
  * Cargo.toml (Cargo.toml file for dependencies)
* experiment_75/src/
  * main.rs (rust script)
  * main.nf (nextflow script)
  * merged.vcf (merged vcf file as output file after running main.rs)
  * output.json (output json file)
  * sample1.vcf (sample 1 vcf file as input file)
  * sample2.vcf (sample 2 vcf file as input file)
  * vcf_list.txt (text file contain vcf file list) 
  * output.txt (text file output)
* experiment_75/src/results/
  * merged_vcf.bcf (merged vcf bcf file as output file after running main.nf)
  * pipeline_report.html (pipeline report html file as output file after running main.nf)
* experiment_75/src/work/0a/2d77eea602eedc6fe33279113344e1/
  * pipeline_report.html (pipeline report html file as output file after running main.nf)
* experiment_75/src/work/85/5ef3143a9d19c98e673a43cae7e651/
  * local_vcf_list.txt (text file contain local vcf list)
  * merged_vcf.bcf (merged vcf bcf file as output file after running main.nf)
* experiment_75/target/debug/
  * rust_vcf_merge_tool.rar (compressed rust_vcf_merge_tool execution file output from running main.rs)

#### How to run:

run main.rs in wsl:

```wsl
cargo run -- --vcf-list vcf_list.txt --out merged.vcf --threads 4 --format bcf | tee output.txt
```

(run main.rs with vcf_list.txt using 4 threads as input parameter and merged.vcf as output file and formatted to bcf file and save the text output in output.txt)

run main.nf in wsl:

```wsl
nextflow run main.nf
```

run main.nf with this parameters:
params.sample_list = 'vcf_list.txt'
params.output_vcf = 'merged.vcf'
params.threads = 4
params.format = 'bcf'
params.tool_path = "/mnt/c/Users/trian/BGVR/chapter_07/experiment_74/target/debug/rust_vcf_merge_tool"

#### [dependencies]

```toml
anyhow = "1.0"
clap = { version = "4.0", features = ["derive"] }
rayon = "1.5"
rust-htslib = "0.49"
env_logger = "0.11.8"
log = "0.4"
```

#### Explanation of the Output
##### 🦀 main.rs (Rust): Output Explanation
This file contains a Rust program designed to merge multiple VCF files into a single, multi-sample VCF or BCF. When executed, it produces:

###### 🔹 merged.vcf or merged.bcf
* Created by the rust_vcf_merge_tool.
* The merged output of all VCF files listed in vcf_list.txt.
* Stored in the working directory or wherever the --out parameter points to.
* The example you gave used BCF format (--format bcf), resulting in:

```text
merged_vcf.bcf
```

###### 🔹 output.txt (stdout)
* Contains a log message from the tool:

```text
Merge completed in 28 ms
```

* Shows successful completion and performance timing.

##### 🧬 main.nf (Nextflow): Output Explanation
This file defines a Nextflow pipeline that automates the execution of the Rust merging tool and then generates a simple HTML report.

###### 🔹 Workflow Breakdown
###### ✅ Step 1: mergeVCF process
* Takes:
  * vcf_list.txt (VCF paths)
  * All .vcf files in the same directory
* Writes local_vcf_list.txt (a new list with local paths) to ensure compatibility with the Rust tool.
  * Runs the compiled Rust binary rust_vcf_merge_tool
  * Produces:
    * merged_vcf.bcf

###### ✅ Step 2: generateReport process
* Takes:
  * merged_vcf.bcf from the previous step
  * Runs a reporting command (simulated with fallback HTML generation)
* Produces:
  * pipeline_report.html

###### 📂 Output Directory: /mnt/c/Users/trian/BGVR/chapter_07/experiment_74/src/results/
You will find:
* merged_vcf.bcf – The merged variant file
* pipeline_report.html – A simple HTML summary

```html
<html><body><h1>Report for merged_vcf.bcf</h1><p>Processing complete.</p></body></html>
```

#### ✅ Conclusion
You have successfully:

* 📦 Written a Rust CLI tool to efficiently merge VCF files.

* 🔄 Wrapped the tool in a Nextflow pipeline to automate and scale the merging + reporting process.

* 🧪 Verified correct functionality with structured output:

  * The merge was completed (merged_vcf.bcf)
  * A simple report was generated (pipeline_report.html)

This setup is now modular, reproducible, and ready to handle more complex VCF workflows or plug into larger bioinformatics pipelines.

