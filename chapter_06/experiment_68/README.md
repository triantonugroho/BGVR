## 6.8. Summary and Future Directions

### experiment_68

The following Rust program demonstrates how to create a command-line binary for computing coverage in a BAM file. It uses clap (successor to structopt) for argument parsing, anyhow for robust error handling, and rust-htslib for BAM/CRAM indexing and iteration. This design is ideal for HPC or cloud pipelines orchestrated by Nextflow. Each workflow step can be containerized with Docker or Singularity to guarantee reproducible environments, ensuring that the Rust toolchain and associated system libraries remain consistent across clusters or cloud providers.

This revised tool reads the user-specified BAM file and region, calculates the number of aligned reads in that region, and writes the final coverage count to a text file. It relies on indexed access, so a corresponding BAI file is necessary for random queries. The usage of anyhow and env_logger ensures that errors are contextualized and logged appropriately, which aids debugging in large-scale HPC or cloud runs. When tasks are parallelized—either by chunking the genome into multiple regions or distributing sample sets across HPC nodes—each ephemeral container or compute node runs this binary and writes its partial result to a file or object store, minimizing memory usage and total runtime.

If more advanced machine learning operations are needed, the tch-rs crate can inject PyTorch-based neural networks directly into the workflow, while ndarray enables vectorized numeric computations on coverage matrices. Projects that require on-the-fly data querying can integrate libraries like polar or polars for structured queries, and linfa for classical machine learning algorithms. Regardless of additional complexity, Rust’s strong type system and concurrency model help keep large genomic pipelines reliable and maintainable, while containerization ensures consistent environments across HPC clusters or cloud platforms.

In the following Nextflow code, each task processes one BAM file and one genomic region in an ephemeral container, writes a partial coverage result, and then merges these partial outputs in a final stage. This design is typical for HPC or cloud pipelines managed by Nextflow, ensuring each container only handles a subset of the data while logging and versioning can be standardized for robust industrial usage.

In this workflow, the coverageCalc process uses the rust_coverage_tool compiled from the code shown previously. Each ephemeral container processes a unique BAM file and genomic region, writing coverage results to coverage_<BAM>_<REGION>.txt. Nextflow’s channel mechanics collect the outputs, passing them to mergeCoverage, which concatenates all partial results into merged_coverage.txt. This pattern allows high-throughput computation on large cohorts or extensive genomic intervals, as Nextflow can scale out tasks onto multiple HPC nodes or cloud instances.

In production pipelines, one often adds more robust error handling (e.g., retries for transient failures), logging strategies (potentially with structured logs that can be aggregated by tools like the ELK stack), and versioning (ensuring that each container runs a fixed version of Rust and the coverage tool). Such practices foster reproducibility and traceability, whether you are working on a single workstation or a massive cloud cluster.

#### Files contents:
* experiment_68/
  * Cargo.toml (Cargo.toml file for dependencies)
* experiment_68/src/
  * main.rs (rust script)
  * main.nf (nextflow script)
  * ref.fa (reference fasta file input)
  * ref.fa.fai (indexed ref.fa)
  * annotations.gff (annotation gff file input)
  * bams.txt (text file contain bam file name list)
  * test.vcf (test vcf file)
  * test1.bam (test 1 bam file input file)
  * test1.bam.bai (indexed test1.bam)
  * test1.sam (sam file to make test1.bam file)
  * test1.sorted.bam.bai (sorted indexed test1.bam file)
  * output.txt (text file output)
* experiment_68/target/debug/
  * rust_integrate_tool.rar (compressed rust_integrate_tool execution file output from running main.rs)
* experiment_68/src/work/5a/a130bb1fc0ec5c8c10aaeb3f5e1308/
  * integrated_test1.bam_split_chr1.bcf.json (integrated tes1.bam_split_chr1.bcf json output file)

#### How to run:

run main.rs in wsl:

```wsl
cargo run -- --bam /mnt/c/Users/trian/BGVR/chapter_06/experiment_67/src/test1.bam --bcf /mnt/c/Users/trian/BGVR/chapter_06/experiment_67/src/cohort.bcf --gff /mnt/c/Users/trian/BGVR/chapter_06/experiment_67/src/annotations.gff | tee output.txt
```

(run main.rs with test1.bam, cohort.bcf and annotations.gff as input parameter and save the output in output.txt)

run main.nf in wsl:

```wsl
nextflow run main.nf
```

run main.nf with this parameters:
params.bam_list = 'bams.txt'
params.bcf_file = 'cohort.bcf'
params.gff_file = 'annotations.gff'
params.bam_dir = '.' // Default to current directory

#### [dependencies]

```toml
anyhow = "1.0"
clap = { version = "4.4", features = ["derive"] }
env_logger = "0.11.8"
log = "0.4"
rayon = "1.8"
rust-htslib = "0.49.0"
```

#### ✅ Explanation of Output and Workflow Execution
Workflow executed successfully in Nextflow DSL2, and each process produced the expected outputs. Let's walk through what happened and explain the final result.

##### 🧩 Workflow Breakdown
###### 1. splitBcf Process

* Input: cohort.bcf
* Action:
  * Indexed the BCF file.
  * Extracted each chromosome separately using bcftools view.
* Output: Files like split_chr1.bcf.

###### 2. integrateData Process

* Inputs:
  * split_chr1.bcf (from splitBcf)
  * test1.bam (from bams.txt)
  * annotations.gff
* Action:
  * Ran your Rust tool rust_integrate_tool with the provided BAM, BCF, and GFF.
  * Inside the tool:
    * Parsed GFF annotations into a lookup table.
    * Read all variants from the BCF file.
    * For each variant:
      * Queried the BAM file for reads overlapping the variant position to compute coverage.
      * Matched the variant’s position against annotated regions in the GFF.
   * Produced a JSON output summarizing coverage and annotations per variant.

* Output:
integrated_test1.bam_split_chr1.bcf.json, containing this:

```txt
Loading GFF annotations from annotations.gff
Starting integrative analysis on BAM: test1.bam, BCF: split_chr1.bcf
Processed 1 variants.
Variant: chr1:4 ref=A alt=T coverage=1 annotation=Some("ID=gene1;Name=GeneA")
```

✔️ This confirms that:
* The variant at chr1:4 had 1 read covering it.
* It matched an annotated region: "ID=gene1;Name=GeneA".

###### 3. mergeIntegrations Process
* Input: One or more integrated_*.json files
* Action:
  * Merged the files into a single JSON array.
* Output: final_integration.json — in this case, just wrapping the single JSON object in brackets [...].

##### 📘 Final Output Structure
The working directory for integration:

```psql
work/5a/a130bb1fc0ec5c8c10aaeb3f5e1308/
├── integrated_test1.bam_split_chr1.bcf.json
├── split_chr1.bcf
├── annotations.gff
├── test1.bam
├── test1.bam.bai
```

And merged result:

```text
final_integration.json  → contains:
[
  {
    chrom: "chr1",
    pos: 4,
    ref_allele: "A",
    alt_allele: "T",
    coverage: 1,
    annotation: "ID=gene1;Name=GeneA"
  }
]
```

#### ✅ Conclusion
🔄 The pipeline correctly split the BCF, indexed the BAM, annotated and calculated coverage using the Rust tool.

📦 The final JSON output is suitable for downstream analysis, such as visualization or machine learning input.

🧪 This pipeline is now modular and reproducible for any combination of BAM + BCF + GFF files.

