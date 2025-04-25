## 6.7. Integrative Analyses with Rust-HTSlib

### experiment_67

The following Rust program integrates read coverage and variant annotation and it demonstrates how developers might handle concurrency, logging, and error handling while reading data from a BAM file with rust-htslib, parsing variants from a BCF file, and annotating each variant with gene information loaded from a GFF. Crates like ndarray or tch-rs can be included for advanced numerical or deep-learning tasks, while concurrency and safety are guided by Rust’s type system, preventing the data corruption issues common in large genomic pipelines.

This code begins by parsing command-line arguments using clap, enabling flexible configuration of the BAM file, BCF file, and GFF annotation path. The load_gff_annotations function demonstrates a simple method to read annotation lines into a HashMap, although real-world usage often calls for more elaborate data structures such as interval trees or suffix arrays for efficient querying.

After loading all variant records from the BCF file, the code uses rayon’s parallel iterator to process each variant concurrently, opening a new IndexedReader in process_variant to retrieve local coverage counts. Because Rust ensures immutable data by default, concurrency errors like race conditions are largely avoided. If a particular variant fails to process, the error is logged but does not terminate the entire pipeline, making the system more robust for large or imperfect datasets.

In a true production environment, developers often go further by adding HPC orchestration (e.g., via Nextflow or Snakemake), advanced numeric libraries for analyzing large coverage matrices (ndarray), and deep learning solutions (tch-rs) for tasks such as variant prioritization. Regardless of complexity, Rust’s memory and concurrency guarantees help maintain reliability and performance even for multi-terabyte genomic datasets spread across hundreds of nodes in a cluster.

In the following Nextflow scode, each pipeline stage runs in an ephemeral environment, whether on an HPC cluster or in the cloud, providing scalability and efficient use of resources. The snippet splits a BCF file by chromosome, processes each subset alongside BAM files and a GFF annotation, and then merges partial integration results.

In the splitBcf step, the large BCF file is indexed and then partitioned by chromosome. This approach leverages bcftools for subsetting, ensuring that each ephemeral container handles only a manageable slice of the data. The integrateData process invokes the Rust-based rust_integrate_tool, which follows the logic outlined in the previous Rust code. Each container receives one BCF subset, one BAM file, and the GFF annotation, enabling parallel analysis. Because Rust’s concurrency primitives ensure thread safety, many such containers can run simultaneously without risking data corruption.

Finally, the mergeIntegrations stage uses a hypothetical rust_merge_integration tool to combine partial JSON outputs into a single final_integration.json. This consolidated file can be further analyzed using Python, R, or additional Rust pipelines. For large-scale deployments, Nextflow can automatically spin up more containers based on the number of chromosomes or BAM files, while Rust’s performance characteristics keep memory usage and execution times predictable. This design helps maintain efficiency even when analyzing massive datasets in distributed computing environments.

#### Files contents:
* experiment_67/
  * Cargo.toml (Cargo.toml file for dependencies)
* experiment_66/src/
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
* experiment_67/target/debug/
  * rust_integrate_tool.rar (compressed rust_integrate_tool execution file output from running main.rs)
* experiment_67/src/work/5a/a130bb1fc0ec5c8c10aaeb3f5e1308/
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
