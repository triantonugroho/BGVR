## 6.5. Advanced Data Structures for HTS Analysis

### experiment_65

The following Rust program uses an interval tree for genomic coverage queries in Rust. It includes a command-line interface, logging, and concurrency considerations, illustrating how developers can integrate this data structure into HPC or cloud-based pipelines. The tree construction sorts intervals, identifies a center, and divides intervals into left and right subtrees, thereby enabling efficient queries for overlapping intervals. Once the tree is built, queries can safely be performed in parallel thanks to Rust’s immutability by default, eliminating the need for manual synchronization.

In this example, each interval includes a coverage value, but you can store any additional data or metadata as needed, such as read depth or variant calling statistics. The IntervalTree struct partitions intervals around a chosen center, placing strictly smaller intervals to the left and larger ones to the right. Intervals that overlap the center remain in the current node, avoiding duplication. This design significantly reduces query times, as the search visits only those subtrees potentially intersecting with the query range.

Once the tree is built, queries are read-only, so you can safely run multiple queries in parallel. This concurrency is handled by rayon’s parallel iterators, which create a thread pool and distribute the work evenly. The immutability of the interval tree ensures that no locks or synchronization overhead is needed for queries, even in large-scale HPC or cloud environments. In more advanced pipelines, you might integrate crate features like rust-htslib to stream intervals directly from BAM or CRAM files, tch-rs for deep learning models that predict coverage anomalies, or ndarray for vectorized numerical transformations.

Because Rust produces a statically linked binary with minimal runtime overhead, this approach is ideal for container-based deployments on HPC clusters or in the cloud. Each container can construct and query interval trees for a portion of the data, then shut down, all under the control of orchestrators like Nextflow or Kubernetes. By relying on Rust’s safety guarantees, you reduce the likelihood of concurrency bugs and segmentation faults, which become especially important for large genomic datasets that can run into the hundreds of gigabytes or terabytes.

The Nextflow workflow below aligned with the previously discussed Rust code that creates and queries interval trees for coverage data. Each process corresponds to a distinct step in the pipeline, running in an ephemeral container that spins up, performs the task, and shuts down automatically. This approach scales efficiently on HPC clusters or cloud platforms, especially for large genomic datasets where batch processing and parallel querying are critical to minimizing runtime.

In this workflow, the coverageComputation step runs rust_coverage_tool on each BAM file, leveraging Rust’s concurrency features (via rayon or similar) to compute coverage intervals quickly. The mergeCoverage step consolidates all coverage TSV files, preparing a single dataset for interval queries. Finally, the intervalQuery step invokes rust_interval_query_tool, which reads the merged file, constructs an interval tree, and answers coverage queries in parallel.

Such ephemeral task orchestration saves both time and computational resources. Whenever additional nodes or containers become available, Nextflow can launch more coverageComputation tasks concurrently, making it straightforward to handle large cohorts. On HPC clusters, Nextflow can interface with schedulers like SLURM, PBS, or LSF. In cloud environments, it can scale dynamically on platforms like AWS Batch or Kubernetes. Developers often wrap these Rust tools in Docker or Singularity containers, ensuring a consistent runtime environment (including the Rust standard library, rust-htslib dependencies, or any required GPU libraries for tch-rs). Because Rust compiles to a single statically linked binary, containers can remain minimal, leading to faster spin-up times and lower overhead when dealing with high-throughput genomic workflows.

Many AI engineers in pharmaceutical R&D apply these Rust workflows to large-scale population studies, often investigating how certain intervals or variants correlate with disease risk or drug response. By combining ephemeral container usage with well-crafted interval trees or graph-based references, they can analyze multi-terabyte genomic data sets in a fraction of the time once required by more conventional scripting. Reports from major consortia confirm that dynamic load balancing and concurrency safety prevent the data corruption or performance bottlenecks sometimes seen with more ad-hoc approaches (Di Tommaso et al. (2017)).

#### Files contents:
* experiment_65/
  * Cargo.toml (Cargo.toml file for dependencies)
* experiment_65/src/
  * main.rs (rust script)
  * main.nf (nextflow script)
  * bams.txt (bams dataset list in text file)
  * genome_intervals.txt (genome intervals in text file)
  * sample1.bam (sample 1 bam file)
  * sample2.bam (sample 2 bam file)
  * output.txt (text file output)
* experiment_65/target/debug/
  * interval_query_tool.rar (compressed interval_query_tool execution file output from running main.rs)
* experiment_64/src/work/1b/90098e3428a3e74a207b0d7487cc61/
  * coverage_sample2.bam.tsv (tsv coverage sample 2 output file)
* experiment_64/src/work/1e/e68580ed5054b0a8c2004c5451626c/
  * merged_coverage.tsv (tsv merged coverage output file)
  * query_result_50-60.tsv (tsv query result region 50-60 output file)
* experiment_64/src/work/5e/4032832140ba89a1a6f61fd225a275/
  * merged_coverage.tsv (tsv merged coverage output file)
  * query_result_40-70.tsv (tsv query result region 50-60 output file)
* experiment_64/src/work/20/8f63f34a73e7b282c8b7ceddf3a14a/
  * coverage_sample1.bam.tsv (tsv coverage sample 1 output file) 
* experiment_64/src/work/63/4004f9ab588f91e2b3c2b0cf759ccf/
  * merged_coverage.tsv (tsv merged coverage output file)
  * query_result_40-70.tsv (tsv query result region 40-70 output file)
* experiment_64/src/work/b9/c791cf844bb5d13891840c0ad99942/
  * merged_coverage.tsv (tsv merged coverage output file)
  * query_result_20-40.tsv (tsv query result region 20-40 output file)
* experiment_64/src/work/bc/2072d2366055ffe43d52fc4fd916ef/
  * merged_coverage.tsv (tsv merged coverage output file)
  * query_result_1-10.tsv (tsv query result region 1-10 output file) 
* experiment_64/src/work/e3/68fcb7fdfa6d9564a7f1799986a01b/
  * merged_coverage.tsv (tsv merged coverage output file)
  * query_result_0-5.tsv (tsv query result region 0-5 output file) 
* experiment_64/src/work/ed/a4a1d8967a56e3e735289c62ec5c37/
  * merged_coverage.tsv (tsv merged coverage output file)
  * coverage_sample1.bam.tsv (tsv coverage sample 1 output file) 
  * coverage_sample2.bam.tsv (tsv coverage sample 2 output file) 
* experiment_64/src/work/f4/04206d570d36efc492253ba1bffe2f/
  * merged_coverage.tsv (tsv merged coverage output file)
  * query_result_5-25.tsv (tsv query result region 5-25 output file)

#### How to run:

run main.rs in wsl:

```wsl
cargo run -- \
  --intervals 1-10:3 5-25:7 20-40:6 30-32:4 50-60:1 0-5:2 40-70:5 \
  --queries 1-10 5-25 20-40 30-32 50-60 0-5 40-70 | tee output.txt
```

(run main.rs with intervals and queries parameters and save the output in output.txt)

run main.nf in wsl:

```wsl
nextflow run main.nf
```

run main.nf with this parameters:
params.bam_list       = 'bams.txt'
params.ref_intervals  = 'genome_intervals.txt'
params.parallel_chunck_size = 50000 
params.mock = true 

#### [dependencies]

```toml
anyhow = "1.0"
clap = { version = "4.4", features = ["derive"] }
env_logger = "0.11.8"
log = "0.4"
rayon = "1.7"
```

#### Explanation of the Output
The output directory structure and files shown reflect the following stages in your Nextflow pipeline:

##### 1. Coverage Computation for Each BAM File (coverage_sampleX.bam.tsv):

* The pipeline starts by processing each BAM file (e.g., coverage_sample1.bam.tsv, coverage_sample2.bam.tsv).

* For each BAM file, the coverageComputation process generates a corresponding coverage file in TSV format. This file contains genomic coverage data for different regions (e.g., chromosomes, start and end positions) along with the depth of coverage.

* These coverage files represent the genomic regions covered by reads from the BAM files.

##### 2. Merging Coverage Files (merged_coverage.tsv):

* The mergeCoverage process combines the coverage files generated in the previous step. It outputs a single merged_coverage.tsv file, which consolidates coverage data across all input BAM files.

* This file serves as the base for subsequent querying of intervals.

##### 3. Interval Querying (query_result_X-Y.tsv):

* The intervalQuery process runs the interval query tool (implemented in Rust) on the merged coverage data, performing queries for different genomic regions (e.g., query_result_50-60.tsv, query_result_40-70.tsv).

* The queries are based on the intervals specified in genome_intervals.txt. For each query, the tool searches for overlapping intervals within the merged coverage file.

* Each query result file (query_result_X-Y.tsv) contains the genomic regions (chromosomes, start, end, coverage) from the merged coverage data that overlap with the query interval (e.g., 50-60, 40-70, etc.).

#### Summary of the Output:
* The output directories correspond to different stages of the workflow, with each step having its own unique files:

  * Coverage Files (coverage_sampleX.bam.tsv): Represent the computed coverage data for each BAM file.

  * Merged Coverage (merged_coverage.tsv): The consolidated coverage data across all BAM files.

  * Query Results (query_result_X-Y.tsv): The results of the interval queries, showing the overlaps with specific genomic ranges.

Here’s a breakdown of the output files:

* Experiment 64 has several subdirectories that correspond to the different stages of the workflow.

  * Each subdirectory contains the output files from coverage computation, merging, and interval querying.

  * The query result files (query_result_X-Y.tsv) contain the intervals that match the query ranges and overlap with the merged coverage.

For example, the file query_result_50-60.tsv contains the results of querying the merged coverage file for intervals that overlap with the range 50-60. Similarly, the other query_result_X-Y.tsv files contain the results for other genomic intervals (e.g., 1-10, 20-40, etc.).

#### Conclusion:
* The pipeline works as expected by processing BAM files to compute coverage, merging the coverage data, and running interval queries on the merged data.

* Each output directory corresponds to a specific task in the pipeline, and the generated files provide useful data for genomic analysis.

  * Coverage files are essential for understanding how different regions of the genome are covered by reads in the input BAM files.

  * Merged coverage combines these files into a single dataset that facilitates easier querying and analysis.

  * Query results provide the intervals from the merged coverage that match the genomic intervals specified in the queries, giving insight into the regions of the genome that were covered at different depths.

* This approach efficiently handles large genomic data sets by leveraging Nextflow's parallelism and the Rust-based interval tree querying tool, ensuring scalability and fast querying.

You can now use these results to further analyze coverage patterns, genomic region overlaps, or any other bioinformatics tasks related to the coverage data.
