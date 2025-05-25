## 9.4 Differential Expression Analysis

### experiment_9_4

Differential expression analysis in a production environment typically begins with a normalized count matrix. The following Rust code illustrates how one might compute a simplistic negative binomial test for each gene, then apply Benjamini-Hochberg correction across all tests. This example uses ndarray for matrix operations, serde for data input, and statrs for statistical distributions.

Although this example aims to illustrate a fundamental approach, genuine production environments use more sophisticated statistical techniques. The statrs crate can implement various distributions, and the negative binomial approach would typically involve maximum likelihood estimation rather than simplistic placeholders. Additional crates such as argmin or ncvx might be employed for robust parameter estimation. Logging frameworks and concurrency patterns (e.g., via rayon) also ensure scalability and fault tolerance for large experiments.

Below is a Nextflow script that orchestrates the differential expression step. It is intended to follow earlier pipeline stages such as normalization or batch correction, reflecting a typical analysis flow.

This pipeline demonstrates how Nextflow automates the read-process-write sequence, ensuring that each step is reproducible and containerized for consistent execution. For industrial-scale usage, one might integrate advanced reporting features, environment modules, or distributed computing backends. The synergy of Rust’s performance and type safety with Nextflow’s orchestration capabilities is well-suited to large bioinformatics workflows.

In real-world pipelines, AI engineers and bioinformaticians often tailor the choice of model, statistical parameters, and coding patterns to the needs of specific projects. Teams may integrate advanced shrinkage methods or Bayesian frameworks to achieve higher confidence in small-sample scenarios. A global pharmaceutical consortium recently reported a dramatic reduction in computational overhead when swapping an older R-based differential expression step for a Rust-based approach orchestrated via Nextflow. The streamlined pipeline allowed them to pivot quickly to new hypotheses, showing how engineering optimizations can directly impact scientific progress in drug discovery and precision medicine.

#### Project Structure:

```plaintext
experiment_9_4/
├── Cargo.toml                              # Rust package configuration and dependencies
├── main.nf                                 # Nextflow pipeline script
├── generate_sample_data.py                 # Python script to generate sample data
├── rnaseq-normalizer                       # Compiled Rust executable binary
├── README.md                               # Project documentation
├── nextflow.config                         # Nextflow configuration (optional)
│
├── src/                                    # Rust source code
│   └── main.rs                            # Main Rust normalization implementation
│
├── target/                                 # Rust build artifacts
│   └── release/
│       └── rnaseq-normalizer              # Compiled Rust executable binary
│
├── test_data_small/                        # Small test dataset (100 genes × 6 samples)
│   ├── raw_counts.tsv                     # Raw count data
│   ├── batch_metadata.tsv                 # Sample batch information
│   └── README.md                          # Dataset description
│
├── test_data_medium/                       # Medium test dataset (1,000 genes × 12 samples)
│   ├── raw_counts.tsv                     # Raw count data
│   ├── batch_metadata.tsv                 # Sample batch information
│   └── README.md                          # Dataset description
│
├── test_data_large/                        # Large test dataset (5,000 genes × 24 samples)
│   ├── raw_counts.tsv                     # Raw count data
│   ├── batch_metadata.tsv                 # Sample batch information
│   └── README.md                          # Dataset description
│
├── results/                                # Final pipeline outputs
│   ├── normalized_counts.tsv              # DESeq2-style normalized counts
│   ├── normalization_stats.txt            # Size factors and normalization statistics
│   ├── corrected_counts.tsv               # Batch-corrected counts (if batch correction enabled)
│   ├── batch_correction_log.txt           # Batch correction processing log
│   └── summary.txt                        # Comprehensive pipeline summary report
│
└── work/                                   # Nextflow working directory (temporary files)
    ├── 07/
    │   └── 4315c6682960910f75c0ee92c4e3ee/
    │       └── summary.txt                # Process-specific summary
    ├── 34/
    │   └── d047061e8076a874b3b2ea2adc430d/
    │       ├── normalized_counts.tsv      # Intermediate normalized counts
    │       ├── normalization_stats.txt    # Intermediate statistics
    │       └── rnaseq-normalizer          # Process-local binary copy
    ├── 54/
    │   └── 8616295b7ea49b858dcf0b8acd7407/
    │       └── validation_report.txt      # Input validation report
    └── 81/
        └── bf1320833832674e387716c46199e5/
            └── corrected_counts.tsv       # Intermediate batch-corrected counts
```

#### Cargo.toml

```toml
[package]
name = "differential-expression-analyzer"
version = "1.0.0"
edition = "2021"
authors = ["Bioinformatics Pipeline <pipeline@example.com>"]
description = "A robust differential expression analysis tool for RNA-seq data"
license = "MIT"
repository = "https://github.com/username/differential-expression-analyzer"
keywords = ["bioinformatics", "rnaseq", "differential-expression", "genomics", "statistics"]
categories = ["science", "command-line-utilities"]

[[bin]]
name = "diff-expr-analyzer"
path = "src/main.rs"

[dependencies]
clap = "4.4"
serde = { version = "1.0", features = ["derive"] }
serde_json = "1.0"
ndarray = "0.15"
statrs = "0.16"
log = "0.4"
env_logger = "0.10"
csv = "1.3"
anyhow = "1.0"
thiserror = "1.0"

[dev-dependencies]
tempfile = "3.8"
assert_cmd = "2.0"
predicates = "3.0"
approx = "0.5"

[profile.release]
opt-level = 3
lto = true
codegen-units = 1
panic = "abort"

[profile.dev]
opt-level = 0
debug = true
```

#### How to run:

Run main.rs in wsl:

```wsl
(base) trian@triantoharyo:/mnt/c/Users/trian/BGVR/chapter_09/experiment_9_4$ python3 generate_the_data.py --create-test-sets

============================================================
Creating small differential expression dataset
============================================================
Generating normalized counts for 100 genes and 6 samples
Control samples: ['Control_Rep1', 'Control_Rep2', 'Control_Rep3']
Treatment samples: ['Treatment_Rep1', 'Treatment_Rep2', 'Treatment_Rep3']
Normalized count data saved to de_test_data_small/normalized_counts.tsv
Total entries: 600

Sample summary statistics:
                count    mean     std
sample_id
Control_Rep1      100  112.17  148.09
Control_Rep2      100  113.63  160.47
Control_Rep3      100  117.08  171.10
Treatment_Rep1    100  165.20  340.56
Treatment_Rep2    100  167.97  274.89
Treatment_Rep3    100  140.25  233.00

Gene categories:
- Upregulated: 10 genes (10.0%)
- Downregulated: 10 genes (10.0%)
- Unchanged: 80 genes (80.0%)

Sample metadata saved to de_test_data_small/sample_metadata.tsv
        sample_id      group  condition replicate   batch
0    Control_Rep1    Control    Control         1  Batch1
1    Control_Rep2    Control    Control         2  Batch2
2    Control_Rep3    Control    Control         3  Batch3
3  Treatment_Rep1  Treatment  Treatment         1  Batch1
4  Treatment_Rep2  Treatment  Treatment         2  Batch2
5  Treatment_Rep3  Treatment  Treatment         3  Batch3
Dataset created in: de_test_data_small
Files: normalized_counts.tsv, sample_metadata.tsv, README.md, expected_results.txt

============================================================
Creating medium differential expression dataset
============================================================
Generating normalized counts for 1000 genes and 12 samples
Control samples: ['Control_Rep1', 'Control_Rep2', 'Control_Rep3', 'Control_Rep4', 'Control_Rep5', 'Control_Rep6']
Treatment samples: ['Treatment_Rep1', 'Treatment_Rep2', 'Treatment_Rep3', 'Treatment_Rep4', 'Treatment_Rep5', 'Treatment_Rep6']
Normalized count data saved to de_test_data_medium/normalized_counts.tsv
Total entries: 12000

Sample summary statistics:
                count    mean     std
sample_id
Control_Rep1     1000  179.78  381.15
Control_Rep2     1000  177.78  391.13
Control_Rep3     1000  181.06  364.63
Control_Rep4     1000  179.62  386.98
Control_Rep5     1000  176.90  373.53
Control_Rep6     1000  183.00  391.46
Treatment_Rep1   1000  226.51  554.02
Treatment_Rep2   1000  240.46  778.09
Treatment_Rep3   1000  230.01  572.71
Treatment_Rep4   1000  242.19  624.56
Treatment_Rep5   1000  231.34  584.88
Treatment_Rep6   1000  233.34  593.73

Gene categories:
- Upregulated: 100 genes (10.0%)
- Downregulated: 100 genes (10.0%)
- Unchanged: 800 genes (80.0%)

Sample metadata saved to de_test_data_medium/sample_metadata.tsv
         sample_id      group  condition replicate   batch
0     Control_Rep1    Control    Control         1  Batch1
1     Control_Rep2    Control    Control         2  Batch2
2     Control_Rep3    Control    Control         3  Batch3
3     Control_Rep4    Control    Control         4  Batch1
4     Control_Rep5    Control    Control         5  Batch2
5     Control_Rep6    Control    Control         6  Batch3
6   Treatment_Rep1  Treatment  Treatment         1  Batch1
7   Treatment_Rep2  Treatment  Treatment         2  Batch2
8   Treatment_Rep3  Treatment  Treatment         3  Batch3
9   Treatment_Rep4  Treatment  Treatment         4  Batch1
10  Treatment_Rep5  Treatment  Treatment         5  Batch2
11  Treatment_Rep6  Treatment  Treatment         6  Batch3
Dataset created in: de_test_data_medium
Files: normalized_counts.tsv, sample_metadata.tsv, README.md, expected_results.txt

============================================================
Creating large differential expression dataset
============================================================
Generating normalized counts for 5000 genes and 24 samples
Control samples: ['Control_Rep1', 'Control_Rep2', 'Control_Rep3', 'Control_Rep4', 'Control_Rep5', 'Control_Rep6', 'Control_Rep7', 'Control_Rep8', 'Control_Rep9', 'Control_Rep10', 'Control_Rep11', 'Control_Rep12']
Treatment samples: ['Treatment_Rep1', 'Treatment_Rep2', 'Treatment_Rep3', 'Treatment_Rep4', 'Treatment_Rep5', 'Treatment_Rep6', 'Treatment_Rep7', 'Treatment_Rep8', 'Treatment_Rep9', 'Treatment_Rep10', 'Treatment_Rep11', 'Treatment_Rep12']
Normalized count data saved to de_test_data_large/normalized_counts.tsv
Total entries: 120000

Sample summary statistics:
                 count    mean      std
sample_id
Control_Rep1      5000  175.96   452.21
Control_Rep10     5000  177.29   468.23
Control_Rep11     5000  177.88   449.15
Control_Rep12     5000  176.90   459.43
Control_Rep2      5000  176.98   437.14
Control_Rep3      5000  172.79   422.01
Control_Rep4      5000  176.05   448.80
Control_Rep5      5000  176.47   459.86
Control_Rep6      5000  181.18   518.27
Control_Rep7      5000  174.65   446.55
Control_Rep8      5000  172.88   426.98
Control_Rep9      5000  173.57   440.68
Treatment_Rep1    5000  243.27   865.59
Treatment_Rep10   5000  248.40   915.80
Treatment_Rep11   5000  241.65   977.39
Treatment_Rep12   5000  232.99   781.73
Treatment_Rep2    5000  241.29  1056.03
Treatment_Rep3    5000  242.66   921.16
Treatment_Rep4    5000  249.78  1036.59
Treatment_Rep5    5000  244.13  1038.77
Treatment_Rep6    5000  248.49  1083.69
Treatment_Rep7    5000  254.31  1137.29
Treatment_Rep8    5000  252.52  1169.84
Treatment_Rep9    5000  250.26   995.91

Gene categories:
- Upregulated: 500 genes (10.0%)
- Downregulated: 500 genes (10.0%)
- Unchanged: 4000 genes (80.0%)

Sample metadata saved to de_test_data_large/sample_metadata.tsv
          sample_id      group  condition replicate   batch
0      Control_Rep1    Control    Control         1  Batch1
1      Control_Rep2    Control    Control         2  Batch2
2      Control_Rep3    Control    Control         3  Batch3
3      Control_Rep4    Control    Control         4  Batch1
4      Control_Rep5    Control    Control         5  Batch2
5      Control_Rep6    Control    Control         6  Batch3
6      Control_Rep7    Control    Control         7  Batch1
7      Control_Rep8    Control    Control         8  Batch2
8      Control_Rep9    Control    Control         9  Batch3
9     Control_Rep10    Control    Control        10  Batch1
10    Control_Rep11    Control    Control        11  Batch2
11    Control_Rep12    Control    Control        12  Batch3
12   Treatment_Rep1  Treatment  Treatment         1  Batch1
13   Treatment_Rep2  Treatment  Treatment         2  Batch2
14   Treatment_Rep3  Treatment  Treatment         3  Batch3
15   Treatment_Rep4  Treatment  Treatment         4  Batch1
16   Treatment_Rep5  Treatment  Treatment         5  Batch2
17   Treatment_Rep6  Treatment  Treatment         6  Batch3
18   Treatment_Rep7  Treatment  Treatment         7  Batch1
19   Treatment_Rep8  Treatment  Treatment         8  Batch2
20   Treatment_Rep9  Treatment  Treatment         9  Batch3
21  Treatment_Rep10  Treatment  Treatment        10  Batch1
22  Treatment_Rep11  Treatment  Treatment        11  Batch2
23  Treatment_Rep12  Treatment  Treatment        12  Batch3
Dataset created in: de_test_data_large
Files: normalized_counts.tsv, sample_metadata.tsv, README.md, expected_results.txt

# Build the Rust application
cargo build --release

(base) trian@triantoharyo:/mnt/c/Users/trian/BGVR/chapter_09/experiment_9_4$ cargo build --release
   Compiling autocfg v1.4.0
   Compiling proc-macro2 v1.0.95
   Compiling unicode-ident v1.0.18
   Compiling libm v0.2.15
   Compiling libc v0.2.172
   Compiling zerocopy v0.8.25
   Compiling memchr v2.7.4
   Compiling cfg-if v1.0.0
   Compiling paste v1.0.15
   Compiling syn v1.0.109
   Compiling num-traits v0.2.19
   Compiling matrixmultiply v0.3.10
   Compiling bytemuck v1.23.0
   Compiling safe_arch v0.7.4
   Compiling utf8parse v0.2.2
   Compiling serde v1.0.219
   Compiling quote v1.0.40
   Compiling rawpointer v0.2.1
   Compiling syn v2.0.101
   Compiling getrandom v0.2.16
   Compiling typenum v1.18.0
   Compiling rand_core v0.6.4
   Compiling wide v0.7.32
   Compiling num-integer v0.1.46
   Compiling num-complex v0.4.6
   Compiling ppv-lite86 v0.2.21
   Compiling approx v0.5.1
   Compiling anstyle-parse v0.2.6
   Compiling aho-corasick v1.1.3
   Compiling is_terminal_polyfill v1.70.1
   Compiling rand_chacha v0.3.1
   Compiling anstyle-query v1.1.2
   Compiling regex-syntax v0.8.5
   Compiling serde_derive v1.0.219
   Compiling rand v0.8.5
   Compiling colorchoice v1.0.3
   Compiling anstyle v1.0.10
   Compiling anstream v0.6.18
   Compiling simba v0.6.0
   Compiling rand_distr v0.4.3
   Compiling regex-automata v0.4.9
   Compiling num-rational v0.4.2
   Compiling itoa v1.0.15
   Compiling anyhow v1.0.98
   Compiling nalgebra-macros v0.1.0
   Compiling clap_lex v0.7.4
   Compiling thiserror v1.0.69
   Compiling strsim v0.11.1
   Compiling serde_json v1.0.140
   Compiling ryu v1.0.20
   Compiling clap_builder v4.5.38
   Compiling thiserror-impl v1.0.69
   Compiling regex v1.11.1
   Compiling is-terminal v0.4.16
   Compiling csv-core v0.1.12
   Compiling log v0.4.27
   Compiling termcolor v1.4.1
   Compiling lazy_static v1.5.0
   Compiling humantime v2.2.0
   Compiling env_logger v0.10.2
   Compiling ndarray v0.15.6
   Compiling clap v4.5.38
   Compiling nalgebra v0.29.0
   Compiling csv v1.3.1
   Compiling statrs v0.16.1
   Compiling differential-expression-analyzer v1.0.0 (/mnt/c/Users/trian/BGVR/chapter_09/experiment_9_4)
    Finished `release` profile [optimized] target(s) in 1m 30s
```

Run main.nf in wsl:

```wsl
nextflow run main.nf \
    --input de_test_data_medium/normalized_counts.tsv \
    --metadata de_test_data_medium/sample_metadata.tsv \
    --output_dir results

(base) trian@triantoharyo:/mnt/c/Users/trian/BGVR/chapter_09/experiment_9_4$ nextflow run main.nf \
    --input de_test_data_medium/normalized_counts.tsv \
    --metadata de_test_data_medium/sample_metadata.tsv \
    --output_dir results
Nextflow 25.04.2 is available - Please consider updating your version to it

 N E X T F L O W   ~  version 24.10.4

Launching `main.nf` [happy_fourier] DSL2 - revision: 6d6b63180f


Differential Expression Analysis Pipeline
=========================================
Input file:           de_test_data_medium/normalized_counts.tsv
Metadata file:        de_test_data_medium/sample_metadata.tsv
Output directory:     results
Control group:        Control
Treatment group:      Treatment
Significance level:   0.05
Min count threshold:  10.0

executor >  local (5)
[e9/d42e33] BUILD_ANALYZER (building_diff_expr_analyzer)    [100%] 1 of 1 ✔
[36/1d39b0] VALIDATE_INPUTS (validating_inputs)             [100%] 1 of 1 ✔
[f7/c201d7] DIFFERENTIAL_EXPRESSION (differential_analysis) [100%] 1 of 1 ✔
[2a/d22051] GENERATE_PLOTS_DATA (generating_plot_data)      [100%] 1 of 1 ✔
[29/302176] GENERATE_SUMMARY (generating_summary)           [100%] 1 of 1 ✔

    🎯 Differential Expression Analysis Summary
    ==========================================
    Completed at: 2025-05-25T03:09:22.018902827+07:00
    Duration:     1m 26s
    Success:      true
    Work dir:     /mnt/c/Users/trian/BGVR/chapter_09/experiment_9_4/work
    Exit status:  0


        🎉 SUCCESS! Differential expression analysis completed successfully.

        📁 Results available in: results/

        📊 Key outputs:
        ├── differential_expression.tsv    📈 Main DE results
        ├── de_analysis_stats.txt          📋 Detailed statistics
        ├── volcano_plot_data.tsv          🌋 Volcano plot data
        ├── ma_plot_data.tsv               📊 MA plot data
        ├── input_validation_report.txt    ✅ Validation details
        └── summary.txt                    📄 Comprehensive summary

        💡 Next steps:
        1. Review summary.txt for overview
        2. Check de_analysis_stats.txt for detailed statistics
        3. Use plot data files for visualization
        4. Filter significant genes for pathway analysis

        🔬 Your differential expression analysis is ready for biological interpretation!

Completed at: 25-May-2025 03:09:22
Duration    : 1m 25s
CPU hours   : (a few seconds)
Succeeded   : 5
```

#### Output

**differential_expression.tsv:**

```
gene_id	control_mean	treatment_mean	log2_fold_change	p_value	adjusted_p_value	significant
GENE_00001	138.467454	417.670217	1.585886	1.31e-2	7.48e-2	false
GENE_00002	81.288552	335.530279	2.031973	3.38e-3	2.83e-2	true
GENE_00003	165.433666	740.769551	2.156024	1.41e-3	2.00e-2	true
GENE_00004	32.640307	149.826006	2.164622	1.34e-4	5.76e-3	true
GENE_00005	108.540265	530.577967	2.278820	1.75e-3	2.16e-2	true
GENE_00006	90.931687	342.379864	1.901171	1.43e-2	7.97e-2	false
GENE_00007	875.935766	4991.133094	2.509113	7.21e-3	4.69e-2	true
GENE_00008	7.332544	31.236137	1.951850	8.61e-3	5.41e-2	false
GENE_00009	184.906541	1066.630124	2.521762	1.55e-3	2.08e-2	true
GENE_00010	21.294963	107.200333	2.278915	4.16e-3	3.22e-2	true
GENE_00011	69.005172	375.868384	2.428527	1.82e-3	2.16e-2	true
GENE_00012	5.719729	30.779486	2.241621	7.72e-5	5.16e-3	true
GENE_00013	185.890597	563.264005	1.594176	1.83e-2	9.40e-2	false
GENE_00014	263.146946	1545.420130	2.549520	1.05e-2	6.38e-2	false
GENE_00015	122.007745	576.636616	2.231413	6.04e-3	4.18e-2	true
GENE_00016	129.945227	1173.108818	3.164531	1.42e-2	7.97e-2	false
GENE_00017	74.739039	369.367418	2.289848	1.52e-2	8.28e-2	false
GENE_00018	94.670368	479.851140	2.329446	4.20e-3	3.22e-2	true
GENE_00019	151.407812	932.022395	2.613975	2.87e-3	2.55e-2	true
...
```

**corrected_counts.tsv:**

```
gene_id	sample_id	corrected_count
GENE_00001	Control_Batch1_Rep1	416.107304
GENE_00001	Control_Batch1_Rep4	399.222937
GENE_00001	Control_Batch2_Rep3	364.508815
GENE_00001	Control_Batch2_Rep6	368.810845
GENE_00001	Control_Batch3_Rep2	426.077123
GENE_00001	Control_Batch3_Rep5	402.664813
GENE_00001	Treatment_Batch1_Rep2	221.811384
GENE_00001	Treatment_Batch1_Rep5	203.863980
GENE_00001	Treatment_Batch2_Rep1	182.375362
GENE_00001	Treatment_Batch2_Rep4	197.560530
GENE_00001	Treatment_Batch3_Rep3	205.459044
GENE_00001	Treatment_Batch3_Rep6	220.200622
GENE_00002	Control_Batch1_Rep1	116.309511
GENE_00002	Control_Batch1_Rep4	114.638324
GENE_00002	Control_Batch2_Rep3	100.619620
GENE_00002	Control_Batch2_Rep6	105.283373
GENE_00002	Control_Batch3_Rep2	125.580626
GENE_00002	Control_Batch3_Rep5	95.446474
GENE_00002	Treatment_Batch1_Rep2	210.515433
GENE_00002	Treatment_Batch1_Rep5	193.619558
GENE_00002	Treatment_Batch2_Rep1	212.988369
GENE_00002	Treatment_Batch2_Rep4	221.448568
GENE_00002	Treatment_Batch3_Rep3	231.720726
GENE_00002	Treatment_Batch3_Rep6	263.316828
...
```

**ma_plot_data.tsv:**

```
gene_id	average_expression	log2_fold_change	significant	regulation
GENE_00001	8.124477	1.585886	False	none
GENE_00002	7.710183	2.031973	True	up
GENE_00003	8.826871	2.156024	True	up
GENE_00004	6.527214	2.164622	True	up
GENE_00005	8.324447	2.278820	True	up
GENE_00006	7.765904	1.901171	False	none
GENE_00007	11.518916	2.509113	True	up
GENE_00008	4.342294	1.951850	False	none
GENE_00009	9.291788	2.521762	True	up
GENE_00010	6.027854	2.278915	True	up
GENE_00011	7.803723	2.428527	True	up
GENE_00012	4.266757	2.241621	True	up
GENE_00013	8.552966	1.594176	False	none
GENE_00014	9.822226	2.549520	False	none
GENE_00015	8.452539	2.231413	True	up
GENE_00016	9.349894	3.164531	False	none
GENE_00017	7.801244	2.289848	False	none
GENE_00018	8.171231	2.329446	True	up
GENE_00019	9.084051	2.613975	True	up
...
```

**summary.txt:**

```
Differential Expression Analysis Pipeline Summary
================================================
Generated: 2025-05-25 03:09:21

PIPELINE STATUS: SUCCESS
========================
✅ Input Validation: COMPLETED
✅ Differential Analysis: COMPLETED
✅ Plot Data Generation: COMPLETED
✅ Statistical Reporting: COMPLETED

INPUT PARAMETERS:
================
- Normalized counts file: de_test_data_medium/normalized_counts.tsv
- Sample metadata file: de_test_data_medium/sample_metadata.tsv
- Control group: Control
- Treatment group: Treatment
- Significance threshold: 0.05
- Minimum count threshold: 10.0
- Output directory: results

DATA PROCESSING RESULTS:
=======================
- Input count entries: 12000
- Metadata samples: 12
- Genes analyzed: 905
- Significant genes: 141
- Significance rate: 15.5%

OUTPUT FILES:
============
✓ differential_expression.tsv (905 genes analyzed)
✓ de_analysis_stats.txt (detailed statistics)
✓ volcano_plot_data.tsv (volcano plot data)
✓ ma_plot_data.tsv (MA plot data)
✓ input_validation_report.txt (validation details)
✓ summary.txt (this summary)

NEXT STEPS:
===========
1. Review de_analysis_stats.txt for detailed statistics
2. Use differential_expression.tsv for downstream analysis
3. Create volcano plot using volcano_plot_data.tsv
4. Create MA plot using ma_plot_data.tsv
5. Filter significant genes for pathway analysis

STATISTICAL SUMMARY:
===================
From detailed analysis:
- Total genes analyzed: 905
- Significant genes: 141 (15.6%)
- Upregulated genes: 61 (6.7%)
- Downregulated genes: 80 (8.8%)

VALIDATION REPORT:
=================
Input Validation Report
=======================
Count data: 12000 entries
Count header: ['gene_id', 'sample_id', 'count']
Metadata samples: 12
Groups found: ['Treatment', 'Control']
✓ Required groups found
Common samples: 12
✓ Sufficient samples for analysis
PIPELINE COMPLETED SUCCESSFULLY
Data is ready for downstream analysis and visualization.

For visualization in R:
- Use volcano_plot_data.tsv for ggplot2 volcano plots
- Use ma_plot_data.tsv for MA plots
- Filter significant genes: awk -F'\t' '=="true"' differential_expression.tsv

For pathway analysis:
- Extract significant gene lists from differential_expression.tsv
- Use tools like GSEA, DAVID, or Enrichr for functional annotation
```

**the_analysis_stats.txt:**

```
Differential Expression Analysis Statistics
==========================================
Total genes analyzed: 905
Total samples: 12
Control samples: 6
Treatment samples: 6

Significant genes: 141 (15.6%)
Upregulated genes: 61 (6.7%)
Downregulated genes: 80 (8.8%)

Top 10 Most Significant Genes:
gene_id	log2FC	adj_p_value
GENE_00105	-1.754	1.52e-3
GENE_00156	-1.657	1.52e-3
GENE_00116	-1.769	1.55e-3
GENE_00112	-2.134	1.66e-3
GENE_00130	-1.779	1.78e-3
GENE_00176	-2.067	1.78e-3
GENE_00108	-1.342	3.07e-3
GENE_00157	-1.678	3.07e-3
GENE_00081	2.563	3.72e-3
GENE_00115	-1.850	3.72e-3
```
**volcano_plot_data.tsv:**

```
gene_id	log2_fold_change	neg_log10_p_value	significant	regulation
GENE_00001	1.585886	1.126098	False	none
GENE_00002	2.031973	1.548214	True	up
GENE_00003	2.156024	1.698970	True	up
GENE_00004	2.164622	2.239578	True	up
GENE_00005	2.278820	1.665546	True	up
GENE_00006	1.901171	1.098542	False	none
GENE_00007	2.509113	1.328827	True	up
GENE_00008	1.951850	1.266803	False	none
GENE_00009	2.521762	1.681937	True	up
GENE_00010	2.278915	1.492144	True	up
GENE_00011	2.428527	1.665546	True	up
GENE_00012	2.241621	2.287350	True	up
GENE_00013	1.594176	1.026872	False	none
GENE_00014	2.549520	1.195179	False	none
GENE_00015	2.231413	1.378824	True	up
GENE_00016	3.164531	1.098542	False	none
GENE_00017	2.289848	1.081970	False	none
GENE_00018	2.329446	1.492144	True	up
GENE_00019	2.613975	1.593460	True	up
...
```

**expected_results.txt:**

```
Expected Differential Expression Results
========================================

Dataset: Large
Genes: 5,000
Samples per group: 12

Expected significant genes: ~1000 (20%)
Expected upregulated: ~500 (10%)
Expected downregulated: ~500 (10%)

Statistical Power:
- With 12 samples per group
- Expected to detect 2+ fold changes
- At 5% FDR with 80%+ power

Visualization Expectations:
- Volcano plot: Clear separation of significant genes
- MA plot: Even distribution around zero for non-DE genes
- P-value histogram: Enrichment of small p-values
```

**normalized_counts.tsv:**

```
gene_id	sample_id	count
GENE_00001	Control_Rep1	110.34280612766383
GENE_00001	Control_Rep2	139.68331029539823
GENE_00001	Control_Rep3	181.63155034179445
GENE_00001	Control_Rep4	107.21382631080921
GENE_00001	Control_Rep5	107.21435437582568
GENE_00001	Control_Rep6	184.71887489659503
GENE_00001	Control_Rep7	144.79249496594912
GENE_00001	Control_Rep8	99.90591706060374
GENE_00001	Control_Rep9	135.34664695268415
GENE_00001	Control_Rep10	100.08761192534743
GENE_00001	Control_Rep11	100.01821340420719
GENE_00001	Control_Rep12	123.67534058569831
GENE_00001	Treatment_Rep1	297.4752335065515
GENE_00001	Treatment_Rep2	318.0718961999217
GENE_00001	Treatment_Rep3	474.3121968266527
GENE_00001	Treatment_Rep4	749.5063791247022
GENE_00001	Treatment_Rep5	343.73716275961493
GENE_00001	Treatment_Rep6	496.77211763726615
GENE_00001	Treatment_Rep7	671.2662643146607
...
```

**sample_metadata.tsv:**

```
sample_id	group	condition	replicate	batch
Control_Rep1	Control	Control	1	Batch1
Control_Rep2	Control	Control	2	Batch2
Control_Rep3	Control	Control	3	Batch3
Control_Rep4	Control	Control	4	Batch1
Control_Rep5	Control	Control	5	Batch2
Control_Rep6	Control	Control	6	Batch3
Treatment_Rep1	Treatment	Treatment	1	Batch1
Treatment_Rep2	Treatment	Treatment	2	Batch2
Treatment_Rep3	Treatment	Treatment	3	Batch3
Treatment_Rep4	Treatment	Treatment	4	Batch1
Treatment_Rep5	Treatment	Treatment	5	Batch2
Treatment_Rep6	Treatment	Treatment	6	Batch3
```



