## 4.7. eQTL Mapping and Functional Variant Discovery

### experiment_47

An AI engineer tasked with implementing an eQTL pipeline in Rust might structure it as follows: first parse genotype data into in-memory arrays or memory-mapped structures (using memmap2) for large data sets; next, run association tests in parallel (split by gene or SNP) using rayon or MPI-based approaches; then, merge partial results in a final step. The Rust code below demonstrates a parallel solution for computing eQTL (expression quantitative trait loci) associations between SNP data and gene expression data. Each SNP is processed on a separate thread, while the linear regression-based eQTL analysis is carried out locally for all genes. The results are then merged into a single collection, providing a scalable approach for large genomic datasets in HPC or cloud environments.

After defining structs to represent SNP genotypes and gene expressions, the program uses Rayon’s par_iter to distribute SNP processing across multiple CPU cores. Within each thread, all genes are iterated over to perform a basic linear regression (linear_eqtl) that estimates the slope and p-value (using a Student’s t-distribution for significance testing). A parallel reduction (reduce_with) then concatenates the locally computed results from each thread into a final output vector, which is written to a CSV file. This design eliminates the complexity of nested parallel iterators and ensures straightforward scalability across large numbers of SNPs.

#### Files contents:
* experiment_47/
  * Cargo.toml (Cargo.toml file for dependencies)
* experiment_47/src/
  * main.rs (rust script)
  * output.txt (output file)
  * partial_eqtl.csv

#### How to run:

run main.rs in powershell:

```powershell
cargo run | tee output.txt
```
(run main.rs and get the partial_eqtl.csv output and output.txt)

#### [dependencies]

```toml
rayon = "1.10.0"
statrs = "0.18.0"
```

#### Explanation of the Output
This Rust program performs parallel eQTL (expression Quantitative Trait Loci) analysis using linear regression for each SNP-gene pair. The results are stored in partial_eqtl.csv.

##### Output Breakdown
The program analyzes 2 SNPs (rs1, rs2) and 2 genes (GeneA, GeneB), resulting in 4 eQTL tests.

Each eQTL test provides:

1. SNP ID – The SNP being tested.
2. Gene ID – The gene whose expression is modeled.
3. Slope – The estimated effect size of the SNP on gene expression.
4. p-value – Statistical significance of the association.

The results:

```csv
rs1    GeneA    2.000    0.04052
rs1    GeneB    0.150    0.43486
rs2    GeneA   -0.750    0.41953
rs2    GeneB   -0.125    0.34433
```

##### Interpretation

* rs1 → GeneA:

  * Slope = 2.000 → Strong positive association.
  * p-value = 0.04052 → Statistically significant (p < 0.05).

* rs1 → GeneB:

  * Slope = 0.150 → Weak association.
  * p-value = 0.43486 → Not significant (p > 0.05).

* rs2 → GeneA:
  * Slope = -0.750 → Negative association.
  * p-value = 0.41953 → Not significant.

* rs2 → GeneB:
  * Slope = -0.125 → Almost no effect.
  * p-value = 0.34433 → Not significant.

#### Conclusion
* The program correctly runs parallel eQTL analysis using linear regression and writes results to a CSV file.
* rs1 is significantly associated with GeneA, suggesting a genetic influence on its expression.
* Other SNP-gene pairs show no significant association, meaning the SNPs likely do not affect those genes.
