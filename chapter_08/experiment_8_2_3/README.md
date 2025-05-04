# 8.2. Data Structures for Variant Representation

## 8.1.3 Practical Perspectives and HPC Concurrency

### experiment_8_1_3

Below is an illustrative Rust code snippet that computes genotype frequencies, performs a chi-square–based HW p-value calculation, and exploits concurrency to handle separate chunks of the genome. This example is adapted for industrial-scale usage by incorporating crates that enhance numerical stability, data handling, and concurrency.

In this code, rust-htslib is employed for reading and parsing VCF/BCF files. The rayon crate can be used for parallel iterators if multiple shards of the genome need processing concurrently. ndarray assists with numerical operations for genotype frequency calculations, while statrs ensures robust statistical distributions. polars offers efficient DataFrame manipulation to keep results organized in memory. With thoughtful adjustments to memory usage, chunk sizes, and concurrency levels, this approach scales up to industrial-size datasets typically encountered in large consortia or AI-driven pharmaceutical pipelines.

#### Files contents:
* experiment_8_1_3/
  * Cargo.toml (Cargo.toml file for dependencies)
* experiment_8_1_3/src/
  * main.rs (rust script)
  * synthetic.vcf (synthetic vcf file for input file)
  * synthetic.vcf.hw_results.csv (synthetic.vcf result csv file)
  * output.txt (text file output)
* experiment_8_1_3/src/work/65/54b3be71f96f81cbb7987c87cd42f1/
  * test.coverage.json (test coverage json file)
* experiment_8_1_3/src/work/f0/c029dbc59728cf12ca4ea10d38edb5/
  * test.coverage.json (test coverage json file)
* experiment_8_1_3/src/work/fb/a19c2d0203adcbff8bc1a8c54dc6c6/
  * merged_coverage.json (merged coverage json file)

#### How to run:

run main.rs in wsl:

```wsl
cargo run -- synthetic.vcf 0 1000000
```

(run main.rs and create synthetic.vcf.hw_results.csv output)

#### [dependencies]

```toml
rust-htslib = "0.49.0"
rayon = "1.5.1"
ndarray = "0.16.1"
statrs = "0.18.0"
polars = { version = "0.46", features = ["lazy"] }
```

#### Explanation of the Output

##### ✨ What Happens in the Code (main.rs)?

Overall:
This Rust program reads a VCF file, parses each variant, calculates a Hardy-Weinberg equilibrium p-value for each position using a chi-square test, then saves the results to a CSV file.

##### 🔍 Step-by-Step Explanation:
###### 1. Libraries:
* statrs: To do chi-square test and compute p-value.
* polars: To create a DataFrame (like a table) and easily save it as a CSV.
* rust-htslib is mentioned but not actually used (this is manual parsing).

###### 2. Key Functions:
###### (A) chi_square_hw(aa, ab, bb, p)
* Inputs:
  * aa = count of homozygous reference (0/0) samples
  * ab = count of heterozygous (0/1 or 1/0) samples
  * bb = count of homozygous alternate (1/1) samples
  * p = allele frequency (reference allele)
* Calculates expected counts under Hardy-Weinberg equilibrium.
* Calculates a chi-square statistic.
* Converts chi-square statistic to a p-value.
  * A high p-value (near 1) = good fit to HW equilibrium.
  * A low p-value (near 0) = deviation from HW.

###### (B) process_vcf_file(vcf_path, start_pos, end_pos)
* Opens the VCF file line-by-line.
* Skips meta-information lines (##...).
* Parses the #CHROM header to detect sample columns.
* For each variant line:
  * Filters variants based on position (start_pos to end_pos).
  * Extracts genotypes from samples.
  * Counts 0/0, 0/1, 1/1 occurrences.
  * Calculates allele frequency p.
  * Calculates Hardy-Weinberg p-value.
  * Stores info in a list.
* Creates a DataFrame from the list.
* Returns the DataFrame.

###### (C) main()
* Reads command line arguments:
  * VCF path
  * Start and End positions (optional)
* Calls process_vcf_file.
* Prints the DataFrame to terminal.
* Saves the DataFrame to a .csv file (same name as VCF + .hw_results.csv).

##### 📂 Your Input Dataset: synthetic.vcf

```text
##fileformat=VCFv4.2
##contig=<ID=1,length=249250621>
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Sample1 Sample2 Sample3
1       12345   .       A       G       50.0    PASS    NS=3    GT      0/0     0/1     1/1
1       67890   .       T       C       40.0    PASS    NS=3    GT      0/1     1/1     0/0
```

Two variants (positions 12345 and 67890) across 3 samples (Sample1, Sample2, Sample3).

##### 🧪 Calculation Details per Variant

###### Variant 1:
* Position: 12345
* Genotypes: 0/0, 0/1, 1/1
* Counts:
  * 0/0 → 1 sample
  * 0/1 → 1 sample
  * 1/1 → 1 sample

* Allele Frequency (p):

### Allele Frequency (p)

$$
p = \frac{(2 \times \text{count}(AA) + \text{count}(AB))}{2 \times \text{total samples}}
$$

Substituting values:

$$
p = \frac{(2 \times 1) + 1}{2 \times 3} = \frac{3}{6} = 0.5
$$

---

### Expected Genotype Counts

- Homozygous Reference (AA, `0/0`):

$$
\text{Expected}(AA) = p^2 \times 3 = (0.5)^2 \times 3 = 0.75
$$

- Heterozygous (AB, `0/1` or `1/0`):

$$
\text{Expected}(AB) = 2 \times p \times (1 - p) \times 3 = 2 \times 0.5 \times 0.5 \times 3 = 1.5
$$

- Homozygous Alternate (BB, `1/1`):

$$
\text{Expected}(BB) = (1 - p)^2 \times 3 = (0.5)^2 \times 3 = 0.75
$$

---

### Chi-Square Statistic

Chi-Square formula:

$$
\chi^2 = \sum \frac{(\text{Observed} - \text{Expected})^2}{\text{Expected}}
$$

Substituting values:

$$
\chi^2 = \frac{(1 - 0.75)^2}{0.75} + \frac{(1 - 1.5)^2}{1.5} + \frac{(1 - 0.75)^2}{0.75}
$$

$$
\chi^2 = \frac{0.0625}{0.75} + \frac{0.25}{1.5} + \frac{0.0625}{0.75}
$$

$$
\chi^2 = 0.0833 + 0.1667 + 0.0833 = 0.3333
$$

---

### P-Value Calculation

Using the Chi-Square distribution with 1 degree of freedom:

$$
p\text{-value} = 1 - \text{CDF}(\chi^2, df=1)
$$

Substituting the value:

$$
p\text{-value} = 1 - \text{CDF}(0.3333, 1) \approx 0.5637
$$

---

### Final Interpretation

- If \( p\text{-value} > 0.05 \), there is **no significant deviation** from Hardy-Weinberg Equilibrium (HWE).
- In this case:

$$
p\text{-value} = 0.5637 > 0.05
$$

✅ **Conclusion**: This variant **is in Hardy-Weinberg Equilibrium**.



##### 📝 Output CSV: synthetic.vcf.hw_results.csv

| Chromosome | Position | Reference Allele | Alternate Allele | HWE p-value |
|:----------:|:--------:|:----------------:|:----------------:|:-----------:|
| 1          | 12345    | A                | G                | 0.5637      |
| 1          | 67890    | T                | C                | 0.5637      |

##### 🛠 Executable: vcf_analysis

* Compiled main.rs into a binary.
* Can be called as:

```text
./vcf_analysis synthetic.vcf
```

It will generate and print the CSV automatically.

#### ✅ Conclusion
* Purpose: This program checks if genetic variants are in Hardy-Weinberg equilibrium.

* Result interpretation:

* Both variants have high p-values (~0.56) → No significant deviation from HW equilibrium.

* Why useful?

  * Hardy-Weinberg deviations may indicate genotyping errors, population stratification, selection, etc.

In your case, the dataset is tiny and idealized, so the p-values are expectedly moderate.


