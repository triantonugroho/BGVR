## 1.5. Acquiring and Cleaning Data

### experiment_1_5_4

For more advanced data transformations, developers might adopt polars for tabular data. Suppose a user wants to parse metadata from a CSV file that describes sample conditions (e.g., diseased vs. healthy) and merge it with alignment statistics. polars can load large CSV files quickly while Rust’s strong type checking ensures columns are properly named and typed. The next code block outlines a minimal example:

This code demonstrates two functions (read_csv_basic and read_csv_advanced) for reading CSV files using the Polars library in Rust. The first function applies relatively simple settings (e.g., inferring the schema, assuming a header row, and attempting to parse date columns), while the second function showcases more advanced parsing options (such as customizing the separator, quote character, encoding, and handling specific null value markers). Both functions ultimately read the CSV file into a DataFrame and return it, and the main function prints the resulting data frames. It’s crucial to always consult the latest Polars documentation (or any crate’s documentation) before writing and generating code for GenAI, as APIs and best practices can change over time, ensuring that your code remains accurate and up-to-date.

#### Project Structure:

```plaintext
experiment_1_5_4/
├── Cargo.toml                     # Rust project configuration and dependencies
└── src/
    ├── main.rs                    # Main Rust script containing program logic
    ├── output.csv                 # CSV output file
    ├── Synthesize CSV File.ipynb  # Python notebook to synthesize output.csv
    └── output.txt                 # Text output file
```

#### Cargo.toml

```toml
[package]
name = "experiment_1_5_4"
version = "0.1.0"
edition = "2021"

[dependencies]
polars = { version = "0.46.0", features = ["csv", "lazy"] }
```

#### How to run:

run in powershell:

```powershell
cargo run | tee output.txt
```

(run main.rs and save the output in output.txt)
  

#### Explanation of the Output

##### Overview of the Rust Program
This Rust program reads a CSV file (output.csv), processes it using the Polars library, and performs data summarization. The program consists of three main steps:

1. Reading the CSV file using both basic and advanced configurations.
2. Displaying the first 5 rows of the CSV file for both configurations.
3. Summarizing data by grouping it based on disease_state and computing:
   * Average quality score (avg_quality)
   * Average age (avg_age)
   * Sample count (sample_count)

##### Step-by-Step Execution and Analysis

###### 1. Reading the CSV File (Basic & Advanced)

* The program attempts to read output.csv using two methods:
  * Basic CSV Reader: Reads the file with default options.
  * Advanced CSV Reader: Uses additional configurations like:
    * Comma separator (with_separator(b','))
    * Handling missing values ("NA" and "null")
* Output (First 5 Rows of output.csv)
  * The CSV file contains 14 columns, including sample_id, disease_state, treatment, patient_sex, quality_flag, and age_group.
  * The column headers are truncated due to formatting.
  * Each row represents a sample with attributes such as collection date, treatment type, and age group.

* Key Observations: 
  * Both basic and advanced readers return the same first 5 rows of the dataset.
  * The column names are truncated in the output due to formatting issues.
  * The notes column contains long text descriptions, which cause misalignment.

###### 2. Data Summarization: Grouping by Disease State

After reading the CSV file, the program groups the data by disease_state and calculates:
* avg_quality (mean of quality_score)
* avg_age (mean of patient_age)
* sample_count (number of samples in each group)

Output Summary Table:
┌───────────────┬─────────────┬───────────┬──────────────┐
│ disease_state ┆ avg_quality ┆ avg_age   ┆ sample_count │
│ ---           ┆ ---         ┆ ---       ┆ ---          │
│ str           ┆ f64         ┆ f64       ┆ u32          │
╞═══════════════╪═════════════╪═══════════╪══════════════╡
│ Healthy       ┆ 7.642857    ┆ 51.714286 ┆ 21           │
│ Diseased      ┆ 6.2         ┆ 56.55     ┆ 20           │
│ Remission     ┆ 7.788889    ┆ 33.444444 ┆ 9            │
└───────────────┴─────────────┴───────────┴──────────────┘

* Key Observations:
  * The dataset contains 3 disease states: Healthy, Diseased, and Remission.
  * Healthy individuals have the highest sample count (21).
  * Remission patients have the lowest average age (33.4 years).
  * Diseased patients have the lowest average quality score (6.2).

#### Conclusion

* Functionality: The program successfully reads, processes, and summarizes a CSV dataset using Polars in Rust.
* Data Integrity: The basic and advanced CSV readers produce consistent results.
* Summarization Insights:
  * Healthy individuals dominate the dataset.
  * Diseased individuals tend to be older (average age: 56.55).
  * Remission patients are younger and have higher quality scores.
* Future Improvements
  *  Fix column name truncation in printed output.
  *  Handle long text fields (notes) properly.
  *  Extend analysis with more statistical insights (e.g., standard deviations, distributions).

Overall, the program successfully processes biomedical data for analysis. 


