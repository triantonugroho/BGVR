## 1.5. Acquiring and Cleaning Data

### experiment_1_5_3

An AI engineer tasked with data preprocessing in Rust may develop specialized utilities that parse metadata, filter reads by quality, or convert between genomic file formats. The following code snippet shows a streamlined Rust tool that demonstrates how to read a FASTQ file using bio and apply a simple filter. It also includes logging for HPC environments, using crates such as env_logger to track progress and potential errors. In production, these steps might be integrated into a Nextflow pipeline, ensuring that each HPC job references pinned versions of Rust and its crates.

In this snippet, the bio crate handles FASTQ parsing and writing. The env_logger plus log approach yields HPC-friendly logs, with environment variables controlling verbosity. For industrial-scale usage, engineers commonly add robust error handling with matching patterns, ensuring that each read record is either processed or logged if corrupted. They may also incorporate concurrency through Rayon if the dataset is large, spawning multiple threads to parse different file segments concurrently.

#### Project Structure:
```plaintext
experiment_1_5_3/
├── Cargo.toml                     # Rust project configuration and dependencies
└── src/
    ├── main.rs                    # Main Rust script containing program logic
    ├── reads.fastq.rar            # Compressed FASTQ reads file
    ├── filtered_reads.fastq.rar   # Compressed filtered FASTQ reads file
    └── output.txt                 # Output file
```

#### Cargo.toml

```toml
[package]
name = "experiment_1_5_3"
version = "0.1.0"
edition = "2021"

[dependencies]
bio = "2.2.0"        # For FASTQ/FASTA parsing
log = "0.4.20"       # For logging facade
env_logger = "0.11.6" # For environment-based logging configuration
```

#### How to run:

run in powershell:

```powershell
cargo run | tee output.txt
```

(run main.rs and save the output in output.txt)
  

#### Explanation of the Output

##### Overview of the Rust Program
This Rust program reads a FASTQ file (reads.fastq), applies a basic filtering rule, and writes the filtered reads to a new FASTQ file (filtered_reads.fastq).

* FASTQ Files contain biological sequence data, commonly used in genomics.
* The filtering rule applied is "keep only reads with length > 75".
* Logging (log crate) is used to handle errors and information messages.

##### Step-by-Step Execution

######1. Logging Initialization

```rust
env_logger::init();
```

* Enables logging (e.g., info!, error!) based on environment variables.

###### 2. Reading the FASTQ File

```rust
let reader = fastq::Reader::from_file("reads.fastq")?
```

* Attempts to open reads.fastq for reading.
* If the file is missing or inaccessible, logs an error and exits.

###### 3. Creating an Output File

```rust
let mut writer = fastq::Writer::to_file("filtered_reads.fastq")?;
```

* Opens filtered_reads.fastq for writing the filtered reads.

###### 4.Filtering Reads

```rust
for record_result in reader.records() {
    let record = record_result?;
    if record.seq().len() > 75 {
        writer.write_record(&record)?;
    }
}
```

* Iterates over each read in reads.fastq.
* If a read’s length is greater than 75, it is written to filtered_reads.fastq.
* If an error occurs while reading, it is logged.

###### 5. Completion Message

```rust
info!("Filtering complete. Output saved to filtered_reads.fastq");
println!("reads.fastq have successfully filtered with a trivial filter: keep only reads with length > 75 and saved to filtered_reads.fastq");
```

* Logs the completion message.
* Prints the success message to the terminal.

##### Output Analysis

```rust
reads.fastq have successfully filtered with a trivial filter: keep only reads with length > 75 and saved to filtered_reads.fastq
```
* This confirms that:
  * The program successfully opened reads.fastq.
  * Reads were filtered based on the length criterion (>75).
  * The filtered data was saved to filtered_reads.fastq.

####Conclusion

* Functionality: The script successfully filters biological sequence reads from a FASTQ file and saves them to a new file.
* Reproducibility: Can be integrated into HPC workflows (e.g., Nextflow) for large-scale genomics pipelines.
* Logging & Error Handling: Uses Rust’s log crate to handle errors gracefully.
* Scalability: The simple filtering logic can be extended for more complex bioinformatics processing. 

Overall, the program performed as expected without errors.

