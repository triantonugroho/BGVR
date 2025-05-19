## 1.1. Introduction to Rust Programming Language

### experiment_1_1_10

An example demonstrating how Rust’s built-in data structures, such as Vectors, HashMaps, and Strings, can be combined with the bio crate for typical bioinformatics tasks. In this snippet, we parse a FASTQ file, store each record’s length in a HashMap keyed by read ID, and then summarize the collected statistics.

This code reads each record from a FASTQ file via the Reader provided by the bio crate, storing the read ID and its length in a HashMap. Rust’s Strings and HashMaps offer a straightforward and memory-safe way to handle complex genomic data sets, since the compiler ensures that references to these data structures are correctly managed at runtime. The result is a type-safe approach to file parsing and data manipulation in large-scale bioinformatics tasks, aligning with HPC needs while preserving code clarity and reliability.

#### Project Structure:

```plaintext
experiment_1_1_10/
├── Cargo.toml                     # Rust project configuration and dependencies
└── src/
    ├── main.rs                    # Main Rust script containing program logic
    ├── reads.fastq.rar            # Compressed FASTQ reads file
    └── output.txt                 # Output text file
 ```

#### Cargo.toml

```toml
[package]
name = "experiment_1_1_10"
version = "0.1.0"
edition = "2021"

[dependencies]
bio = "2.2.0"
```

#### How to run:

```powershell
cargo run | tee output.txt
```

(run main.rs and save the output in output.txt)
  

#### Explanation of the Output

##### Understanding the Code
This Rust program parses a FASTQ file using the bio crate, counts the number of reads, and calculates the total number of bases across all reads. It stores read lengths in a HashMap where keys are read IDs, and values are their sequence lengths.

##### Key Code Components

###### 1. Reading a FASTQ File

```rust
let reader = fastq::Reader::from_file("reads.fastq")?;
```

* Uses the bio crate to open and read a FASTQ file (reads.fastq).
* The ? operator propagates errors if the file cannot be read.

###### 2. Storing Read Lengths in a HashMap

```rust
let mut length_map: HashMap<String, usize> = HashMap::new();
```

* A HashMap stores the read ID as the key and its length as the value.

###### 3. Iterating Over FASTQ Records

```rust
for record_result in reader.records() {
    let record = record_result?;
    let seq_id = record.id().to_string();
    let seq_length = record.seq().len();
    length_map.insert(seq_id, seq_length);
}
```

* Reads each FASTQ record.
* Extracts the sequence ID (record.id()).
* Calculates the sequence length (record.seq().len()).
* Stores the data in length_map.

###### 4. Computing Total Reads and Bases

```rust
let total_reads = length_map.len();
let total_bases: usize = length_map.values().sum();
```

* Total reads = number of entries in length_map.
* Total bases = sum of all sequence lengths.

###### 5. Optional: Checking a Specific Read Length

```rust
if let Some(length) = length_map.get("my_read_id") {
    println!("Read 'my_read_id' has length: {}", length);
}
```

* If "my_read_id" exists in the FASTQ file, it retrieves and prints its length.

###### Expected Output

```sh
Total reads: 54892
Total bases in all reads: 12297218
```

* 54892 reads were found in the FASTQ file.
* The total number of nucleotides across all reads is 12,297,218 bases.

#### Conclusion
This program efficiently counts the number of reads and calculates the total bases in a FASTQ file.
Using a HashMap allows quick lookups of individual read lengths.
The approach is scalable and can process large sequencing datasets efficiently.
Future improvements could include GC-content calculation, length distribution analysis, or filtering reads based on criteria.
















