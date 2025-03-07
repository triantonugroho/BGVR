## 3.1. Introduction to Data Structures and Algorithms

### experiment_31_6

In this example code, the program reads each FASTA record, computes basic statistics such as GC content, and then outputs the aggregated results in JSON format. The snippet demonstrates a straightforward yet robust approach, where the code loops over each record from the FASTA file as needletail parses it and converts that record into a Rust structure. The structure stores sequence length, GC content, and the record identifier. Because needletail returns each record as a copy-on-write slice (Cow<'_, [u8]>), the code simply borrows the slice for computations instead of allocating new buffers. Once all records have been collected into a vector, they are written to a file in JSON form using serde_json.

The program calculates each sequence’s GC content by counting the occurrences of G and C nucleotides, then dividing by the total sequence length. If n is the total number of sequences and \ell is the average length of a sequence, the computation of GC content is O(\ell) for each record. Summing over all records, the overall time complexity becomes O(n \times \ell). For large data sets, this linear time approach is generally efficient, assuming that reading from disk is also managed in a streaming fashion. If one were to parallelize the GC content calculation, it could further reduce wall-clock time on machines with multiple CPU cores, but the asymptotic complexity remains O(n \times \ell).

#### Files contents:
* experiment_31_2/Cargo.toml (Cargo.toml file for dependencies)
* src/
  * main.rs (rust script)
  * example.fasta (fasta file)
  * output.txt (output file)
  * output.json

#### How to run:

run in powershell:

```powershell
cargo run | tee output.txt
```

(run main.rs and save the output in output.txt)
  
#### [dependencies]

```toml
needletail = "0.6.3"  # For parsing FASTA/FASTQ files
serde = { version = "1.0", features = ["derive"] }  # For serializing/deserializing data
serde_json = "1.0"  # For writing JSON output
```

#### Explanation of the Output:
The program processes a FASTA file containing genomic sequence data and performs the following tasks:

##### 1. Reads the FASTA file:

* The program expects a FASTA or FASTQ file as input. It either takes the file path from the command-line argument or defaults to "example.fasta".
* It uses the needletail crate’s parse_fastx_file function to read sequences from the file.

##### 2. Extracts Sequence Information:

* For each sequence in the FASTA file, the program calculates:
  * GC content (the percentage of G and C nucleotides in the sequence).
  * Sequence length.
  * Sequence ID (the name or identifier of the sequence).
* The information for each sequence is stored in a struct SeqRecord, which includes:
  * id: The sequence ID.
  * seq_len: The length of the sequence.
  * gc_content: The GC content percentage of the sequence.

##### 3. GC Content Calculation:

* The program defines the calc_gc_content function, which counts the number of G and C nucleotides in the sequence and calculates the ratio of these bases to the total length of the sequence.

##### 4. Writing Results to Output:

* After processing all sequences in the FASTA file, the program serializes the sequence records into a JSON file (output.json) using the serde_json library.
* It also prints the number of sequences processed.

##### Output Breakdown:

###### 1. Console Output:

```rust
Processed 1 sequences
```

* This message indicates that the program has successfully processed 1 sequence from the FASTA file.
* This suggests that the input FASTA file contains only 1 sequence or that only 1 sequence was successfully parsed.

##### 2. output.json Content:
```json
[
  {
    "id": "Synthetic Genome",
    "seq_len": 10000000,
    "gc_content": 0.4999007
  }
]
```

* This JSON file contains an array of one object representing the sequence data:
  * id: "Synthetic Genome" - This is the identifier of the sequence.
  * seq_len: 10000000 - The sequence length is 10 million nucleotides.
  * gc_content: 0.4999007 - The GC content of the sequence is approximately 50%.

#### Conclusion:
* The program successfully processed a single sequence from the FASTA file (example.fasta or a provided file).
* The sequence has a length of 10 million bases, and its GC content is about 50%.
* The sequence's details were serialized into a JSON file (output.json), which includes the sequence's ID, length, and GC content.
* The message "Processed 1 sequences" indicates that only one sequence was present or successfully processed.

#### Key Takeaways:
* The program can process FASTA files, calculate GC content, and output the results in JSON format.
* If multiple sequences were present in the input file, the program would process each sequence and include all their details in the output JSON.
* The current output suggests that only one sequence was processed, which could be due to either the input file containing just one sequence or potential issues in parsing the file (though the latter seems unlikely since the sequence was processed correctly).



