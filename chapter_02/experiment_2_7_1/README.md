## 2.7. Data and Computational Foundations

### experiment 2_7_1 

A concise example of downloading a FASTQ file from NCBI using SRAToolkit), such as NCBI. By making an HTTP request to a specified URL, the program retrieves the raw sequencing data and writes it to a local

The following Rust code provides a concise example of downloading a FASTQ file from an online repository, such as NCBI by using SRAToolkit. The program retrieves the raw sequencing data and writes it to a local file. This capability is valuable for pipelines that integrate external data sources, enabling researchers to automate the acquisition of public genomic resources directly within their computational workflows.

#### Project Structure:

```plaintext
experiment_2_7_1/
├── Cargo.toml                     # Rust project configuration and dependencies
└── src/
    ├── main.rs                    # Main Rust script containing program logic
    ├── SRR11192680_1.rar          # SRR11192680_1.fastq compressed in RAR file
    ├── SRR11192680_2.rar          # SRR11192680_2.fastq compressed in RAR file
    └── output.txt                 # Output file
```

#### Cargo.toml

```toml
[package]
name = "experiment_2_7"
version = "0.1.0"
edition = "2021"

[dependencies]

```

#### How to run:

```powershell
cargo run | tee output.txt
```

(run main.rs and save the output in output.txt)
  

SRAToolkit : fasterq-dump tool need to be installed 

#### Explanation of the Output

This program downloads a FASTQ file from the SRA (Sequence Read Archive) using the fasterq-dump tool. Below is a breakdown of the process and the output:

##### 1. Process Overview

* The program specifies an SRA accession number:

  SRR11192680
  
  * This represents 16S rRNA sequencing data from a rectal sample of a prostate cancer patient.

* The program executes the fasterq-dump command using Rust’s Command API to download the sequencing data from NCBI SRA.

* The output is saved to the directory:

  C:\Users\trian\BGVR\chapter_01\experiment_2_7_1\src

##### 2. Breakdown of the Output

FASTQ file has been successfully downloaded to 'C:\Users\trian\BGVR\chapter_01\experiment_17_1\src'.

* This confirms that the fasterq-dump command executed successfully.
* The downloaded FASTQ file contains raw sequencing reads from the SRA dataset (SRR11192680).
* The program also logs the success message in output.txt.

##### 3. Key Insights

* Automated Data Retrieval: The program automates the download of sequencing data from NCBI SRA.
* FASTQ Format: The downloaded file is in FASTQ format, commonly used in bioinformatics for storing raw nucleotide sequences.
* Error Handling: If the download failed, the program would log an error message in output.txt.

#### Conclusion

The program successfully retrieved the FASTQ sequencing file from NCBI SRA using fasterq-dump. This data can now be used for further bioinformatics analysis, such as quality control, taxonomic classification, or microbial community analysis.

