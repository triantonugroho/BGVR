use bio::io::fastq; 
// 'bio' crate: provides parsers for FASTQ/FASTA and alignment functions.

use log::{info, error}; 
// 'log' crate: standard logging facade in Rust.
use env_logger; 
// 'env_logger': enables logging via environment variables (e.g., RUST_LOG=info).

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Initialize logging, typically configured via environment variables on HPC
    env_logger::init();

    // Read from a FASTQ file (potentially downloaded via Nextflow process)
    let reader = fastq::Reader::from_file("reads.fastq")
        .map_err(|e| {
            error!("Failed to open reads.fastq: {:?}", e);
            e
        })?;

    // We create an output writer for processed reads
    let mut writer = fastq::Writer::to_file("filtered_reads.fastq")?;

    // For demonstration, we apply a trivial filter: keep only reads with length > 75
    for record_result in reader.records() {
        let record = record_result.map_err(|e| {
            error!("Error reading FASTQ record: {:?}", e);
            e
        })?;

        if record.seq().len() > 75 {
            writer.write_record(&record)?;
        }
    }

    info!("Filtering complete. Output saved to filtered_reads.fastq");
    println!("reads.fastq have successfully filtered with a trivial filter: keep only reads with length > 75 and saved to filtered_reads.fastq");
    Ok(())
}
