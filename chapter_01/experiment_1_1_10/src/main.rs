use std::collections::HashMap;
use std::error::Error;
use bio::io::fastq;

fn main() -> Result<(), Box<dyn Error>> {
    // Open a FASTQ file using the 'bio' crate
    let reader = fastq::Reader::from_file("reads.fastq")?;

    // Create a HashMap to track the lengths of reads keyed by their IDs
    let mut length_map: HashMap<String, usize> = HashMap::new();

    // Iterate over each record in the FASTQ file
    for record_result in reader.records() {
        let record = record_result?;
        let seq_id = record.id().to_string();
        let seq_length = record.seq().len();
        length_map.insert(seq_id, seq_length);
    }

    // Summarize and print the stored information
    let total_reads = length_map.len();
    let total_bases: usize = length_map.values().sum();
    println!("Total reads: {}", total_reads);
    println!("Total bases in all reads: {}", total_bases);

    // Optionally, you can check specific read lengths
    if let Some(length) = length_map.get("my_read_id") {
        println!("Read 'my_read_id' has length: {}", length);
    }

    Ok(())
}