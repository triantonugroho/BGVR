use std::process::{Command, Stdio};
use std::error::Error;
use std::fs;
use bincode::{serialize, deserialize}; // ✅ Import serialize & deserialize
use serde::{Serialize, Deserialize};

#[derive(Serialize, Deserialize, Debug)]
struct PartialSuffixArray {
    data: Vec<u8>,
}

// Function to serialize data
fn serialize_data(partial_sa: &PartialSuffixArray) -> Vec<u8> {
    bincode::serialize(partial_sa).expect("Failed to serialize PartialSuffixArray")
}

// Function to deserialize data
fn deserialize_data(recv_bytes: &[u8]) -> PartialSuffixArray {
    bincode::deserialize(recv_bytes).expect("Failed to deserialize PartialSuffixArray")
}

fn run_nextflow_process(process: &str) -> Result<String, Box<dyn Error>> {
    let script_path = r"C:\Users\trian\BGVR\chapter_03\experiment_36_2\src\main.nf";

    // Check if the script exists
    if !fs::metadata(script_path).is_ok() {
        return Err(format!("Nextflow script not found: {}", script_path).into());
    }

    // Run the Nextflow command and capture its output
    println!("Running Nextflow process: {}", process);
    let output = Command::new("nextflow")
        .arg("run")
        .arg(script_path)
        .arg("-process")
        .arg(process)
        .arg("-with-docker") // Use Docker if specified
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .output()?; // Capture the output

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        let stdout = String::from_utf8_lossy(&output.stdout);
        return Err(format!(
            "Nextflow failed for process `{}`.\nstdout: {}\nstderr: {}",
            process, stdout, stderr
        )
        .into());
    }

    // Return expected output file names based on process
    match process {
        "BUILD_INDEX" => Ok("partial_fm_indexes.json".to_string()),
        "ALIGN_READS" => Ok("alignment_results.json".to_string()),
        "SUMMARIZE" => Ok("final_summary.json".to_string()),
        _ => Err(format!("Unknown Nextflow process: {}", process).into()),
    }
}

fn example_serialization() -> Result<(), Box<dyn Error>> {
    let partial_sa = PartialSuffixArray { data: vec![1, 2, 3, 4, 5] };

    // Serialize data
    let serialized_data = serialize_data(&partial_sa);
    println!("Serialized data: {:?}", serialized_data);

    // Deserialize data
    let deserialized_data = deserialize_data(&serialized_data);
    println!("Deserialized data: {:?}", deserialized_data);

    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    // Step 1: Build index
    println!("Building index...");
    let index_path = run_nextflow_process("BUILD_INDEX")?;
    println!("Index built at: {}", index_path);

    // Step 2: Align reads
    println!("Aligning reads...");
    let alignment_result = run_nextflow_process("ALIGN_READS")?;
    println!("Alignment results saved in: {}", alignment_result);

    // Step 3: Summarize results
    println!("Summarizing results...");
    let summary_result = run_nextflow_process("SUMMARIZE")?;
    println!("Summary saved in: {}", summary_result);

    println!("Pipeline complete.");

    // Run serialization example
    example_serialization()?;

    Ok(())
}
