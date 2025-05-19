use std::process::{Command, Stdio};
use std::error::Error;
use std::fs::File;
use std::io::Write;

fn main() -> Result<(), Box<dyn Error>> {
    // SRA accession number
    let sra_accession = "SRR11192680"; // 16S rRNA seq of Homo sapiens: prostate cancer patient rectal sample

    // Paths for the downloaded FASTQ file and output log
    let output_dir = "C:\\Users\\trian\\BGVR\\chapter_01\\experiment_17_1\\src";
    let output_path = format!("{}\\output.txt", output_dir);
    let fasterq_dump_path = "C:\\SRAToolkit\\bin\\fasterq-dump.exe"; // Full path to fasterq-dump.exe

    // Open output file
    let mut output_file = File::create(&output_path)?;

    // Construct the fasterq-dump command with the full path
    let output = Command::new(fasterq_dump_path)
        .arg(sra_accession)
        .arg("--outdir")
        .arg(output_dir)
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .output()?;

    // Check if the command was successful
    if !output.status.success() {
        let error_message = String::from_utf8_lossy(&output.stderr);
        writeln!(output_file, "Error: {}", error_message)?;
        return Err("Failed to download FASTQ file using fasterq-dump.".into());
    }

    // Write success message to output file
    writeln!(output_file, "FASTQ file has been successfully downloaded to '{}'.", output_dir)?;
    println!("FASTQ file has been successfully downloaded to '{}'.", output_dir);

    Ok(())
}
