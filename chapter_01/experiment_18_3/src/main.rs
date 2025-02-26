use std::fs::File;
use std::io::{self, Write};
use rust_htslib::bam;
use rust_htslib::bam::Read;
use ndarray::Array1;
use nalgebra::DVector;
use rayon::prelude::*;

// This function processes a BAM file in parallel using Rayon
fn parallel_bam_processing(bam_path: &str) -> (Array1<f64>, DVector<f64>) {
    let mut bam_reader = bam::Reader::from_path(bam_path)
        .expect("Failed to open BAM file");
    
    // Example placeholder: compute a coverage vector, parallelized
    let coverage: Vec<f64> = bam_reader.records().par_bridge().map(|record_res| {
        let record = record_res.expect("Invalid BAM record");
        // Simple placeholder for coverage computation
        record.seq_len() as f64
    }).collect();

    let nd_coverage = Array1::from_vec(coverage);
    let na_coverage = DVector::from_fn(nd_coverage.len(), |i, _| nd_coverage[i]);

    (nd_coverage, na_coverage)
}

// Function to write output to a file
fn write_output_to_file(filename: &str, array_cov: &Array1<f64>, dvec_cov: &DVector<f64>) -> io::Result<()> {
    let mut file = File::create(filename)?;
    
    // Write array coverage length to file
    writeln!(file, "Array coverage length: {}", array_cov.len())?;
    
    // Write DVector coverage length to file
    writeln!(file, "DVector coverage length: {}", dvec_cov.len())?;

    // Optionally, write the actual coverage values
    writeln!(file, "Coverage values (Array1): {:?}", array_cov)?;
    writeln!(file, "Coverage values (DVector): {:?}", dvec_cov)?;

    Ok(())
}

fn main() {
    let path = "/mnt/c/Users/trian/BGVR/chapter_01/experiment_18_3/src/reads.bam";
    let (array_cov, dvec_cov) = parallel_bam_processing(path); // Menggunakan variabel path

    // Specify the output file path
    let output_file = "output.txt";

    // Write the results to the output file
    if let Err(e) = write_output_to_file(output_file, &array_cov, &dvec_cov) {
        eprintln!("Error writing to file: {}", e);
    } else {
        println!("Output written to {}", output_file);
    }
}