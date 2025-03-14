#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.synthetic_fastq = '/mnt/c/Users/trian/BGVR/chapter_04/experiment_41/src/synthetic_reads.fastq'
params.rust_project_path = '/mnt/c/Users/trian/BGVR/chapter_04/experiment_41/src'
params.output_dir = './results'

workflow {
    runAnalysis()
}

process runAnalysis {
    publishDir params.output_dir, mode: 'copy'
    
    output:
        path 'pwm_results.txt', optional: true
        path 'mrf_results.txt', optional: true
    
    script:
    """
    #!/bin/bash
    set -e
    
    echo "Current directory: \$(pwd)"
    echo "Starting compilation..."
    
    # Kompilasi program Rust
    cd ${params.rust_project_path}
    cargo build --release
    
    echo "Compilation complete. Running program..."
    
    # Jalankan program yang sudah dikompilasi
    ./target/release/experiment_41 \
        ${params.synthetic_fastq} \
        ${params.rust_project_path}/pwm_results.txt \
        ${params.rust_project_path}/mrf_results.txt
    
    # Salin hasil kembali ke direktori kerja
    cp ${params.rust_project_path}/pwm_results.txt ./
    cp ${params.rust_project_path}/mrf_results.txt ./
    
    echo "Analysis complete."
    """
}