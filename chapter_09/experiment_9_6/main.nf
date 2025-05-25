#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.input_matrix = "data/raw/sparse_counts.tsv"
params.output_dir = "results"

log.info """
================================
Simple Single-Cell Pipeline
================================
Input: ${params.input_matrix}
Output: ${params.output_dir}
================================
"""

process VALIDATE_INPUT {
    publishDir "${params.output_dir}/validation", mode: 'copy'
    memory '1.GB'
    cpus 1
    
    input:
    path input_matrix
    
    output:
    path "validated_matrix.tsv", emit: matrix
    path "validation_report.txt", emit: report
    
    script:
    """
    python3 ${projectDir}/scripts/validate_input.py ${input_matrix} ./
    """
}

process RUST_ANALYSIS {
    publishDir "${params.output_dir}/coordinates", mode: 'copy'
    memory '2.GB'
    cpus 2
    
    input:
    path matrix
    
    output:
    path "cell_coordinates.tsv", emit: coordinates
    
    script:
    """
    # Find and run Rust binary
    if [ -f "${projectDir}/target/release/scrna-analyzer" ]; then
        ${projectDir}/target/release/scrna-analyzer --input ${matrix} --output cell_coordinates.tsv
    else
        echo "ERROR: Rust binary not found at ${projectDir}/target/release/scrna-analyzer"
        exit 1
    fi
    """
}

workflow {
    input_ch = Channel.fromPath(params.input_matrix, checkIfExists: true)
    
    VALIDATE_INPUT(input_ch)
    RUST_ANALYSIS(VALIDATE_INPUT.out.matrix)
}

workflow.onComplete {
    log.info "Pipeline completed: ${workflow.success}"
}
