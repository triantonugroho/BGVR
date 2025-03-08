#!/usr/bin/env nextflow

params.ref   = file("${params.ref ?: 'reference.fa'}")
params.reads = file("${params.reads ?: 'reads/*.fq'}")

process BUILD_INDEX {
    container 'myregistry.io/bioinformatics/genome-reference-builder:latest'
    
    input:
    path params.ref

    output:
    path "partial_fm_indexes.json"

    script:
    """
    genome-reference-builder ${params.ref} --output partial_fm_indexes.json --size 1000000
    """
}

process ALIGN_READS {
    container 'myregistry.io/bioinformatics/aligner:latest'
    
    input:
    path read_file
    path partial_index

    output:
    path "alignment_results.json"

    script:
    """
    aligner ${read_file} --index ${partial_index} --output alignment_results.json
    """
}

process SUMMARIZE {
    container 'myregistry.io/bioinformatics/summarizer:latest'
    
    input:
    path results

    output:
    path "final_summary.json"

    script:
    """
    summarizer ${results.join(' ')} --output final_summary.json
    """
}

workflow {
    // Generate the index
    index = BUILD_INDEX(params.ref)
    
    // Collect all FASTQ files
    read_files = Channel.fromPath(params.reads, checkIfExists: true)

    // Run alignment separately for each read
    aligned = ALIGN_READS(read_file: read_files, partial_index: index)

    // Collect all results before summarizing
    summarized = SUMMARIZE(results: aligned.collect())

    summarized.view()
}
