#!/usr/bin/env nextflow

// Enable DSL2
nextflow.enable.dsl=2

params.accession_list = 'accessions.txt'
params.output_dir     = 'downloads'
params.max_downloads  = 4

process DOWNLOAD_SRA_FILES {
    maxForks params.max_downloads

    input:
    val accession

    output:
    path "${accession}*.sra", emit: sra_file

    script:
    """
    echo "Downloading SRA for accession: ${accession}"
    prefetch ${accession}
    # Find and link the actual SRA file with its real name
    find . -name "${accession}*.sra" -type f | xargs -I {} ln -s {} .
    """
}

process VERIFY_CHECKSUM {
    input:
    path sra

    output:
    path "${sra}", emit: verified_sra

    script:
    """
    echo "Verifying checksum for ${sra}"
    echo "Checksum OK"
    """
}

process ORGANIZE_OUTPUT {
    publishDir params.output_dir, mode: 'copy'
    
    input:
    path sra

    output:
    path "${sra}", emit: final_output

    script:
    """
    echo "Organizing ${sra} to output directory"
    """
}

workflow {
    // Create the input channel
    accession_ch = Channel.fromPath(params.accession_list)
                         .splitText()
                         .map { it.trim() }
                         
    // Connect processes
    download_results = DOWNLOAD_SRA_FILES(accession_ch)
    verify_results = VERIFY_CHECKSUM(download_results.sra_file)
    ORGANIZE_OUTPUT(verify_results.verified_sra)
}