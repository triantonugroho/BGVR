#!/usr/bin/env nextflow

/*
 * Bioinformatics RNA-seq Processing Pipeline
 * Supports pseudo-alignment with Rust tool
 */

nextflow.enable.dsl = 2

// Help message
def helpMessage() {
    log.info"""
    =======================================================
    Bioinformatics RNA-seq Processing Pipeline
    =======================================================
    
    Usage:
    nextflow run main.nf --method pseudo --reads 'data/*.fastq*'
    
    Mandatory arguments:
      --reads [file]                Path to input reads
      --method [str]                Analysis method: 'pseudo'
    
    Pseudo-alignment options:
      --kmer_index [file]           Path to k-mer index JSON file
      --kmer_length [int]           K-mer length (default: 31)
      --min_read_length [int]       Minimum read length (default: 50)
    
    Other options:
      --outdir [path]               Output directory (default: 'results')
      --skip_qc                     Skip quality control
      --skip_multiqc                Skip MultiQC report
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// Validate parameters
if (!params.reads) {
    log.error "Please provide input reads with --reads"
    exit 1
}

if (!params.method.matches('pseudo')) {
    log.error "Currently only 'pseudo' method is supported"
    exit 1
}

if (!file(params.kmer_index).exists()) {
    log.error "K-mer index file not found: ${params.kmer_index}"
    exit 1
}

// Header
log.info """\
    =======================================================
    RNA-seq Pseudo-alignment Pipeline
    =======================================================
    reads        : ${params.reads}
    method       : ${params.method}
    outdir       : ${params.outdir}
    kmer_index   : ${params.kmer_index}
    kmer_length  : ${params.kmer_length}
    skip_qc      : ${params.skip_qc}
    =======================================================
    """.stripIndent()

// Processes
process RUST_PSEUDOALIGN {
    tag "Pseudo-align $meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/pseudoalignment", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    path kmer_index

    output:
    tuple val(meta), path("*.tsv"), emit: quantification
    tuple val(meta), path("*.log"), emit: log
    path "versions.yml", emit: versions

    script:
    def prefix = meta.id
    """
    # Check if rust_pseudoalign is available, if not use cargo run
    if command -v rust_pseudoalign &> /dev/null; then
        RUST_CMD="rust_pseudoalign"
    else
        RUST_CMD="${projectDir}/target/release/rust_pseudoalign"
    fi
    
    # If binary doesn't exist, try to build it
    if [ ! -f "\$RUST_CMD" ]; then
        echo "Building Rust pseudo-alignment tool..."
        cd ${projectDir}
        cargo build --release
        RUST_CMD="${projectDir}/target/release/rust_pseudoalign"
    fi
    
    # Run pseudo-alignment
    export RUST_LOG=info
    
    \$RUST_CMD \\
        --index ${kmer_index} \\
        --reads ${reads} \\
        --output ${prefix}_quantification.tsv \\
        --kmer-length ${params.kmer_length} \\
        --min-read-length ${params.min_read_length} \\
        --threads ${task.cpus} \\
        --verbose \\
        2>&1 | tee ${prefix}_pseudoalign.log

    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rust_pseudoalign: \$(\$RUST_CMD --version 2>/dev/null || echo "0.1.0")
    END_VERSIONS
    """
}

process CUSTOM_DUMPSOFTWAREVERSIONS {
    label 'process_single'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    path versions

    output:
    path "software_versions.yml", emit: yml
    path "software_versions_mqc.yml", emit: mqc_yml

    script:
    """
    cat $versions > software_versions.yml
    
    cat <<-END_VERSIONS > software_versions_mqc.yml
    id: 'software_versions'
    section_name: 'Software Versions'
    section_href: 'https://github.com/your-org/bio-pseudo-align'
    plot_type: 'html'
    description: 'are collected at run time from the software output.'
    data: |
        <dl class="dl-horizontal">
            <dt>Nextflow</dt><dd>${workflow.nextflow.version}</dd>
            <dt>Pipeline</dt><dd>${workflow.manifest.version}</dd>
        </dl>
    END_VERSIONS
    """
}

// Main workflow
workflow {
    ch_versions = Channel.empty()

    // Create input channel from reads
    Channel
        .fromPath(params.reads, checkIfExists: true)
        .map { file -> 
            def meta = [:]
            meta.id = file.getSimpleName()
            return [meta, file]
        }
        .set { ch_reads }

    // Run pseudo-alignment
    kmer_index = Channel.fromPath(params.kmer_index, checkIfExists: true)
    
    RUST_PSEUDOALIGN(ch_reads, kmer_index)
    ch_versions = ch_versions.mix(RUST_PSEUDOALIGN.out.versions)
    
    // Software versions
    CUSTOM_DUMPSOFTWAREVERSIONS(ch_versions.collectFile(name: 'collated_versions.yml'))
}

// Completion
workflow.onComplete {
    log.info """
    =======================================================
    Pipeline execution summary
    =======================================================
    Completed at : ${workflow.complete}
    Duration     : ${workflow.duration}
    Success      : ${workflow.success}
    Results Dir  : ${file(params.outdir).toAbsolutePath()}
    Work Dir     : ${workflow.workDir.toAbsolutePath()}
    Exit status  : ${workflow.exitStatus}
    =======================================================
    """.stripIndent()
}
