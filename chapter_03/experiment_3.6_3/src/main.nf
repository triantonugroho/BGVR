#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reads          = '/mnt/c/Users/trian/BGVR/chapter_03/experiment_36_3/src/data/sample.fastq'
params.kmerSize       = 21
params.numCpus        = 8
params.outputDir      = 'results'
params.abyssExecutable = 'abyss-pe'  // Or 'mpirun -n 4 abyss-pe' for MPI

process preprocessReads {
    input:
    file fastq

    output:
    file "filtered_*"

    script:
    """
    /mnt/c/Users/trian/BGVR/chapter_03/experiment_36_3/target/release/rust_preprocess --input ${fastq} --output filtered_${fastq}
    """
}

process runAbyss {
    input:
    path filtered_reads

    output:
    file "contigs.fa"

    script:
    """
    ${params.abyssExecutable} k=${params.kmerSize} B=128M name=assembly in='${filtered_reads}'
mv assembly-contigs.fa contigs.fa
    """
}

process assembleStats {
    input:
    file contigs

    output:
    file "stats.txt"

    script:
    """
    abyss-fac ${contigs} > stats.txt
    """
}

workflow {
    reads_ch = Channel.fromPath(params.reads)

    filtered_reads_ch = preprocessReads(reads_ch) 
    assembled_contigs_ch = runAbyss(filtered_reads_ch.collect()) 
    assembleStats(assembled_contigs_ch)
}
