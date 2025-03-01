nextflow.enable.dsl=2

workflow {
    fastq_files = Channel.fromPath('reads/*.fastq')
    trimmed_fastqs = RUST_TRIM(fastq_files)
    trimmed_fastqs.view()
}

process RUST_TRIM {
    container 'my-rust-trimmer:1.0.0'
    executor 'local'   // ✅ Use local execution instead of SLURM

    input:
    file readFile

    output:
    file "trimmed_*.fastq"

    script:
    """
    /mnt/c/Users/trian/BGVR/chapter_02/experiment_26_1/target/release/my_rust_trim --input ${readFile} --output trimmed_${readFile.simpleName}.fastq
    """
}
