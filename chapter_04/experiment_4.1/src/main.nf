#!/usr/bin/env nextflow

params.synthetic_fastq = params.synthetic_fastq ?: 'synthetic_reads.fastq'

workflow {
    // Buat channel dari parameter
    def synthetic_fastq_ch = Channel.of(params.synthetic_fastq)

    main:
        buildPWMandMRF(synthetic_fastq_ch)
}

process buildPWMandMRF {
    input:
    value file_in from synthetic_fastq_ch

    output:
    file 'pwm_results.txt'
    file 'mrf_results.txt'

    """
    set -e
    set -x
    
    # Masuk ke direktori kode sumber
    cd src
    
    # Compile kode Rust
    cargo build --release
    
    # Periksa apakah binary hasil build tersedia
    if [ ! -f target/release/experiment_41 ]; then
        echo "Error: Binary target/release/experiment_41 tidak ditemukan!" >&2
        exit 1
    fi
    
    # Jalankan hasil binary dengan path yang benar
    ./target/release/experiment_41 ${file_in} pwm_results.txt mrf_results.txt
    """
}
