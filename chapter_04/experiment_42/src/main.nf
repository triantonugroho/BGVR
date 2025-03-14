#!/usr/bin/env nextflow

/*
 * Usage example:
 *   nextflow run main.nf --num_genes 1000 --num_samples 50 --output_file partial_adjacency.bin
 */

params.num_genes    = (params.num_genes    ?: 1000)
params.num_samples  = (params.num_samples  ?: 50)
params.output_file  = (params.output_file  ?: "partial_adjacency.bin")

workflow {
    // Kirim langsung ke proses (gunakan file() untuk output)
    buildAdjacencyMatrix(params.num_genes, params.num_samples, file(params.output_file))
}

process buildAdjacencyMatrix {
    input:
        val genes
        val samples
        path outputFile

    output:
        path outputFile

    script:
    """
    set -e
    set -x
    
    # Masuk ke direktori kode sumber
    cd /mnt/c/Users/trian/BGVR/chapter_04/experiment_42
    
    # Compile kode Rust
    cargo build --release
    
    # Jalankan hasil binary dengan output di direktori kerja Nextflow
    /mnt/c/Users/trian/BGVR/chapter_04/experiment_42/target/release/experiment_42 \
        --num-genes ${genes} \
        --num-samples ${samples} \
        --output ${outputFile}
    """
}
