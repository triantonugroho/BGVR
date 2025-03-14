#!/usr/bin/env nextflow

/*
 * Usage example:
 *   nextflow run pipeline.nf --num_genes 1000 --num_samples 50
 */

params.num_genes    = (params.num_genes    ?: 1000)
params.num_samples  = (params.num_samples  ?: 50)
params.output_file  = (params.output_file  ?: "partial_adjacency.bin")

process buildAdjacencyMatrix {
    input:
    val genes      = params.num_genes
    val samples    = params.num_samples
    val outputFile = params.output_file

    output:
    file outputFile

    """
    # Compile the Rust code
    cd rust_code
    cargo build --release
    
    # Run the compiled binary with the specified parameters
    ./target/release/parallel_correlation \\
        --num-genes ${genes} \\
        --num-samples ${samples} \\
        --output ../${outputFile}
    """
}

workflow {
    main:
        buildAdjacencyMatrix()
}