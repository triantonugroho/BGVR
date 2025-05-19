#!/usr/bin/env nextflow

/*
 * nextflow.config could set nextflow.enable.dsl=2
 * For example:
 *   nextflow.enable.dsl=2
 *   process.executor = 'local'
 *
 * This pipeline demonstrates pinned dependencies, data versioning,
 * environment variables, and numeric parameters for reproducibility.
 */

// ---------------- PARAMETERS ----------------
params.input_data       = 'dataset_v1.0'    // Version label or path
params.rust_version     = '1.85.0'          // Pinned Rust compiler version
params.samtools_version = '1.19.2'          // Pinned Samtools version
params.numeric_param    = 21                // Example numeric/hyperparameter
params.env_variable     = 'MYAPP_DEBUG=1'   // Example environment variable
params.publish_results  = 'publish_output'  // Directory to save final artifacts

// Step 1: Compile a small Rust application using the pinned Rust version
process BUILD_RUST_APP {
    publishDir "${params.publish_results}", mode: 'copy'

    input:
    val data

    output:
    path "myapp"

    script:
    """
    echo "== Building with Rust version: ${params.rust_version}"
    echo "== Input dataset label: ${data}"
    echo "== Environment variable: ${params.env_variable}"
    export ${params.env_variable}

    rustc --version
    echo 'fn main() { println!("Hello from pinned Rust!"); }' > app.rs
    rustc app.rs -o myapp
    """
}

// Step 2: Run the compiled Rust app alongside Samtools
process RUN_TOOL {
    publishDir "${params.publish_results}", mode: 'copy'

    input:
    path app
    val numeric_param

    output:
    path "results.txt"

    script:
    """
    echo "== Using pinned Samtools version: ${params.samtools_version}"
    samtools --version

    echo "Running the compiled Rust app..."
    ./myapp

    echo "Numeric parameter is: ${numeric_param}"
    echo "Results for param = ${numeric_param}" > results.txt
    """
}

// Main Workflow Block
workflow {
    data_ch = Channel.value(params.input_data)
    numeric_ch = Channel.value(params.numeric_param)
    
    built_app = BUILD_RUST_APP(data_ch)
    results_ch = RUN_TOOL(built_app, numeric_ch)
}