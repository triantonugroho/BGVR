nextflow.enable.dsl=2

workflow {  
    mainFlow() 
}

workflow mainFlow {
    def read_files = Channel.fromPath('reads/*.fastq')  

    def aligned_bams = RUST_ALIGN(read_files)  
    def coverage_tables = RUST_PARSE(aligned_bams)  
    def ml_output = RUST_ML(coverage_tables)  

    ml_output.view()
}

process RUST_ALIGN {
    container 'my-rust-alignment:1.0.0'
    executor 'local'

    input:
    file readFile

    output:
    file 'output.bam'  // ✅ No 'into' here

    script:
    """
    my_rust_align --input ${readFile} --output output.bam
    """
}

process RUST_PARSE {
    container 'my-rust-parse:1.0.0'
    executor 'local'

    input:
    file bamFile

    output:
    file 'coverage.csv'  // ✅ No 'into' here

    script:
    """
    my_rust_parse --bam ${bamFile} --out coverage.csv
    """
}

process RUST_ML {
    container 'my-rust-ml:cuda-1.0.0'
    executor 'local'

    input:
    file covFile

    output:
    file 'ml_results.json'  // ✅ No 'into' here

    script:
    """
    my_rust_ml --inputs ${covFile.join(' ')} --output ml_results.json
    """
}
