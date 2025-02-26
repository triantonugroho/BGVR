params.input_path = 'SRR11192680'
params.output_path = 'de_bruijn_graph.txt'
params.k = 21
params.sra_toolkit_path = '/mnt/c/SRAToolkit/bin/fasterq-dump.exe'
params.output_dir = '/mnt/c/Users/trian/BGVR/chapter_01/experiment_18_1_2/src' // Tambahkan parameter untuk output direktori

inputPathChannel = Channel.of(params.input_path)

process DOWNLOAD_FASTQ {
    input:
    val input_path

    output:
    file 'reads.fastq'

    script:
    """
    echo "Downloading SRA file ${input_path}..."
    ${params.sra_toolkit_path} ${input_path} --outdir ./ --split-files
    mv ${input_path}_1.fastq reads.fastq
    """
}

process BUILD_DEBRUIJN {
    input:
    file reads_fastq
    val k

    output:
    file params.output_path

    script:
    """
    echo "Building De Bruijn graph..."
    rustc main.rs -o build_de_bruijn
    ./build_de_bruijn --k ${k} --input reads.fastq --output ${params.output_path}
    """
}

workflow {
    reads_fastq = DOWNLOAD_FASTQ(inputPathChannel)
    de_bruijn_graph = BUILD_DEBRUIJN(reads_fastq, params.k)

    // Gunakan parameter output_dir untuk publishDir
    de_bruijn_graph.first().publishDir(params.output_dir, mode: 'copy')
}