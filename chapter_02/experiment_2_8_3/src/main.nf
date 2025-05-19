process runRust {
    output:
        path 'output.txt'

    script:
    """
    cargo run > output.txt
    """
}

workflow {
    runRust()
}
