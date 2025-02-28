process runRust {
    output:
        path 'output.txt'

    script:
    """
    cargo run > output.txt 2>&1
    """
}

workflow {
    runRust()
}
