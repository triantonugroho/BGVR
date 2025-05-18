nextflow.enable.dsl=2

process CALL_VARIANTS {
  tag "$chrom"
  
  input:
    path bam
    val chrom
  
  output:
    path "${chrom}.vcf.bgz"
    
  // Disable container by default since it might not be available
  // Uncomment if you have the container available
  // container 'ghcr.io/acme/rust_variants:2.0'
  
  script:
    """
    echo "Would run: rust_variants --config ../data/pipeline.toml call --bam $bam --out_vcf ${chrom}.vcf.bgz --region $chrom"
    # Create dummy output file
    dd if=/dev/urandom of=${chrom}.vcf.bgz bs=1k count=5
    """
}

workflow {
  // Use absolute or relative paths to the data directory
  bam_ch = channel.fromPath("../data/sample.bam")
  chroms_ch = channel.fromPath("../data/chroms.txt").splitText().map{it.trim()}
  
  CALL_VARIANTS(bam_ch, chroms_ch)
}