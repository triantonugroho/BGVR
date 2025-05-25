#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reads = 'data/reads/*_{1,2}.fastq.gz'
params.outdir = 'results'
params.rust_tool = './target/release/rust_expression_tool'

log.info """
=== Simple RNA-seq Test Pipeline ===
reads: ${params.reads}
outdir: ${params.outdir}
"""

Channel
    .fromFilePairs(params.reads, checkIfExists: true)
    .set { read_pairs_ch }

process SIMPLE_QC {
    tag "$sample_id"
    publishDir "${params.outdir}/qc", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_qc.txt")
    
    script:
    """
    echo "=== QC for ${sample_id} ===" > ${sample_id}_qc.txt
    echo "Date: \$(date)" >> ${sample_id}_qc.txt
    echo "Files: ${reads}" >> ${sample_id}_qc.txt
    
    for file in ${reads}; do
        reads_count=\$((\$(zcat \$file | wc -l) / 4))
        echo "  \$file: \$reads_count reads" >> ${sample_id}_qc.txt
    done
    """
}

process GENERATE_EXPRESSION {
    tag "$sample_id"
    publishDir "${params.outdir}/expression", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path("${sample_id}_expression.tsv")
    
    script:
    """
    echo -e "sample_id\\tgene_id\\tnormalized_count\\traw_count" > ${sample_id}_expression.tsv
    
    for gene in ENSG00000001 ENSG00000002 ENSG00000003 ENSG00000004 ENSG00000005; do
        raw_count=\$(shuf -i 100-5000 -n 1)
        normalized=\$(echo "scale=2; \$raw_count * 20" | bc -l)
        echo -e "${sample_id}\\t\$gene\\t\$normalized\\t\$raw_count" >> ${sample_id}_expression.tsv
    done
    """
}

process COMBINE_AND_ANALYZE {
    publishDir "${params.outdir}/final", mode: 'copy'
    
    input:
    path(expression_files)
    
    output:
    path("combined_matrix.tsv")
    path("analysis_summary.txt")
    
    script:
    """
    echo -e "sample_id\\tgene_id\\tnormalized_count\\traw_count" > combined_matrix.tsv
    cat ${expression_files} | grep -v "sample_id" >> combined_matrix.tsv
    
    echo -e "sample_id\\tcondition\\treplicate\\tbatch" > sample_metadata.tsv
    cut -f1 combined_matrix.tsv | sort | uniq | tail -n +2 | while read sample; do
        if [[ "\$sample" == *"control"* ]]; then
            condition="control"
        else
            condition="treatment"
        fi
        replicate=\$(echo "\$sample" | grep -o '[0-9]\\+\$' || echo "1")
        echo -e "\$sample\\t\$condition\\t\$replicate\\tbatch1" >> sample_metadata.tsv
    done
    
    if [ -f "${params.rust_tool}" ]; then
        echo "Running Rust analysis..."
        if [ -f "${projectDir}/data/gene_info.tsv" ]; then
            cp ${projectDir}/data/gene_info.tsv ./
        fi
        
        ${params.rust_tool} summarize \\
            --input combined_matrix.tsv \\
            --output analysis_summary.txt \\
            --threshold 1.0
    else
        echo "=== Analysis Summary ===" > analysis_summary.txt
        echo "Generated: \$(date)" >> analysis_summary.txt
        echo "Total genes: \$(tail -n +2 combined_matrix.tsv | cut -f2 | sort | uniq | wc -l)" >> analysis_summary.txt
        echo "Total samples: \$(tail -n +2 combined_matrix.tsv | cut -f1 | sort | uniq | wc -l)" >> analysis_summary.txt
        echo "Total records: \$(tail -n +2 combined_matrix.tsv | wc -l)" >> analysis_summary.txt
    fi
    """
}

workflow {
    SIMPLE_QC(read_pairs_ch)
    GENERATE_EXPRESSION(read_pairs_ch)
    COMBINE_AND_ANALYZE(GENERATE_EXPRESSION.out.collect())
}

workflow.onComplete {
    log.info "Pipeline completed! Status: ${workflow.success ? 'SUCCESS' : 'FAILED'}"
}
