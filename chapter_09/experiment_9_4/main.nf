#!/usr/bin/env nextflow

// Pipeline parameters
params.input = 'normalized_counts.tsv'
params.metadata = 'sample_metadata.tsv'
params.output_dir = 'de_results'
params.control_group = 'Control'
params.treatment_group = 'Treatment'
params.alpha = 0.05
params.min_count = 10.0
params.help = false

// Print help message
if (params.help) {
    log.info """
    Differential Expression Analysis Pipeline
    ========================================
    
    Usage:
        nextflow run main.nf [options]
    
    Required parameters:
        --input             Normalized count file (TSV format) [default: normalized_counts.tsv]
        --metadata          Sample metadata file [default: sample_metadata.tsv]
    
    Optional parameters:
        --output_dir        Output directory [default: de_results]
        --control_group     Name of control group [default: Control]
        --treatment_group   Name of treatment group [default: Treatment]
        --alpha             Significance threshold [default: 0.05]
        --min_count         Minimum mean count threshold [default: 10.0]
        --help              Show this help message
    """
    exit 0
}

// Validate input parameters
if (!file(params.input).exists()) {
    log.error "Input file does not exist: ${params.input}"
    exit 1
}

if (!file(params.metadata).exists()) {
    log.error "Metadata file does not exist: ${params.metadata}"
    exit 1
}

// Print pipeline parameters
log.info """
Differential Expression Analysis Pipeline
=========================================
Input file:           ${params.input}
Metadata file:        ${params.metadata}
Output directory:     ${params.output_dir}
Control group:        ${params.control_group}
Treatment group:      ${params.treatment_group}
Significance level:   ${params.alpha}
Min count threshold:  ${params.min_count}
"""

process BUILD_ANALYZER {
    tag "building_diff_expr_analyzer"
    
    input:
    path cargo_toml, name: 'project_cargo.toml'
    path main_rs, name: 'project_main.rs'
    
    output:
    path "diff-expr-analyzer"
    
    script:
    """
    # Copy source files to correct locations
    mkdir -p src
    cp project_main.rs src/main.rs
    cp project_cargo.toml Cargo.toml
    
    # Build the Rust application
    cargo build --release
    cp target/release/diff-expr-analyzer .
    
    # Verify the binary
    ./diff-expr-analyzer --help
    """
    
    stub:
    """
    touch diff-expr-analyzer
    chmod +x diff-expr-analyzer
    """
}

process VALIDATE_INPUTS {
    tag "validating_inputs"
    
    input:
    path counts
    path metadata
    
    output:
    path counts, emit: validated_counts
    path metadata, emit: validated_metadata
    path "input_validation_report.txt", emit: report
    
    script:
    """
    #!/usr/bin/env python3
    
    import csv
    from collections import defaultdict
    
    print("Validating differential expression analysis inputs...")
    
    validation_report = []
    validation_report.append("Input Validation Report")
    validation_report.append("=======================")
    
    # Validate count data
    try:
        with open('${counts}', 'r') as f:
            reader = csv.reader(f, delimiter='\\t')
            header = next(reader)
            count_entries = sum(1 for _ in reader)
        
        validation_report.append(f"Count data: {count_entries} entries")
        validation_report.append(f"Count header: {header}")
        
    except Exception as e:
        validation_report.append(f"Error reading count data: {e}")
        raise
    
    # Validate metadata
    try:
        samples_in_metadata = set()
        groups_found = set()
        
        with open('${metadata}', 'r') as f:
            reader = csv.DictReader(f, delimiter='\\t')
            for row in reader:
                if 'sample_id' in row and 'group' in row:
                    samples_in_metadata.add(row['sample_id'])
                    groups_found.add(row['group'])
        
        validation_report.append(f"Metadata samples: {len(samples_in_metadata)}")
        validation_report.append(f"Groups found: {list(groups_found)}")
        
        # Check for required groups
        required_groups = {'${params.control_group}', '${params.treatment_group}'}
        missing_groups = required_groups - groups_found
        if missing_groups:
            validation_report.append(f"WARNING: Missing required groups: {missing_groups}")
        else:
            validation_report.append("✓ Required groups found")
            
    except Exception as e:
        validation_report.append(f"Error reading metadata: {e}")
        raise
    
    # Cross-validate samples
    try:
        samples_in_counts = set()
        with open('${counts}', 'r') as f:
            reader = csv.reader(f, delimiter='\\t')
            next(reader)  # Skip header
            for row in reader:
                if len(row) >= 2:
                    samples_in_counts.add(row[1])
        
        common_samples = samples_in_metadata & samples_in_counts
        validation_report.append(f"Common samples: {len(common_samples)}")
        
        if len(common_samples) < 4:
            validation_report.append("WARNING: Very few samples for differential analysis")
        else:
            validation_report.append("✓ Sufficient samples for analysis")
            
    except Exception as e:
        validation_report.append(f"Error cross-validating samples: {e}")
    
    # Write validation report
    with open('input_validation_report.txt', 'w') as f:
        f.write('\\n'.join(validation_report))
    
    print("Input validation completed successfully")
    """
    
    stub:
    """
    echo "Input validation passed" > input_validation_report.txt
    """
}

process DIFFERENTIAL_EXPRESSION {
    tag "differential_analysis"
    publishDir "${params.output_dir}", mode: 'copy'
    
    input:
    path analyzer_binary
    path counts
    path metadata
    
    output:
    path "differential_expression.tsv", emit: results
    path "de_analysis_stats.txt", emit: stats
    
    script:
    """
    # Make the binary executable
    chmod +x ${analyzer_binary}
    
    # Run differential expression analysis
    echo "Starting differential expression analysis..."
    echo "Input: ${counts}"
    echo "Metadata: ${metadata}"
    echo "Control group: ${params.control_group}"
    echo "Treatment group: ${params.treatment_group}"
    
    ./${analyzer_binary} \\
        --input ${counts} \\
        --metadata ${metadata} \\
        --output differential_expression.tsv \\
        --stats de_analysis_stats.txt \\
        --control ${params.control_group} \\
        --treatment ${params.treatment_group} \\
        --alpha ${params.alpha} \\
        --min-count ${params.min_count}
    
    # Verify outputs
    if [[ ! -s differential_expression.tsv ]]; then
        echo "ERROR: Differential expression analysis failed - no results generated"
        exit 1
    fi
    
    if [[ ! -s de_analysis_stats.txt ]]; then
        echo "ERROR: Statistics file not generated"
        exit 1
    fi
    
    echo "Differential expression analysis completed successfully"
    echo "Results entries: \$(tail -n +2 differential_expression.tsv | wc -l)"
    """
    
    stub:
    """
    echo -e "gene_id\\tcontrol_mean\\ttreatment_mean\\tlog2_fold_change\\tp_value\\tadjusted_p_value\\tsignificant" > differential_expression.tsv
    echo -e "GENE_00001\\t100.5\\t200.3\\t1.2\\t0.001\\t0.01\\ttrue" >> differential_expression.tsv
    echo "Differential Expression Analysis Statistics" > de_analysis_stats.txt
    echo "Total genes analyzed: 1000" >> de_analysis_stats.txt
    """
}

process GENERATE_PLOTS_DATA {
    tag "generating_plot_data"
    publishDir "${params.output_dir}", mode: 'copy'
    
    input:
    path de_results
    
    output:
    path "volcano_plot_data.tsv", emit: volcano
    path "ma_plot_data.tsv", emit: ma
    
    script:
    """
    #!/usr/bin/env python3
    
    import csv
    import math
    
    print("Generating plot data from differential expression results...")
    
    # Read DE results
    results = []
    with open('${de_results}', 'r') as f:
        reader = csv.DictReader(f, delimiter='\\t')
        for row in reader:
            try:
                result = {
                    'gene_id': row['gene_id'],
                    'control_mean': float(row['control_mean']),
                    'treatment_mean': float(row['treatment_mean']),
                    'log2_fold_change': float(row['log2_fold_change']),
                    'p_value': float(row['p_value']),
                    'adjusted_p_value': float(row['adjusted_p_value']),
                    'significant': row['significant'].lower() == 'true'
                }
                results.append(result)
            except (ValueError, KeyError) as e:
                print(f"Warning: Skipping malformed row for gene {row.get('gene_id', 'unknown')}: {e}")
    
    print(f"Processed {len(results)} genes for plotting")
    
    # Generate volcano plot data
    with open('volcano_plot_data.tsv', 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\\t')
        writer.writerow(['gene_id', 'log2_fold_change', 'neg_log10_p_value', 'significant', 'regulation'])
        
        for result in results:
            neg_log10_p = -math.log10(max(result['adjusted_p_value'], 1e-10))
            regulation = 'none'
            if result['significant']:
                if result['log2_fold_change'] > 0:
                    regulation = 'up'
                elif result['log2_fold_change'] < 0:
                    regulation = 'down'
            
            writer.writerow([
                result['gene_id'],
                f"{result['log2_fold_change']:.6f}",
                f"{neg_log10_p:.6f}",
                result['significant'],
                regulation
            ])
    
    # Generate MA plot data
    with open('ma_plot_data.tsv', 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\\t')
        writer.writerow(['gene_id', 'average_expression', 'log2_fold_change', 'significant', 'regulation'])
        
        for result in results:
            avg_expr = math.log2((result['control_mean'] + result['treatment_mean']) / 2 + 1)
            regulation = 'none'
            if result['significant']:
                if result['log2_fold_change'] > 0:
                    regulation = 'up'
                elif result['log2_fold_change'] < 0:
                    regulation = 'down'
            
            writer.writerow([
                result['gene_id'],
                f"{avg_expr:.6f}",
                f"{result['log2_fold_change']:.6f}",
                result['significant'],
                regulation
            ])
    
    print("Plot data generation completed successfully")
    """
    
    stub:
    """
    echo -e "gene_id\\tlog2_fold_change\\tneg_log10_p_value\\tsignificant\\tregulation" > volcano_plot_data.tsv
    echo -e "GENE_00001\\t1.2\\t2.5\\ttrue\\tup" >> volcano_plot_data.tsv
    echo -e "gene_id\\taverage_expression\\tlog2_fold_change\\tsignificant\\tregulation" > ma_plot_data.tsv
    echo -e "GENE_00001\\t7.5\\t1.2\\ttrue\\tup" >> ma_plot_data.tsv
    """
}

process GENERATE_SUMMARY {
    tag "generating_summary"
    publishDir "${params.output_dir}", mode: 'copy'
    
    input:
    path input_counts
    path metadata
    path de_results
    path de_stats
    path validation_report
    
    output:
    path "summary.txt"
    
    script:
    """
    #!/bin/bash
    
    # Count entries in files
    count_entries=0
    metadata_entries=0
    de_genes=0
    significant_genes=0
    
    if [ -f "${input_counts}" ]; then
        count_entries=\$(tail -n +2 "${input_counts}" | wc -l)
    fi
    
    if [ -f "${metadata}" ]; then
        metadata_entries=\$(tail -n +2 "${metadata}" | wc -l)
    fi
    
    if [ -f "${de_results}" ]; then
        de_genes=\$(tail -n +2 "${de_results}" | wc -l)
        significant_genes=\$(tail -n +2 "${de_results}" | awk -F'\\t' '\$7=="true" {count++} END {print count+0}')
    fi
    
    timestamp=\$(date '+%Y-%m-%d %H:%M:%S')
    significance_rate=\$(echo "scale=1; \${significant_genes} * 100 / \${de_genes}" | bc -l 2>/dev/null || echo "N/A")
    
    # Create comprehensive summary
    cat > summary.txt << SUMMARY_EOF
Differential Expression Analysis Pipeline Summary
================================================
Generated: \${timestamp}

PIPELINE STATUS: SUCCESS
========================
✅ Input Validation: COMPLETED
✅ Differential Analysis: COMPLETED
✅ Plot Data Generation: COMPLETED
✅ Statistical Reporting: COMPLETED

INPUT PARAMETERS:
================
- Normalized counts file: ${params.input}
- Sample metadata file: ${params.metadata}
- Control group: ${params.control_group}
- Treatment group: ${params.treatment_group}
- Significance threshold: ${params.alpha}
- Minimum count threshold: ${params.min_count}
- Output directory: ${params.output_dir}

DATA PROCESSING RESULTS:
=======================
- Input count entries: \${count_entries}
- Metadata samples: \${metadata_entries}
- Genes analyzed: \${de_genes}
- Significant genes: \${significant_genes}
- Significance rate: \${significance_rate}%

OUTPUT FILES:
============
✓ differential_expression.tsv (\${de_genes} genes analyzed)
✓ de_analysis_stats.txt (detailed statistics)
✓ volcano_plot_data.tsv (volcano plot data)
✓ ma_plot_data.tsv (MA plot data)
✓ input_validation_report.txt (validation details)
✓ summary.txt (this summary)

NEXT STEPS:
===========
1. Review de_analysis_stats.txt for detailed statistics
2. Use differential_expression.tsv for downstream analysis
3. Create volcano plot using volcano_plot_data.tsv
4. Create MA plot using ma_plot_data.tsv
5. Filter significant genes for pathway analysis

STATISTICAL SUMMARY:
===================
SUMMARY_EOF

    # Extract key statistics from de_analysis_stats.txt if available
    if [ -f "${de_stats}" ]; then
        echo "From detailed analysis:" >> summary.txt
        grep -E "(Total genes|Significant genes|Upregulated|Downregulated)" "${de_stats}" | sed 's/^/- /' >> summary.txt
    fi

    cat >> summary.txt << SUMMARY_EOF

VALIDATION REPORT:
=================
SUMMARY_EOF

    if [ -f "${validation_report}" ]; then
        cat "${validation_report}" >> summary.txt
    else
        echo "Validation report not available" >> summary.txt
    fi

    cat >> summary.txt << SUMMARY_EOF

PIPELINE COMPLETED SUCCESSFULLY
Data is ready for downstream analysis and visualization.

For visualization in R:
- Use volcano_plot_data.tsv for ggplot2 volcano plots
- Use ma_plot_data.tsv for MA plots
- Filter significant genes: awk -F'\\t' '\$7=="true"' differential_expression.tsv

For pathway analysis:
- Extract significant gene lists from differential_expression.tsv
- Use tools like GSEA, DAVID, or Enrichr for functional annotation
SUMMARY_EOF

    echo "Summary report created successfully"
    echo "Count entries: \${count_entries}"
    echo "DE genes: \${de_genes}"
    echo "Significant genes: \${significant_genes}"
    """
}

// Main workflow
workflow {
    // Create input channels
    input_counts = Channel.fromPath(params.input)
    metadata = Channel.fromPath(params.metadata)
    
    // Prepare source files for building
    cargo_toml = Channel.fromPath("Cargo.toml")
    main_rs = Channel.fromPath("src/main.rs")
    
    // Build the Rust analyzer
    analyzer_binary = BUILD_ANALYZER(cargo_toml, main_rs)
    
    // Validate inputs
    validation_results = VALIDATE_INPUTS(input_counts, metadata)
    
    // Run differential expression analysis
    de_results = DIFFERENTIAL_EXPRESSION(
        analyzer_binary,
        validation_results.validated_counts,
        validation_results.validated_metadata
    )
    
    // Generate plot data
    plot_data = GENERATE_PLOTS_DATA(de_results.results)
    
    // Generate comprehensive summary
    GENERATE_SUMMARY(
        validation_results.validated_counts,
        validation_results.validated_metadata,
        de_results.results,
        de_results.stats,
        validation_results.report
    )
}

// Completion handler
workflow.onComplete {
    log.info """
    🎯 Differential Expression Analysis Summary
    ==========================================
    Completed at: ${workflow.complete}
    Duration:     ${workflow.duration}
    Success:      ${workflow.success}
    Work dir:     ${workflow.workDir}
    Exit status:  ${workflow.exitStatus}
    """
    
    if (workflow.success) {
        log.info """
        🎉 SUCCESS! Differential expression analysis completed successfully.
        
        📁 Results available in: ${params.output_dir}/
        
        📊 Key outputs:
        ├── differential_expression.tsv    📈 Main DE results
        ├── de_analysis_stats.txt          📋 Detailed statistics
        ├── volcano_plot_data.tsv          🌋 Volcano plot data
        ├── ma_plot_data.tsv               📊 MA plot data
        ├── input_validation_report.txt    ✅ Validation details
        └── summary.txt                    📄 Comprehensive summary
        
        💡 Next steps:
        1. Review summary.txt for overview
        2. Check de_analysis_stats.txt for detailed statistics
        3. Use plot data files for visualization
        4. Filter significant genes for pathway analysis
        
        🔬 Your differential expression analysis is ready for biological interpretation!
        """
    } else {
        log.error """
        ❌ Pipeline failed. 
        Check error messages above and .nextflow.log for details.
        
        💡 Common fixes:
        - Ensure diff-expr-analyzer binary is built correctly
        - Check input file format and sample metadata
        - Verify group names match between metadata and parameters
        """
    }
}