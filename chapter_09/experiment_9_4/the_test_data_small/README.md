# Differential Expression Test Dataset: Small

## Dataset Information
- Number of genes: 100
- Samples per group: 3
- Total samples: 6
- Total count entries: 600

## Experimental Design
- **Control group**: 3 samples (Control_Rep1 to Control_Rep3)
- **Treatment group**: 3 samples (Treatment_Rep1 to Treatment_Rep3)

## Simulated Differential Expression
- **Upregulated genes**: ~10% (2-8 fold increase in treatment)
- **Downregulated genes**: ~10% (2-8 fold decrease in treatment)
- **Unchanged genes**: ~80% (minor random variation)

## Files
- `normalized_counts.tsv`: Normalized count data (gene_id, sample_id, count)
- `sample_metadata.tsv`: Sample metadata with group assignments

## Expected Results
When analyzing this dataset, you should expect:
- Significant genes: ~20% of total genes
- False discovery rate: <5% with proper multiple testing correction
- Clear separation between control and treatment groups in PCA

## Usage

### Basic Analysis
```bash
# Build the differential expression analyzer
cargo build --release
cp target/release/diff-expr-analyzer .

# Run analysis
./diff-expr-analyzer \
    --input normalized_counts.tsv \
    --metadata sample_metadata.tsv \
    --output de_results.tsv \
    --stats de_stats.txt \
    --control Control \
    --treatment Treatment
```

### Nextflow Pipeline
```bash
# Run full pipeline
nextflow run main.nf \
    --input normalized_counts.tsv \
    --metadata sample_metadata.tsv \
    --control_group Control \
    --treatment_group Treatment \
    --output_dir results_small
```

## Quality Control Expectations
- Input validation should pass
- Statistical power should be adequate for 3 samples per group
- Multiple testing correction should be applied
- Volcano and MA plots should show clear differential expression patterns

## Troubleshooting
If analysis fails:
1. Check that sample IDs in metadata match those in count data
2. Verify group names are exactly 'Control' and 'Treatment'
3. Ensure sufficient samples per group (minimum 3 recommended)
4. Check for missing or invalid count values
