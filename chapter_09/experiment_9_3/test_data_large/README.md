# Test Dataset: Large

## Dataset Information
- Number of genes: 5,000
- Number of samples: 24
- Total count entries: 120,000

## Files
- `raw_counts.tsv`: Raw count data (gene_id, sample_id, count)
- `batch_metadata.tsv`: Sample metadata with batch information

## Sample Design
- Conditions: Control vs Treatment
- Batches: Batch1, Batch2, Batch3
- Replicates: Multiple replicates per condition/batch combination

## Usage
```bash
# Run normalization
./rnaseq-normalizer --input raw_counts.tsv --output normalized_counts.tsv

# Run full pipeline with batch correction
nextflow run main.nf --input raw_counts.tsv --batch_correct true --batch_metadata batch_metadata.tsv
```
