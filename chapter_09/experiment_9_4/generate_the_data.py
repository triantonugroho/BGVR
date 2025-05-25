#!/usr/bin/env python3
"""
Sample data generator for differential expression analysis testing
"""

import pandas as pd
import numpy as np
import argparse
from pathlib import Path

def generate_normalized_counts(n_genes=1000, n_samples_per_group=6, output_file="normalized_counts.tsv"):
    """
    Generate synthetic normalized RNA-seq count data with differential expression
    """
    np.random.seed(42)  # For reproducibility
    
    # Generate gene IDs
    genes = [f"GENE_{i:05d}" for i in range(1, n_genes + 1)]
    
    # Generate sample IDs for two groups
    control_samples = [f"Control_Rep{i+1}" for i in range(n_samples_per_group)]
    treatment_samples = [f"Treatment_Rep{i+1}" for i in range(n_samples_per_group)]
    all_samples = control_samples + treatment_samples
    
    print(f"Generating normalized counts for {len(genes)} genes and {len(all_samples)} samples")
    print(f"Control samples: {control_samples}")
    print(f"Treatment samples: {treatment_samples}")
    
    # Generate count data with differential expression
    count_data = []
    
    # Define different gene categories
    n_upregulated = int(0.1 * n_genes)      # 10% upregulated
    n_downregulated = int(0.1 * n_genes)    # 10% downregulated  
    n_unchanged = n_genes - n_upregulated - n_downregulated  # 80% unchanged
    
    gene_categories = (['upregulated'] * n_upregulated + 
                      ['downregulated'] * n_downregulated + 
                      ['unchanged'] * n_unchanged)
    
    for i, gene in enumerate(genes):
        category = gene_categories[i]
        
        # Base expression level (log-normal distribution)
        base_expression = np.random.lognormal(mean=4, sigma=1.5)  # Moderate expression
        
        # Generate control samples
        for sample in control_samples:
            # Add biological and technical noise
            noise_factor = np.random.lognormal(mean=0, sigma=0.3)
            control_count = base_expression * noise_factor
            
            count_data.append({
                'gene_id': gene,
                'sample_id': sample,
                'count': max(0.1, control_count)  # Minimum count of 0.1
            })
        
        # Generate treatment samples with differential expression
        for sample in treatment_samples:
            # Apply differential expression based on category
            if category == 'upregulated':
                # 2-8 fold upregulation
                fold_change = np.random.uniform(2.0, 8.0)
                treatment_expression = base_expression * fold_change
            elif category == 'downregulated':
                # 2-8 fold downregulation
                fold_change = np.random.uniform(0.125, 0.5)  # 1/8 to 1/2
                treatment_expression = base_expression * fold_change
            else:  # unchanged
                # Minor random variation around base expression
                fold_change = np.random.uniform(0.8, 1.25)
                treatment_expression = base_expression * fold_change
            
            # Add biological and technical noise
            noise_factor = np.random.lognormal(mean=0, sigma=0.3)
            treatment_count = treatment_expression * noise_factor
            
            count_data.append({
                'gene_id': gene,
                'sample_id': sample,
                'count': max(0.1, treatment_count)  # Minimum count of 0.1
            })
    
    # Convert to DataFrame and save
    df = pd.DataFrame(count_data)
    df.to_csv(output_file, sep='\t', index=False)
    
    print(f"Normalized count data saved to {output_file}")
    print(f"Total entries: {len(count_data)}")
    
    # Generate summary statistics
    summary_stats = df.groupby('sample_id')['count'].agg(['count', 'mean', 'std']).round(2)
    print("\nSample summary statistics:")
    print(summary_stats)
    
    # Print gene category information
    print(f"\nGene categories:")
    print(f"- Upregulated: {n_upregulated} genes ({n_upregulated/n_genes*100:.1f}%)")
    print(f"- Downregulated: {n_downregulated} genes ({n_downregulated/n_genes*100:.1f}%)")
    print(f"- Unchanged: {n_unchanged} genes ({n_unchanged/n_genes*100:.1f}%)")
    
    return df, all_samples, gene_categories

def generate_sample_metadata(samples, output_file="sample_metadata.tsv"):
    """
    Generate sample metadata file with group assignments
    """
    metadata = []
    
    for sample in samples:
        # Determine group based on sample name
        if sample.startswith('Control'):
            group = 'Control'
            condition = 'Control'
        elif sample.startswith('Treatment'):
            group = 'Treatment'
            condition = 'Treatment'
        else:
            group = 'Unknown'
            condition = 'Unknown'
        
        # Extract replicate number
        try:
            replicate = sample.split('_Rep')[1] if '_Rep' in sample else '1'
        except:
            replicate = '1'
        
        metadata.append({
            'sample_id': sample,
            'group': group,
            'condition': condition,
            'replicate': replicate,
            'batch': f"Batch{(int(replicate) - 1) % 3 + 1}"  # Assign to 3 batches cyclically
        })
    
    df_meta = pd.DataFrame(metadata)
    df_meta.to_csv(output_file, sep='\t', index=False)
    
    print(f"\nSample metadata saved to {output_file}")
    print(df_meta)
    
    return df_meta

def create_test_datasets():
    """
    Create multiple test datasets of different sizes for DE analysis
    """
    datasets = [
        {"name": "small", "genes": 100, "samples_per_group": 3},
        {"name": "medium", "genes": 1000, "samples_per_group": 6},
        {"name": "large", "genes": 5000, "samples_per_group": 12}
    ]
    
    for dataset in datasets:
        print(f"\n{'='*60}")
        print(f"Creating {dataset['name']} differential expression dataset")
        print(f"{'='*60}")
        
        output_dir = Path(f"the_test_data_{dataset['name']}")
        output_dir.mkdir(exist_ok=True)
        
        counts_file = output_dir / "normalized_counts.tsv"
        metadata_file = output_dir / "sample_metadata.tsv"
        
        df, samples, gene_categories = generate_normalized_counts(
            n_genes=dataset['genes'],
            n_samples_per_group=dataset['samples_per_group'],
            output_file=str(counts_file)
        )
        
        generate_sample_metadata(samples, str(metadata_file))
        
        # Create a detailed README for this dataset
        readme_content = f"""# Differential Expression Test Dataset: {dataset['name'].title()}

## Dataset Information
- Number of genes: {dataset['genes']:,}
- Samples per group: {dataset['samples_per_group']}
- Total samples: {dataset['samples_per_group'] * 2}
- Total count entries: {dataset['genes'] * dataset['samples_per_group'] * 2:,}

## Experimental Design
- **Control group**: {dataset['samples_per_group']} samples (Control_Rep1 to Control_Rep{dataset['samples_per_group']})
- **Treatment group**: {dataset['samples_per_group']} samples (Treatment_Rep1 to Treatment_Rep{dataset['samples_per_group']})

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
./diff-expr-analyzer \\
    --input normalized_counts.tsv \\
    --metadata sample_metadata.tsv \\
    --output the_results.tsv \\
    --stats the_stats.txt \\
    --control Control \\
    --treatment Treatment
```

### Nextflow Pipeline
```bash
# Run full pipeline
nextflow run main.nf \\
    --input normalized_counts.tsv \\
    --metadata sample_metadata.tsv \\
    --control_group Control \\
    --treatment_group Treatment \\
    --output_dir results_{dataset['name']}
```

## Quality Control Expectations
- Input validation should pass
- Statistical power should be adequate for {dataset['samples_per_group']} samples per group
- Multiple testing correction should be applied
- Volcano and MA plots should show clear differential expression patterns

## Troubleshooting
If analysis fails:
1. Check that sample IDs in metadata match those in count data
2. Verify group names are exactly 'Control' and 'Treatment'
3. Ensure sufficient samples per group (minimum 3 recommended)
4. Check for missing or invalid count values
"""
        
        with open(output_dir / "README.md", 'w') as f:
            f.write(readme_content)
        
        # Create expected results summary
        expected_results = f"""Expected Differential Expression Results
========================================

Dataset: {dataset['name'].title()}
Genes: {dataset['genes']:,}
Samples per group: {dataset['samples_per_group']}

Expected significant genes: ~{int(0.2 * dataset['genes'])} ({0.2 * 100:.0f}%)
Expected upregulated: ~{int(0.1 * dataset['genes'])} ({0.1 * 100:.0f}%)
Expected downregulated: ~{int(0.1 * dataset['genes'])} ({0.1 * 100:.0f}%)

Statistical Power:
- With {dataset['samples_per_group']} samples per group
- Expected to detect 2+ fold changes
- At 5% FDR with 80%+ power

Visualization Expectations:
- Volcano plot: Clear separation of significant genes
- MA plot: Even distribution around zero for non-DE genes
- P-value histogram: Enrichment of small p-values
"""
        
        with open(output_dir / "expected_results.txt", 'w') as f:
            f.write(expected_results)
        
        print(f"Dataset created in: {output_dir}")
        print(f"Files: {counts_file.name}, {metadata_file.name}, README.md, expected_results.txt")

def generate_realistic_counts_from_existing(existing_file, output_file="realistic_normalized_counts.tsv"):
    """
    Convert existing count data to more realistic differential expression format
    """
    print(f"Converting existing count data from {existing_file}")
    
    # Read existing data
    df = pd.read_csv(existing_file, sep='\t')
    print(f"Read {len(df)} entries from existing file")
    
    # Extract unique genes and samples
    genes = df['gene_id'].unique()
    samples = df['sample_id'].unique()
    
    print(f"Found {len(genes)} genes and {len(samples)} samples")
    
    # Create realistic metadata based on sample names
    metadata = []
    for sample in samples:
        if 'Control' in sample or 'control' in sample:
            group = 'Control'
        elif 'Treatment' in sample or 'treatment' in sample or 'Treat' in sample:
            group = 'Treatment'
        else:
            # Assign alternating groups if unclear
            group = 'Control' if len(metadata) % 2 == 0 else 'Treatment'
        
        metadata.append({
            'sample_id': sample,
            'group': group,
            'condition': group,
            'replicate': str(len([m for m in metadata if m['group'] == group]) + 1)
        })
    
    # Save realistic metadata
    pd.DataFrame(metadata).to_csv("realistic_sample_metadata.tsv", sep='\t', index=False)
    
    # Use existing counts as they are (already normalized)
    df.to_csv(output_file, sep='\t', index=False)
    
    print(f"Realistic count data saved to {output_file}")
    print(f"Metadata saved to realistic_sample_metadata.tsv")
    
    return df, metadata

def main():
    parser = argparse.ArgumentParser(description='Generate sample data for differential expression analysis')
    parser.add_argument('--genes', type=int, default=1000, help='Number of genes')
    parser.add_argument('--samples-per-group', type=int, default=6, help='Number of samples per group')
    parser.add_argument('--output', default='normalized_counts.tsv', help='Output count file name')
    parser.add_argument('--create-test-sets', action='store_true', 
                       help='Create multiple test datasets')
    parser.add_argument('--convert-existing', type=str, 
                       help='Convert existing count file to DE format')
    
    args = parser.parse_args()
    
    if args.convert_existing:
        generate_realistic_counts_from_existing(args.convert_existing)
    elif args.create_test_sets:
        create_test_datasets()
    else:
        df, samples, _ = generate_normalized_counts(args.genes, args.samples_per_group, args.output)
        generate_sample_metadata(samples, "sample_metadata.tsv")

if __name__ == "__main__":
    main()