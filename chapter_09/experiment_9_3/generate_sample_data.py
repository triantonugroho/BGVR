#!/usr/bin/env python3
"""
Sample RNA-seq dataset generator for testing the normalization pipeline
"""

import pandas as pd
import numpy as np
import argparse
from pathlib import Path

def generate_count_data(n_genes=1000, n_samples=12, output_file="raw_counts.tsv"):
    """
    Generate synthetic RNA-seq count data
    """
    np.random.seed(42)  # For reproducibility
    
    # Generate gene IDs
    genes = [f"GENE_{i:05d}" for i in range(1, n_genes + 1)]
    
    # Generate sample IDs with different conditions and batches
    samples = []
    conditions = ['Control', 'Treatment']
    batches = ['Batch1', 'Batch2', 'Batch3']
    
    for i in range(n_samples):
        condition = conditions[i % len(conditions)]
        batch = batches[i % len(batches)]
        samples.append(f"{condition}_{batch}_Rep{i//len(conditions) + 1}")
    
    print(f"Generating count data for {len(genes)} genes and {len(samples)} samples")
    print(f"Samples: {samples}")
    
    # Generate count data
    count_data = []
    
    for gene in genes:
        # Base expression level (log-normal distribution)
        base_expression = np.random.lognormal(mean=5, sigma=2)
        
        for sample in samples:
            # Add sample-specific and batch effects
            batch_effect = 1.0
            if 'Batch2' in sample:
                batch_effect = 1.5  # 50% higher expression in Batch2
            elif 'Batch3' in sample:
                batch_effect = 0.7  # 30% lower expression in Batch3
            
            # Treatment effect for some genes
            treatment_effect = 1.0
            if 'Treatment' in sample and gene.endswith(('1', '2', '3', '4', '5')):
                # 50% of genes show treatment effect
                if int(gene.split('_')[1]) % 2 == 0:
                    treatment_effect = 2.0  # Upregulated
                else:
                    treatment_effect = 0.5  # Downregulated
            
            # Final count with noise
            expected_count = base_expression * batch_effect * treatment_effect
            # Add Poisson noise
            actual_count = np.random.poisson(expected_count)
            
            count_data.append({
                'gene_id': gene,
                'sample_id': sample,
                'count': actual_count
            })
    
    # Convert to DataFrame and save
    df = pd.DataFrame(count_data)
    df.to_csv(output_file, sep='\t', index=False, header=False)
    
    print(f"Count data saved to {output_file}")
    print(f"Total entries: {len(count_data)}")
    
    # Generate summary statistics
    summary_stats = df.groupby('sample_id')['count'].agg(['sum', 'mean', 'std']).round(2)
    print("\nSample summary statistics:")
    print(summary_stats)
    
    return df, samples

def generate_batch_metadata(samples, output_file="batch_metadata.tsv"):
    """
    Generate batch metadata file
    """
    metadata = []
    
    for sample in samples:
        # Extract condition and batch from sample name
        parts = sample.split('_')
        condition = parts[0]
        batch = parts[1]
        
        metadata.append({
            'sample_id': sample,
            'condition': condition,
            'batch': batch,
            'replicate': parts[2] if len(parts) > 2 else 'Rep1'
        })
    
    df_meta = pd.DataFrame(metadata)
    df_meta.to_csv(output_file, sep='\t', index=False)
    
    print(f"\nBatch metadata saved to {output_file}")
    print(df_meta)
    
    return df_meta

def create_test_datasets():
    """
    Create multiple test datasets of different sizes
    """
    datasets = [
        {"name": "small", "genes": 100, "samples": 6},
        {"name": "medium", "genes": 1000, "samples": 12},
        {"name": "large", "genes": 5000, "samples": 24}
    ]
    
    for dataset in datasets:
        print(f"\n{'='*50}")
        print(f"Creating {dataset['name']} dataset")
        print(f"{'='*50}")
        
        output_dir = Path(f"test_data_{dataset['name']}")
        output_dir.mkdir(exist_ok=True)
        
        count_file = output_dir / "raw_counts.tsv"
        metadata_file = output_dir / "batch_metadata.tsv"
        
        df, samples = generate_count_data(
            n_genes=dataset['genes'],
            n_samples=dataset['samples'],
            output_file=str(count_file)
        )
        
        generate_batch_metadata(samples, str(metadata_file))
        
        # Create a README for this dataset
        readme_content = f"""# Test Dataset: {dataset['name'].title()}

## Dataset Information
- Number of genes: {dataset['genes']:,}
- Number of samples: {dataset['samples']}
- Total count entries: {dataset['genes'] * dataset['samples']:,}

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
"""
        
        with open(output_dir / "README.md", 'w') as f:
            f.write(readme_content)

def main():
    parser = argparse.ArgumentParser(description='Generate sample RNA-seq datasets')
    parser.add_argument('--genes', type=int, default=1000, help='Number of genes')
    parser.add_argument('--samples', type=int, default=12, help='Number of samples')
    parser.add_argument('--output', default='raw_counts.tsv', help='Output file name')
    parser.add_argument('--create-test-sets', action='store_true', 
                       help='Create multiple test datasets')
    
    args = parser.parse_args()
    
    if args.create_test_sets:
        create_test_datasets()
    else:
        df, samples = generate_count_data(args.genes, args.samples, args.output)
        generate_batch_metadata(samples, "batch_metadata.tsv")

if __name__ == "__main__":
    main()