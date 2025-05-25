#!/usr/bin/env python3
"""
Fixed version of synthetic single-cell data generator
Handles numpy data type serialization properly
"""

import numpy as np
import pandas as pd
import json
from pathlib import Path

def convert_numpy_types(obj):
    """Convert numpy types to native Python types for JSON serialization."""
    if isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, dict):
        return {key: convert_numpy_types(value) for key, value in obj.items()}
    elif isinstance(obj, list):
        return [convert_numpy_types(item) for item in obj]
    else:
        return obj

def generate_simple_synthetic_data(n_cells=1000, n_genes=2000, seed=42):
    """Generate simple synthetic single-cell data."""
    np.random.seed(seed)
    
    print(f"Generating synthetic data:")
    print(f"  Cells: {n_cells:,}")
    print(f"  Genes: {n_genes:,}")
    
    # Generate realistic sparse entries
    target_entries = int(n_cells * n_genes * 0.15)  # 15% fill rate
    
    entries = []
    for _ in range(target_entries):
        gene_idx = np.random.randint(0, n_genes)
        cell_idx = np.random.randint(0, n_cells)
        
        # Generate realistic count using exponential + poisson
        base_rate = np.random.exponential(1.5)
        count = max(1, int(np.random.poisson(base_rate) + 1))
        
        entries.append({
            'gene_idx': int(gene_idx),  # Ensure native Python int
            'cell_idx': int(cell_idx),  # Ensure native Python int
            'count': float(count)       # Ensure native Python float
        })
    
    # Convert to DataFrame and remove duplicates
    df = pd.DataFrame(entries)
    df = df.groupby(['gene_idx', 'cell_idx'])['count'].max().reset_index()
    
    # Add cell type structure
    n_types = 5
    enhanced_entries = []
    
    for _, row in df.iterrows():
        gene_idx = int(row['gene_idx'])
        cell_idx = int(row['cell_idx'])
        count = float(row['count'])
        
        # Assign cell type
        cell_type = cell_idx % n_types
        
        # Boost signature genes for each cell type
        type_signature_size = n_genes // n_types
        signature_start = cell_type * type_signature_size
        signature_end = min((cell_type + 1) * type_signature_size, n_genes)
        
        if signature_start <= gene_idx < signature_end:
            count *= np.random.uniform(2.0, 4.0)
        
        enhanced_entries.append({
            'gene_idx': gene_idx,
            'cell_idx': cell_idx,
            'count': count
        })
    
    df_final = pd.DataFrame(enhanced_entries)
    df_final = df_final.sort_values(['gene_idx', 'cell_idx'])
    
    return df_final

def main():
    print("=== Fixed Single-Cell Data Generator ===")
    
    # Create directories
    Path("data/raw").mkdir(parents=True, exist_ok=True)
    Path("data/synthetic").mkdir(parents=True, exist_ok=True)
    
    # Generate data
    df = generate_simple_synthetic_data()
    
    # Save sparse matrix
    output_file = "data/raw/sparse_counts.tsv"
    df.to_csv(output_file, sep='\t', index=False)
    
    # Calculate statistics (ensure all are native Python types)
    stats = {
        'total_entries': int(len(df)),
        'unique_genes': int(df['gene_idx'].nunique()),
        'unique_cells': int(df['cell_idx'].nunique()),
        'total_counts': float(df['count'].sum()),
        'mean_count': float(df['count'].mean()),
        'median_count': float(df['count'].median()),
        'max_count': float(df['count'].max()),
        'min_count': float(df['count'].min()),
        'generation_parameters': {
            'n_cells': 1000,
            'n_genes': 2000,
            'n_cell_types': 5,
            'seed': 42
        }
    }
    
    # Calculate sparsity
    matrix_size = stats['unique_genes'] * stats['unique_cells']
    stats['sparsity'] = float(1 - (stats['total_entries'] / matrix_size))
    
    # Save metadata (now it should work)
    metadata_file = "data/synthetic/metadata.json"
    try:
        with open(metadata_file, 'w') as f:
            json.dump(stats, f, indent=2)
        print(f"✅ Metadata saved to: {metadata_file}")
    except Exception as e:
        print(f"⚠️  Could not save metadata: {e}")
        print("Continuing without metadata file...")
    
    print(f"\n✅ Data generation completed!")
    print(f"   Output: {output_file}")
    print(f"   Entries: {stats['total_entries']:,}")
    print(f"   Genes: {stats['unique_genes']:,}")
    print(f"   Cells: {stats['unique_cells']:,}")
    print(f"   Total counts: {stats['total_counts']:,.0f}")
    print(f"   Mean count: {stats['mean_count']:.2f}")
    print(f"   Sparsity: {stats['sparsity']:.3f}")
    
    # Show sample data
    print(f"\nSample data preview:")
    print(df.head(10).to_string(index=False))
    
    return True

if __name__ == "__main__":
    main()
