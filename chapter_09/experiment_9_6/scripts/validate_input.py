#!/usr/bin/env python3

import pandas as pd
import sys
import os

def main():
    if len(sys.argv) != 3:
        print("Usage: python validate_input.py <input_file> <output_dir>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_dir = sys.argv[2]
    
    try:
        print(f"Validating input matrix: {input_file}")
        df = pd.read_csv(input_file, sep='\t', header=0)
        print(f"Input shape: {df.shape}")
        
        # Check required columns
        required_cols = ['gene_idx', 'cell_idx', 'count']
        missing_cols = [col for col in required_cols if col not in df.columns]
        
        if missing_cols:
            raise ValueError(f"Missing required columns: {missing_cols}")
        
        # Remove zero counts
        df_filtered = df[df['count'] > 0].copy()
        
        # Save validated matrix
        output_file = os.path.join(output_dir, "validated_matrix.tsv")
        df_filtered.to_csv(output_file, sep='\t', index=False)
        
        # Generate report
        report = f"""Input Validation Report
======================
Original entries: {len(df):,}
Non-zero entries: {len(df_filtered):,}
Unique genes: {df_filtered['gene_idx'].nunique():,}
Unique cells: {df_filtered['cell_idx'].nunique():,}
Total counts: {df_filtered['count'].sum():,.0f}"""
        
        report_file = os.path.join(output_dir, "validation_report.txt")
        with open(report_file, 'w') as f:
            f.write(report)
            
        print("Validation completed successfully")
        print(f"Output files: {output_file}, {report_file}")
        
    except Exception as e:
        print(f"Validation failed: {e}")
        # Make sure output directory exists
        os.makedirs(output_dir, exist_ok=True)
        report_file = os.path.join(output_dir, "validation_report.txt")
        with open(report_file, 'w') as f:
            f.write(f"ERROR: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
