#!/usr/bin/env python3
"""
Script to compute TM-scores, DockQ, ICS, and IPS using USalign and other tools between a reference structure and all processed PDB files.
Adds the computed scores as new columns to the ranking confidences CSV file.
"""

import os
import sys
import subprocess
import pandas as pd
from pathlib import Path
import time

# Import functions from util.py
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from util import compute_DockQ_score, calculate_ics_ips, execute_command_and_get_output

def compute_structure_similarity(pdb_file1, pdb_file2, multimer=True):
    """
    Computes structural similarity between two PDB files using USalign.
    Optimized for protein-peptide complexes with multimeric support.
    
    Args:
        pdb_file1 (str): Path to first PDB file
        pdb_file2 (str): Path to second PDB file  
        multimer (bool): Whether to treat structures as multimeric complexes
        
    Returns:
        float: Similarity score (TM-score)
    """
    try:
        # Build USalign command with appropriate options
        # -mol prot: only align proteins
        # -mm 1: alignment of multi-chain oligomeric structures
        # -ter 1: align all chains of the first model (for asymmetric units)
        # -het 1: align both ATOM and HETATM residues (important for peptides)
        base_cmd = 'USalign'
        options = ['-mol prot']
        
        if multimer:
            options.extend(['-mm 1', '-ter 1'])
        options.append('-het 1')
        
        clistring = f'{base_cmd} {pdb_file1} {pdb_file2} {" ".join(options)}'
        print(f"Running: {clistring}")
        result = subprocess.run(clistring, shell=True, capture_output=True, text=True)
        
        # Parse TM-score from output
        tm_score = None
        for line in result.stdout.splitlines():
            if "TM-score=" in line:
                # USalign outputs multiple TM-scores, we want the first one
                # which is normalized by the length of the first structure
                tm_score = float(line.split()[1])
                break
                
        if tm_score is None:
            print(f"Warning: Failed to parse TM-score from USalign output for {pdb_file2}")
            print("USalign output:")
            print(result.stdout)
            print("USalign error:")
            print(result.stderr)
            return None
            
        return tm_score
        
    except Exception as e:
        print(f"Error computing structure similarity for {pdb_file2}: {e}")
        return None

def main():
    # Define paths
    reference_pdb = "./Files_from_Alberto_7_17/Seq951-AF3_relaxed.pdb"
    csv_file = "./LANA_ET/LANA_ET_SCORES.csv"
    pdb_dir = "./PDB_FILES/LANA_ET_processed/"
    
    # Check if reference PDB exists
    if not os.path.exists(reference_pdb):
        print(f"Error: Reference PDB file not found: {reference_pdb}")
        sys.exit(1)
    
    # Check if CSV file exists
    if not os.path.exists(csv_file):
        print(f"Error: CSV file not found: {csv_file}")
        sys.exit(1)
    
    # Check if PDB directory exists
    if not os.path.exists(pdb_dir):
        print(f"Error: PDB directory not found: {pdb_dir}")
        sys.exit(1)
    
    # Read the CSV file
    print(f"Reading CSV file: {csv_file}")
    df = pd.read_csv(csv_file)
    
    # Check if 'pdb_file' column exists
    if 'pdb_file' not in df.columns:
        print("Error: 'pdb_file' column not found in CSV file")
        sys.exit(1)
    
    # Define new columns
    new_columns = ['tm_score_to_MELD', 'dockq_to_MELD']
    
    # Check if columns already exist, if not add them
    for col in new_columns:
        if col not in df.columns:
            df[col] = None
            print(f"Added new column: {col}")
        else:
            print(f"Column {col} already exists - will preserve existing values")
    
    # Determine output file name
    output_file = csv_file.replace('.csv', '_with_structural_scores.csv')
    print(f"Results will be saved to: {output_file}")
    
    # Get total number of models to process
    total_models = len(df)
    print(f"Processing {total_models} models...")
    
    # Process each model
    for idx, row in df.iterrows():
        pdb_filename = row['pdb_file']
        pdb_file = os.path.join(pdb_dir, pdb_filename)
        
        print(f"Processing {idx+1}/{total_models}: {pdb_filename}")
        
        # Check if all metrics are already calculated for this model
        all_metrics_exist = all(pd.notna(row.get(col, None)) for col in new_columns)
        if all_metrics_exist:
            print(f"  All metrics already calculated - skipping")
            continue
        
        # Check if PDB file exists
        if not os.path.exists(pdb_file):
            print(f"Warning: PDB file not found: {pdb_file}")
            # Set all metrics to None for this model
            for col in new_columns:
                df.at[idx, col] = None
            continue
        
        # Compute TM-score
        tm_score = compute_structure_similarity(reference_pdb, pdb_file, multimer=True)
        df.at[idx, 'tm_score_to_MELD'] = tm_score
        
        # Compute DockQ score
        dockq_score = None
        try:
            iRMS, LRMS, DockQ, Fnat, Fnonnat, F1, clashes = compute_DockQ_score(reference_pdb, pdb_file)
            dockq_score = DockQ
        except Exception as e:
            print(f"  DockQ computation failed: {e}")
        df.at[idx, 'dockq_to_MELD'] = dockq_score
        

        
        # Print progress
        print(f"  TM-score: {tm_score:.4f}" if tm_score is not None else "  TM-score: Failed to compute")
        print(f"  DockQ: {dockq_score:.4f}" if dockq_score is not None else "  DockQ: Failed to compute")
        
        # Save progress after each model (incremental saving)
        df.to_csv(output_file, index=False)
        print(f"  Progress saved to {output_file}")
        
        # Add a small delay to avoid overwhelming the system
        # time.sleep(0.1)
    
    print(f"\nAll processing completed. Final results saved to: {output_file}")
    
    # Print summary
    print(f"\nSummary:")
    print(f"Total models processed: {total_models}")
    
    # Statistics for each metric
    metrics = ['tm_score_to_MELD', 'dockq_to_MELD']
    for metric in metrics:
        successful_computations = df[metric].notna().sum()
        metric_name = metric.replace('_to_MELD', '').upper()
        print(f"\n{metric_name} statistics:")
        print(f"  Successful computations: {successful_computations}")
        print(f"  Failed computations: {total_models - successful_computations}")
        
        if successful_computations > 0:
            scores = df[metric].dropna()
            print(f"  Mean: {scores.mean():.4f}")
            print(f"  Median: {scores.median():.4f}")
            print(f"  Min: {scores.min():.4f}")
            print(f"  Max: {scores.max():.4f}")
            print(f"  Std: {scores.std():.4f}")

if __name__ == "__main__":
    main()
