import os
import json
import pickle
import pandas as pd
import glob
import sys
from tqdm import tqdm

def process_single_pkl(pkl_file):
    if pkl_file.endswith('.json'):
        with open(pkl_file, 'r') as f:
            data = json.load(f)
        return {
            'ptm': data.get('ptm'),
            'iptm': data.get('iptm'), 
            'ranking_confidence': data.get('ranking_confidence')
        }
    elif pkl_file.endswith('.pkl'):
        with open(pkl_file, 'rb') as f:
            data = pickle.load(f)
        return {
            'ptm': data.get('ptm'),
            'iptm': data.get('iptm'), 
            'ranking_confidence': data.get('ranking_confidence')
        }
    else:
        raise Exception(f'{pkl_file} is not a .pkl or .json file')
    return None

def process_pdb(pdb_id):
    # Convert PDB ID to uppercase
    pdb_id = pdb_id.upper()
    
    # Define the directory path
    pkl_dir = f'./PDB_FILES/JSON/{pdb_id}_json/'
    output_dir = f'./'
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Find all .pkl and .json files
    pkl_files = glob.glob(os.path.join(pkl_dir, '*.pkl'))
    json_files = glob.glob(os.path.join(pkl_dir, '*.json'))
    pkl_files = pkl_files + json_files

    if not pkl_files:
        print(f"No .pkl files found in {pkl_dir}")
        return
    
    # Initialize lists to store results
    results = []
    
    # Process each .pkl file
    for pkl_file in tqdm(pkl_files):
        # Get the basename of the file
        basename = os.path.basename(pkl_file).split('.')[0]
        if basename.startswith('result_') == False:
            continue
        
        # Process the PKL file
        metrics = process_single_pkl(pkl_file)
        if metrics and isinstance(metrics, dict):
            # Add model name to metrics
            metrics['model'] = basename.replace('result_', '') + '_af2'
            results.append(metrics)
        elif metrics and isinstance(metrics, list):
            print(metrics)
            raise Exception(f'metrics is a list for {pkl_file}')

    
    if results:
        # Create DataFrame and save to CSV
        df = pd.DataFrame(results)
        output_file = os.path.join(output_dir, f'{pdb_id}_ranking_confidences.csv')
        df.to_csv(output_file, index=False)
        print(f"Results saved to {output_file}")
    else:
        print("No results were generated")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python process_pdb_metrics.py <pdb_id>")
        sys.exit(1)
    
    pdb_id = sys.argv[1]
    process_pdb(pdb_id.upper()) 