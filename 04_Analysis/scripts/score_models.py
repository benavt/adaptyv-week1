import os
import pandas as pd
import sys
import argparse
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Sequence import Sequence, SubSequence
from util import write_to_csv, find_occluded_nh_groups_cleanex, find_buried_residues_sasa

def parse_arguments():
    parser = argparse.ArgumentParser(description='Score models for a given PDB ID')
    parser.add_argument('pdb_id', type=str, help='PDB ID (e.g., LANA_ET, BICRAV7R_ET)')
    parser.add_argument('--pdb_file', type=str, help='Optional: specific PDB file to score (e.g., model_1_multimer__105_af2.pdb)')
    return parser.parse_args()

def get_data_seq(df_exp, df_bool):
    sequence = ['_']*100
    residue_numbers = []
    values = {'CSP': [], 'Intensity Ratio': [], 'Significant CSPs': [], 
              'Significant Intensity Ratios': [], 'CSP_bool': [], 'Intensity Ratio_bool': [], 
              'CleanEx_Holo_bool': [], 'CleanEx_Apo_bool': []}
    for i in range(len(df_exp)):
        sequence[int(df_exp.iloc[i]['Sequence'][1:])] = df_exp.iloc[i]['Sequence'][0]
        residue_numbers.append(int(df_exp.iloc[i]['Sequence'][1:]))
        values['CSP'].append(df_exp.iloc[i]['CSP'])
        values['Intensity Ratio'].append(df_exp.iloc[i]['Intensity Ratio'])
        values['Significant CSPs'].append(df_exp.iloc[i]['Significant CSPs'])
        values['Significant Intensity Ratios'].append(df_exp.iloc[i]['Significant Intensity Ratios'])
        values['CSP_bool'].append(df_bool.iloc[i]['CSP'].strip() == 'T')
        values['Intensity Ratio_bool'].append(df_bool.iloc[i]['Intensity Ratio'].strip() == 'T')
        values['CleanEx_Holo_bool'].append(df_bool.iloc[i]['CleanEx_Holo'].strip() == 'T')
        values['CleanEx_Apo_bool'].append(df_bool.iloc[i]['CleanEx_Apo'].strip() == 'T')

    # Remove placeholder '_' and join to string
    seq_str = ''.join([aa for aa in sequence if aa != '_'])
    # Instantiate SubSequence (assuming chain_id='A', is_protein=True, is_peptide=False)
    subseq = SubSequence(
        sequence=seq_str,
        chain_id='A',
        is_protein=True,
        is_peptide=False,
        residue_numbers=residue_numbers,
        values=values
    )
    # Return Sequence object
    return Sequence([subseq])


def generate_score(pdb_path, data_seq, apo_pdb_file, DEBUG = False):
    data_subseq = data_seq.subsequences[0]

    # get occupancy of pdb-file
    occluded_nh_0A, sequence_obj_0A, delta_sasa_array_0A = find_occluded_nh_groups_cleanex(pdb_path, apo_pdb_file, sasa_threshold=0.0, DEBUG=DEBUG)
    sequence_obj_0A.subsequences[0].values['CleanEx'] = occluded_nh_0A
    sequence_obj_0A.subsequences[0].values['Delta SASA'] = delta_sasa_array_0A
    subseq_0A = sequence_obj_0A.subsequences[0]
    DELTA_SASA_subseq = find_buried_residues_sasa(pdb_path).subsequences[0]

    def get_csp_score(delta_sasa_subseq, data_subseq, DEBUG = False):
        alignment = Sequence.align_subsequences(delta_sasa_subseq, data_subseq)

        sasa_sequence = alignment[0]
        exp_sequence = alignment[1]
        sasa_buried_array = alignment[2]['buried_array']
        csp_bool = alignment[3]['CSP_bool']
        
        # mask out prolines
        for i in range(len(sasa_sequence)):
            if sasa_sequence[i] == 'P':
                sasa_buried_array[i] = False
                
        if DEBUG:
            print("SASA Sequence | Buried Array | Exp Sequence | CSP Bool")
            print("-" * 60)
            for i in range(len(sasa_sequence)):
                print(f"{sasa_sequence[i]:<13} | {str(sasa_buried_array[i]):<12} | {exp_sequence[i]:<12} | {str(csp_bool[i])}")
        TP = 0
        FP = 0
        FN = 0
        TN = 0


        # get the number of significant CSPs
        for i in range(len(sasa_sequence)):
            if sasa_buried_array[i] == True:
                if csp_bool[i] == True:
                    TP += 1
                else:
                    FP += 1
            else:
                if csp_bool[i] == True:
                    FN += 1 
                else:
                    TN += 1
        F1 = 2 * TP / (2 * TP + FP + FN)
        
        if DEBUG:
            print("TP, FP, FN, TN")
            print(TP, FP, FN, TN)
            print("F1")
            print(F1)
            print("-" * 60)


        return F1

    def get_intensity_ratio_score(delta_sasa_subseq, data_subseq, DEBUG = False):
        alignment = Sequence.align_subsequences(delta_sasa_subseq, data_subseq)

        sasa_sequence = alignment[0]
        exp_sequence = alignment[1]
        sasa_buried_array = alignment[2]['buried_array']
        intensity_ratio_bool = alignment[3]['Intensity Ratio_bool']
        
        # mask out prolines
        for i in range(len(sasa_sequence)):
            if sasa_sequence[i] == 'P':
                sasa_buried_array[i] = False
                
        if DEBUG:
            print("SASA Sequence | Buried Array | Exp Sequence | Intensity Ratio Bool")
            print("-" * 60)
            for i in range(len(sasa_sequence)):
                print(f"{sasa_sequence[i]:<13} | {str(sasa_buried_array[i]):<12} | {exp_sequence[i]:<12} | {str(intensity_ratio_bool[i])}")
        TP = 0
        FP = 0
        FN = 0
        TN = 0


        # get the number of significant CSPs
        for i in range(len(sasa_sequence)):
            if sasa_buried_array[i] == True:
                if intensity_ratio_bool[i] == True:
                    TP += 1
                else:
                    FP += 1
            else:
                if intensity_ratio_bool[i] == True:
                    FN += 1 
                else:
                    TN += 1

        F1 = 2 * TP / (2 * TP + FP + FN)
        if DEBUG:
            print("TP, FP, FN, TN")
            print(TP, FP, FN, TN)
            print("F1")
            print(F1)
            print("-" * 60)

        return F1

    def get_cleanex_score(delta_nh_sasa_subseq, data_subseq, DEBUG = False):
        print(delta_nh_sasa_subseq)
        print(data_subseq)
        alignment = Sequence.align_subsequences(delta_nh_sasa_subseq, data_subseq)

        sasa_sequence = alignment[0]
        exp_sequence = alignment[1]
        sasa_cleanex_bool = alignment[2]['CleanEx']
        delta_sasa_array = alignment[2]['Delta SASA']
        exp_cleanex_holo_bool = alignment[3]['CleanEx_Holo_bool']
        exp_cleanex_apo_bool = alignment[3]['CleanEx_Apo_bool']
        exp_cleanex_delta_bool = []
        for i in range(len(exp_cleanex_holo_bool)):
            if exp_cleanex_holo_bool[i] and not(exp_cleanex_apo_bool[i]):
                exp_cleanex_delta_bool.append(True)
            elif not(exp_cleanex_holo_bool[i]) and exp_cleanex_apo_bool[i]:
                exp_cleanex_delta_bool.append(True)
            else:
                exp_cleanex_delta_bool.append(False)
            
        if DEBUG:
            print("SASA Sequence")
            print(sasa_sequence)
            print("Exp Sequence")
            print(exp_sequence)
            print("SASA CleanEx")
            print(sasa_cleanex_bool)
            print("Delta SASA")
            print(delta_sasa_array)
            print("SASA Sequence | Sasa CleanEx | Delta SASA | Exp Sequence | CleanEx Holo Bool | CleanEx Apo Bool | CleanEx Delta Bool")
            print("-" * 60)
            for i in range(len(sasa_sequence)):
                # Check that none of the values for this index are None before printing
                values_to_check = [
                    sasa_sequence[i],
                    sasa_cleanex_bool[i],
                    delta_sasa_array[i] if sasa_sequence[i] != 'P' else "NONE",
                    exp_sequence[i],
                    exp_cleanex_holo_bool[i],
                    exp_cleanex_apo_bool[i],
                    exp_cleanex_delta_bool[i]
                ]
                if all(v is not None for v in values_to_check):
                    if sasa_sequence[i] == 'P':
                        print(f"{sasa_sequence[i]:<13} | {str(sasa_cleanex_bool[i]):<12} | {'NONE':<12} | {exp_sequence[i]:<12} \
                              | {str(exp_cleanex_holo_bool[i])} | {str(exp_cleanex_apo_bool[i])} | {str(exp_cleanex_delta_bool[i])}")
                    else:
                        print(f"{sasa_sequence[i]:<13} | {str(sasa_cleanex_bool[i]):<12} | {delta_sasa_array[i]:.6f} | {exp_sequence[i]:<12} | \
                              {str(exp_cleanex_holo_bool[i])} | {str(exp_cleanex_apo_bool[i])} | {str(exp_cleanex_delta_bool[i])}")

        Trues = 0
        Falses = 0
        for i in range(len(sasa_sequence)):
            if exp_cleanex_delta_bool[i] == True:
                if sasa_cleanex_bool[i] == True:
                    Trues += 1
                else:
                    Falses += 1
        cleanex_precision_score = 2 * Trues / (2 * Trues + Falses)
        return cleanex_precision_score
    
    csp_score = get_csp_score(DELTA_SASA_subseq, data_subseq, DEBUG=DEBUG)
    intensity_ratio_score = get_intensity_ratio_score(DELTA_SASA_subseq, data_subseq, DEBUG=DEBUG)
    cleanex_precision_score = get_cleanex_score(subseq_0A, data_subseq, DEBUG=DEBUG)

    if DEBUG: 
        print(csp_score)
        print(intensity_ratio_score)
        print(cleanex_precision_score)

    score = csp_score + intensity_ratio_score + cleanex_precision_score

    score_dict = {'CSP': csp_score, 'Intensity Ratio': intensity_ratio_score, 'CleanEx': cleanex_precision_score, 'Cumulative Score': score}
    return score_dict

def main():
    args = parse_arguments()
    pdb_id = args.pdb_id.upper()
    specific_pdb_file = args.pdb_file
    
    # Define paths based on PDB ID
    data_dir = f'./{pdb_id}/'
    bool_data_path = f'{data_dir}{pdb_id}_BOOL_DATA.csv'
    exp_data_path = f'{data_dir}{pdb_id}_DATA.csv'
    scores_csv_path = f'{data_dir}{pdb_id}_SCORES.csv'
    pdb_dir = f'./PDB_FILES/{pdb_id}_processed/'
    # pdb_dir = './selected_models/'
    apo_pdb_file = './30782.pdb'
    
    # Check if required files exist
    if not os.path.exists(bool_data_path):
        print(f"Error: Boolean data file not found at {bool_data_path}")
        sys.exit(1)
    
    if not os.path.exists(exp_data_path):
        print(f"Error: Experimental data file not found at {exp_data_path}")
        sys.exit(1)
    
    if not os.path.exists(pdb_dir):
        print(f"Error: PDB directory not found at {pdb_dir}")
        sys.exit(1)
    
    # Load csv files
    df_bool = pd.read_csv(bool_data_path)
    df_exp = pd.read_csv(exp_data_path)
    
    data_seq = get_data_seq(df_exp, df_bool)
    
    # If specific PDB file is provided, score only that file and print to stdout
    if specific_pdb_file:
        pdb_path = os.path.join(pdb_dir, specific_pdb_file)
        
        if not os.path.exists(pdb_path):
            print(f"Error: PDB file not found at {pdb_path}")
            sys.exit(1)
        
        try:
            # Generate a score for the specific pdb file
            score_dict = generate_score(pdb_path, data_seq, apo_pdb_file, DEBUG=True)
            
            # Print scores to stdout
            print(f"PDB File: {specific_pdb_file}")
            print(f"CSP Score: {score_dict['CSP']:.6f}")
            print(f"Intensity Ratio Score: {score_dict['Intensity Ratio']:.6f}")
            print(f"CleanEx Score: {score_dict['CleanEx']:.6f}")
            print(f"Cumulative Score: {score_dict['Cumulative Score']:.6f}")
            
        except Exception as e:
            print(f"Error processing {specific_pdb_file}: {e}")
            sys.exit(1)
        
        return
    
    # Load existing scores if CSV file exists
    existing_scores = {}
    if os.path.exists(scores_csv_path):
        try:
            df_existing = pd.read_csv(scores_csv_path)
            for _, row in df_existing.iterrows():
                existing_scores[row['pdb_file']] = {
                    'CSP': row['CSP'],
                    'Intensity Ratio': row['Intensity Ratio'],
                    'CleanEx': row['CleanEx'],
                    'Cumulative Score': row['Cumulative Score']
                }
            print(f"Loaded {len(existing_scores)} existing scores from {scores_csv_path}")
        except Exception as e:
            print(f"Error loading existing scores: {e}")
            existing_scores = {}
    
    # Loop through pdb files
    pdb_files = [f for f in os.listdir(pdb_dir) if f.endswith('.pdb')]
    
    # Filter out already processed files
    new_pdb_files = [f for f in pdb_files if f not in existing_scores]
    print(f"Found {len(pdb_files)} total PDB files, {len(new_pdb_files)} new files to process")
    
    from tqdm import tqdm
    
    for pdb_file in tqdm(new_pdb_files, desc="Processing PDB files"):
        pdb_path = os.path.join(pdb_dir, pdb_file)
        
        try:
            # Generate a score for the pdb file using exp_data and bool_data
            score_dict = generate_score(pdb_path, data_seq, apo_pdb_file, DEBUG=False)
            
            # Convert dictionary to list format for CSV writing
            score_row = [pdb_file, score_dict['CSP'], score_dict['Intensity Ratio'], 
                         score_dict['CleanEx'], score_dict['Cumulative Score']]
            
            # Add to existing scores
            existing_scores[pdb_file] = score_dict
            
            # Update CSV file after each calculation
            all_scores = []
            for filename, scores in existing_scores.items():
                all_scores.append([filename, scores['CSP'], scores['Intensity Ratio'], 
                                  scores['CleanEx'], scores['Cumulative Score']])
            
            write_to_csv(scores_csv_path, headers=['pdb_file', 'CSP', 'Intensity Ratio', 'CleanEx', 'Cumulative Score'], data=all_scores)
            
        except Exception as e:
            print(f"Error processing {pdb_file}: {e}")
            continue
    
    print(f"Processing complete. Total scores in {scores_csv_path}: {len(existing_scores)}")

if __name__ == "__main__":
    main()
