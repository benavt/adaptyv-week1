#!/usr/bin/env python3
"""
Script to convert PDB files to distograms.
A distogram is an NxN matrix where each element (i,j) represents the distance between residues i and j.
Processes PDB files based on a PDB ID.
"""
import os
import sys
import numpy as np
from Bio import PDB
from Bio.PDB.PDBParser import PDBParser
import argparse
from pathlib import Path
from tqdm import tqdm
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Sequence import Sequence

def get_ca_atoms(pdb_file, keep_chain_resnums=None):
    """
    Extract CA atoms from a PDB file, optionally filtering by (chain_id, residue_number).
    Args:
        pdb_file: Path to the PDB file
        keep_chain_resnums: List of (chain_id, residue_number) tuples to keep (in order)
    Returns:
        List of CA atom coordinates (in the order of keep_chain_resnums if provided)
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    ca_atoms = []
    ca_chain_resnums = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if 'CA' in residue:
                    ca_atoms.append(residue['CA'].get_coord())
                    ca_chain_resnums.append((chain.id, residue.get_id()[1]))
    if keep_chain_resnums is not None:
        # Build a mapping from (chain_id, residue_number) to CA atom
        resmap = {cr: atom for cr, atom in zip(ca_chain_resnums, ca_atoms)}
        filtered_atoms = []
        for cr in keep_chain_resnums:
            if cr in resmap:
                filtered_atoms.append(resmap[cr])
            else:
                print(f"[DEBUG] (chain, resnum) {cr} not found in {pdb_file}, skipping.")
        return np.array(filtered_atoms)
    return np.array(ca_atoms)

def calculate_distogram(ca_atoms):
    """
    Calculate the distogram from CA atom coordinates.
    
    Args:
        ca_atoms: Array of CA atom coordinates
        
    Returns:
        NxN matrix of distances between residues
    """
    n = len(ca_atoms)
    distogram = np.zeros((n, n))
    
    for i in range(n):
        for j in range(n):
            distogram[i, j] = np.linalg.norm(ca_atoms[i] - ca_atoms[j])
    
    return distogram

def process_pdb_file(pdb_file, output_dir, keep_chain_resnums=None):
    """
    Process a single PDB file and save its distogram.
    Args:
        pdb_file: Path to the PDB file
        output_dir: Directory to save the distogram
        keep_chain_resnums: List of (chain_id, residue_number) tuples to keep (in order)
    """
    try:
        # Get CA atoms
        ca_atoms = get_ca_atoms(pdb_file, keep_chain_resnums=keep_chain_resnums)
        if len(ca_atoms) == 0:
            print(f"[DEBUG] No CA atoms found for selected residues in {pdb_file}")
            return
        # Calculate distogram
        distogram = calculate_distogram(ca_atoms)
        # Get basename of PDB file
        basename = os.path.splitext(os.path.basename(pdb_file))[0]
        # Save distogram
        output_file = os.path.join(output_dir, f"{basename}.npy")
        np.save(output_file, distogram)
        print(f"Distogram saved to {output_file}")
    except Exception as e:
        print(f"Error processing {pdb_file}: {e}")

def get_pdb_paths(directory):
    """
    Get paths to PDB files for a given PDB ID.
    Searches in:
    1. All files in directories matching 'directory/*.pdb'
    
    Args:
        directory: The directory to search for PDB files
        
    Returns:
        List of paths to PDB files
    """
    pdb_files = []
    base_dir = Path(directory)
    print(f"[DEBUG] Searching for PDB files in {base_dir}")

    return [str(f) for f in base_dir.glob('*.pdb')]

def get_available_pdb_dirs(pdb_id):
    """
    Get available PDB file directories for a given PDB ID.
    Searches in ./PDB_FILES/{pdb_id.upper()}*/
    
    Args:
        pdb_id: The PDB ID to search for
        
    Returns:
        List of available directory paths
    """
    pdb_id_upper = pdb_id.upper()
    search_pattern = Path('./PDB_FILES') / f'{pdb_id_upper}*'
    available_dirs = list(search_pattern.parent.glob(search_pattern.name))
    return available_dirs

def get_output_dir(pdb_id):
    """
    Get the output directory for a given PDB ID.
    Creates path in format: './<pdb_id>.upper()/<pdb_id>_distograms/'
    
    Args:
        pdb_id: The PDB ID to process
        
    Returns:
        Path to the output directory
    """
    pdb_id_upper = pdb_id.upper()
    output_dir = Path(f'./{pdb_id_upper}/{pdb_id}_distograms')
    return str(output_dir)

def main():
    parser = argparse.ArgumentParser(description='Convert PDB file(s) to distogram(s)')
    parser.add_argument('pdb_id', help='PDB ID to process (e.g., 1abc)')
    args = parser.parse_args()
    
    # Get available PDB directories
    available_dirs = get_available_pdb_dirs(args.pdb_id)
    
    if not available_dirs:
        print(f"No PDB directories found for {args.pdb_id}")
        sys.exit(1)
        
    # Select directory
    if len(available_dirs) == 1:
        selected_dir = available_dirs[0]
        print(f"Using directory: {selected_dir}")
    else:
        print("\nAvailable PDB directories:")
        for i, dir_path in enumerate(available_dirs, 1):
            print(f"{i}. {dir_path}")
        while True:
            try:
                choice = int(input("\nSelect a directory number: "))
                if 1 <= choice <= len(available_dirs):
                    selected_dir = available_dirs[choice - 1]
                    break
                else:
                    print(f"Please enter a number between 1 and {len(available_dirs)}")
            except ValueError:
                print("Please enter a valid number")
    
    # Get input paths and output directory
    pdb_files = get_pdb_paths(selected_dir)
    output_dir = get_output_dir(args.pdb_id)
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    if not pdb_files:
        print(f"No PDB files found in {selected_dir}")
        sys.exit(1)
    
    print(f"Found {len(pdb_files)} PDB files to process")

    # --- Sequence Extraction and Alignment ---
    sequences = []
    seq_objs = []
    for pdb_file in pdb_files:
        # print(f"[DEBUG] Extracting sequence from: {pdb_file}")
        seq_obj = Sequence.from_pdb_single(pdb_file)
        seq_objs.append(seq_obj)
        # Concatenate protein and peptide chain sequences
        seq_str = ''.join([s.sequence for s in seq_obj.protein_sequences + seq_obj.peptide_sequences])
        sequences.append(seq_str)
        # print(f"[DEBUG] Concatenated sequence (protein + peptide chains): {seq_str}")

    # Print all unique sequences
    unique_sequences = set(sequences)
    print("\n[DEBUG] Unique sequences found:")
    for us in unique_sequences:
        print(us)
    # user_input = input("\nIs it okay to continue with these unique sequences? (y/n): ")
    # if user_input.lower() != 'y':
    #     print("Exiting.")
    #     sys.exit(0)

    # --- Alignment ---
    # Build concatenated SubSequence objects for alignment
    concat_subseqs = []
    concat_chain_resnums_list = []
    from tqdm import tqdm
    for seq_obj in tqdm(seq_objs, desc="Building concatenated SubSequence objects"):
        # Get protein and peptide sequences
        protein_seqs = seq_obj.protein_sequences
        peptide_seqs = seq_obj.peptide_sequences
        
        # Build protein concatenated sequence
        protein_concat_seq = ''.join([s.sequence for s in protein_seqs])
        protein_chain_resnums = []
        for s in protein_seqs:
            if s.residue_numbers:
                protein_chain_resnums.extend([(s.chain_id, rn) for rn in s.residue_numbers])
        
        # Build peptide concatenated sequence  
        peptide_concat_seq = ''.join([s.sequence for s in peptide_seqs])
        peptide_chain_resnums = []
        for s in peptide_seqs:
            if s.residue_numbers:
                peptide_chain_resnums.extend([(s.chain_id, rn) for rn in s.residue_numbers])

        # Create SubSequence objects
        from Sequence import SubSequence
        protein_subseq = SubSequence(
            sequence=protein_concat_seq,
            chain_id='PROT',
            is_protein=True,
            is_peptide=False,
            residue_numbers=protein_chain_resnums
        )
        peptide_subseq = SubSequence(
            sequence=peptide_concat_seq, 
            chain_id='PEPT',
            is_protein=False,
            is_peptide=True,
            residue_numbers=peptide_chain_resnums
        )
        
        # Store tuple of protein and peptide sequences
        concat_subseqs.append((protein_subseq, peptide_subseq))
        
        # Store combined chain residue numbers
        concat_chain_resnums = protein_chain_resnums + peptide_chain_resnums
        # print(f"[DEBUG] Concatenated chain_resnums: {concat_chain_resnums}")
        concat_chain_resnums_list.append(concat_chain_resnums)

    
    # user_input = input("\nIs it okay to continue with these concatenated sequences? (y/n): ")
    # if user_input.lower() != 'y':
    #     print("Exiting.")
    #     sys.exit(0)
        
    # Separate protein and peptide sequences
    protein_subseqs = [pair[0] for pair in concat_subseqs]
    peptide_subseqs = [pair[1] for pair in concat_subseqs]

    def align_chain_resnums(aligned_seq1, aligned_seq2, orig_chain_resnums):
        """
        Align chain/residue numbers to aligned sequences, only keeping positions where both sequences match.
        Args:
            aligned_seq1: First aligned sequence
            aligned_seq2: Second aligned sequence
            orig_chain_resnums: Original chain/residue numbers for the first sequence
        Returns:
            List of chain/residue numbers, with None for any position where either sequence has a gap
        """
        aligned = []
        idx = 0
        for c1, c2 in zip(aligned_seq1, aligned_seq2):
            if c1 == '_' or c2 == '_':  # If either sequence has a gap, skip this position
                continue
            aligned.append(orig_chain_resnums[idx])
            idx += 1
        return aligned

    def create_aligned_subseq(orig_subseq, aligned_seq1, aligned_seq2, aligned_chain_resnums):
        """Create a new SubSequence with aligned sequence and chain/residue numbers."""
        # Create a new sequence that only includes matching positions
        new_seq = ''.join(c1 for c1, c2 in zip(aligned_seq1, aligned_seq2) if c1 != '_' and c2 != '_')
        return SubSequence(
            sequence=new_seq,
            chain_id=orig_subseq.chain_id,
            is_protein=orig_subseq.is_protein,
            is_peptide=orig_subseq.is_peptide,
            residue_numbers=aligned_chain_resnums,
            values=orig_subseq.values
        )

    def align_sequences(sequences, seq_type):
        """Align a list of sequences against the first sequence."""
        print(f"\nAligning {seq_type} sequences:")
        aligned_subseqs = [sequences[0]]
        ref_seq = sequences[0]
        
        for i in tqdm(range(1, len(sequences)), desc=f"Aligning {seq_type} sequences"):
            target_seq = sequences[i]
            # print(f"[DEBUG] Aligning {seq_type} sequence 0 and {i}")
            aligned_ref, aligned_target, _, _, score = Sequence.align_subsequences(ref_seq, target_seq)
            # print(f"[DEBUG] Alignment score: {score}")
            # print(f"[DEBUG] Aligned 0: {aligned_ref}")
            # print(f"[DEBUG] Aligned {i}: {aligned_target}")
            
            if i == 1:
                # For the reference sequence, align against itself to get matching positions
                aligned_ref_chain_resnums = align_chain_resnums(aligned_ref, aligned_ref, ref_seq.residue_numbers or [])
                aligned_subseqs[0] = create_aligned_subseq(ref_seq, aligned_ref, aligned_ref, aligned_ref_chain_resnums)
            
            # For target sequence, align against reference to get matching positions
            aligned_target_chain_resnums = align_chain_resnums(aligned_ref, aligned_target, target_seq.residue_numbers or [])
            aligned_subseqs.append(create_aligned_subseq(target_seq, aligned_ref, aligned_target, aligned_target_chain_resnums))
            
            # print(f"[DEBUG] Aligned chain_resnums for {seq_type} sequence {i}: {aligned_target_chain_resnums}")
            
        return aligned_subseqs

    # Align protein and peptide sequences
    aligned_protein_subseqs = align_sequences(protein_subseqs, "protein")
    print("\nUnique protein alignments:")
    unique_protein_seqs = [(i, seq.sequence, seq.residue_numbers) for i, seq in enumerate(aligned_protein_subseqs)]
    unique_seen = set()
    for i, seq, resnums in unique_protein_seqs:
        if seq not in unique_seen:
            print(f"Unique protein alignment {i}:")
            print(f"Sequence:  {seq}")
            print(f"Residues:  {resnums}")
            unique_seen.add(seq)

        
    aligned_peptide_subseqs = align_sequences(peptide_subseqs, "peptide") 
    print("\nUnique peptide alignments:")
    unique_peptide_seqs = [(i, seq.sequence, seq.residue_numbers) for i, seq in enumerate(aligned_peptide_subseqs)]
    unique_seen = set()
    for i, seq, resnums in unique_peptide_seqs:
        if seq not in unique_seen:
            print(f"Unique peptide alignment {i}:")
            print(f"Sequence:  {seq}")
            print(f"Residues:  {resnums}")
            unique_seen.add(seq)

    # Combine aligned protein and peptide sequences
    aligned_subseqs = list(zip(aligned_protein_subseqs, aligned_peptide_subseqs))

    # After aligning sequences, find common regions to trim to
    def find_common_regions(subseqs):
        """
        Find the longest common substring across all sequences of the same type (protein or peptide).
        Returns a list of lists containing a single (start_resnum, end_resnum) tuple, where:
        - The outer list has one entry per sequence in subseqs
        - Each inner list contains a single tuple of (start_resnum, end_resnum) for that sequence's longest common region
        - The residue numbers are taken from each sequence's residue_numbers list
        """
        if not subseqs:
            return []
            
        # Get all sequences and their residue numbers
        all_seqs = [s.sequence for s in subseqs]
        all_resnums = [s.residue_numbers for s in subseqs]
        if not all_seqs:
            return []
            
        # Quick check if all sequences are identical
        ref_seq = all_seqs[0]
        if all(seq == ref_seq for seq in all_seqs):
            # All sequences are identical, return single region for each sequence
            print(f"[DEBUG] All sequences are identical, returning single region for each sequence")
            ret_val = []
            for i, resnums in enumerate(all_resnums):
                print(f"[DEBUG] resnums: {resnums}, i: {i}")
                ret_val.append([(resnums[0], resnums[-1])])
            return ret_val
            
        # Find longest common substring using dynamic programming
        def find_lcs(sequences):
            if not sequences:
                return ""
            
            # Use first sequence as reference
            ref = sequences[0]
            max_len = 0
            max_start = 0
            
            # For each possible start position in reference sequence
            for i in range(len(ref)):
                # For each possible length from this start position
                for j in range(i + 1, len(ref) + 1):
                    substr = ref[i:j]
                    # Check if this substring appears in all other sequences
                    if all(substr in seq for seq in sequences[1:]):
                        if j - i > max_len:
                            max_len = j - i
                            max_start = i
                            
            return ref[max_start:max_start + max_len] if max_len > 0 else ""
        
        # Find longest common substring
        lcs = find_lcs(all_seqs)
        if not lcs:
            print("[ERROR] No common substring found across sequences!")
            return []
            
        print(f"[DEBUG] Found longest common substring: {lcs}")
        
        # Convert LCS to residue number regions for each sequence
        common_regions = []
        for seq, resnums in zip(all_seqs, all_resnums):
            # Find the start and end positions of LCS in this sequence
            start_pos = seq.find(lcs)
            if start_pos == -1:  # Should never happen since LCS is common
                print(f"[ERROR] LCS not found in sequence: {seq}")
                return []
            end_pos = start_pos + len(lcs) - 1
            
            # Convert to residue numbers
            start_resnum = resnums[start_pos]
            end_resnum = resnums[end_pos]
            common_regions.append([(start_resnum, end_resnum)])
            
        return common_regions

    def trim_to_common_regions(subseqs, common_regions):
        """
        Trim SubSequence objects to only include common regions.
        Args:
            subseqs: List of SubSequence objects
            common_regions: List of lists of (start_resnum, end_resnum) tuples
        Returns:
            List of trimmed SubSequence objects
        """
        trimmed_subseqs = []
        for subseq, regions in zip(subseqs, common_regions):
            # Debug: print sequence and its common regions
            print(f"Debug - Processing sequence: {subseq.sequence}")
            print(f"Debug - Common regions: {regions}")
            # Get the sequence and residue numbers
            seq = subseq.sequence
            resnums = subseq.residue_numbers
            
            # Build new sequence and residue numbers from common regions
            new_seq = []
            new_resnums = []
            
            # For each region in this sequence
            for start_resnum, end_resnum in regions:
                # Find the indices in the original sequence
                start_idx = resnums.index(start_resnum)
                end_idx = resnums.index(end_resnum)
                # Add the sequence and residue numbers for this region
                new_seq.extend(seq[start_idx:end_idx + 1])
                new_resnums.extend(resnums[start_idx:end_idx + 1])
                
            # Create new SubSequence with trimmed data
            trimmed_subseqs.append(SubSequence(
                sequence=''.join(new_seq),
                chain_id=subseq.chain_id,
                is_protein=subseq.is_protein,
                is_peptide=subseq.is_peptide,
                residue_numbers=new_resnums,
                values=subseq.values
            ))
            
        return trimmed_subseqs

    def process_sequence_type(subseqs, seq_type):
        """
        Process a list of subsequences of the same type (protein or peptide) to find common regions,
        trim to those regions, and verify the results.
        
        Args:
            subseqs: List of SubSequence objects of the same type
            seq_type: String describing the sequence type ("protein" or "peptide")
            
        Returns:
            List of trimmed SubSequence objects
        """
        # Find common regions
        common_regions = find_common_regions(subseqs)
        if not common_regions or not any(regions for regions in common_regions):
            print(f"[ERROR] No common regions found across {seq_type} sequences!")
            sys.exit(1)
            
        print(f"\nCommon regions found in {seq_type} sequences:")
        for i, (subseq, regions) in enumerate(zip(subseqs, common_regions)):
            print(f"\nSequence {i}:")
            for j, (start_resnum, end_resnum) in enumerate(regions):
                # Find the sequence segment for this region
                start_idx = subseq.residue_numbers.index(start_resnum)
                end_idx = subseq.residue_numbers.index(end_resnum)
                seq_segment = subseq.sequence[start_idx:end_idx + 1]
                print(f"Region {j+1}: residues {start_resnum}-{end_resnum} (length {end_idx-start_idx+1})")
                print(f"Sequence: {seq_segment}")
            # raise Exception("STOP HERE")
        
        # Trim sequences to common regions
        trimmed_subseqs = trim_to_common_regions(subseqs, common_regions)
        
        # Verify all trimmed sequences are identical
        trimmed_seqs = [s.sequence for s in trimmed_subseqs]
        if len(set(trimmed_seqs)) != 1:
            print(f"[ERROR] Not all trimmed {seq_type} sequences are identical!")
            for i, seq in enumerate(trimmed_seqs):
                print(f"{seq_type.capitalize()} sequence {i}: {seq}")
            sys.exit(1)
        #print(f"{trimmed_subseqs[0].residue_numbers}")
        #print(f"{trimmed_subseqs[0].sequence}")
        #raise Exception("STOP HERE")
        return trimmed_subseqs

    # Process protein sequences
    trimmed_protein_subseqs = process_sequence_type(aligned_protein_subseqs, "protein")
    
    # Create final list of trimmed (protein, peptide) pairs
    trimmed_subseqs = trimmed_protein_subseqs
    from tqdm import tqdm
    # Continue with distogram processing using trimmed sequences
    for i, protein_subseq in enumerate(tqdm(trimmed_subseqs, desc="Processing PDB files")):
        # Combine protein and peptide chain/residue numbers
        keep_chain_resnums = protein_subseq.residue_numbers 
        # print(f"[DEBUG] For file {i}, keeping (chain, resnum): {keep_chain_resnums}")
        process_pdb_file(pdb_files[i], output_dir, keep_chain_resnums=keep_chain_resnums)

if __name__ == "__main__":
    main() 