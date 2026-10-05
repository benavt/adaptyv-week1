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
import tkinter as tk
from tkinter import ttk
import readline  # For tab completion in prompt_residue_ranges

def get_selected_atoms(pdb_file, atom_type='CA', keep_chain_resnums=None):
    """
    Extract selected atoms from a PDB file, optionally filtering by (chain_id, residue_number).
    Args:
        pdb_file: Path to the PDB file
        atom_type: 'CA' for alpha carbons, 'H' for all hydrogens
        keep_chain_resnums: List of (chain_id, residue_number) tuples to keep (in order)
    Returns:
        List of atom coordinates (in the order of keep_chain_resnums if provided)
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    atom_coords = []
    atom_chain_resnums = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom_type == 'CA' and atom.get_id() == 'CA':
                        atom_coords.append(atom.get_coord())
                        atom_chain_resnums.append((chain.id, residue.get_id()[1]))
                    elif atom_type == 'H' and atom.element == 'H':
                        atom_coords.append(atom.get_coord())
                        atom_chain_resnums.append((chain.id, residue.get_id()[1]))
    if keep_chain_resnums is not None:
        resmap = {}
        for cr, coord in zip(atom_chain_resnums, atom_coords):
            if cr not in resmap:
                resmap[cr] = []
            resmap[cr].append(coord)
        filtered_atoms = []
        for cr in keep_chain_resnums:
            if cr in resmap:
                if atom_type == 'CA':
                    filtered_atoms.append(resmap[cr][0])
                elif atom_type == 'H':
                    filtered_atoms.extend(resmap[cr])
            else:
                print(f"[DEBUG] (chain, resnum) {cr} not found in {pdb_file}, skipping.")
        return np.array(filtered_atoms)
    return np.array(atom_coords)

def calculate_distogram(atom_coords):
    """
    Calculate the distogram from atom coordinates.
    
    Args:
        atom_coords: Array of atom coordinates
        
    Returns:
        NxN matrix of distances between residues
    """
    n = len(atom_coords)
    distogram = np.zeros((n, n))
    
    for i in range(n):
        for j in range(n):
            distogram[i, j] = np.linalg.norm(atom_coords[i] - atom_coords[j])
    
    return distogram

def process_pdb_file(pdb_file, output_dir, keep_chain_resnums=None, atom_type='CA'):
    """
    Process a single PDB file and save its distogram.
    Args:
        pdb_file: Path to the PDB file
        output_dir: Directory to save the distogram
        keep_chain_resnums: List of (chain_id, residue_number) tuples to keep (in order)
        atom_type: 'CA' for alpha carbons, 'H' for all hydrogens
    """
    try:
        atom_coords = get_selected_atoms(pdb_file, atom_type=atom_type, keep_chain_resnums=keep_chain_resnums)
        if len(atom_coords) == 0:
            print(f"[DEBUG] No {atom_type} atoms found for selected residues in {pdb_file}")
            return
        distogram = calculate_distogram(atom_coords)
        basename = os.path.splitext(os.path.basename(pdb_file))[0]
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
    parser.add_argument('directory', help='Directory to process')
    parser.add_argument('output_dir', help='Output directory')
    parser.add_argument('--atom_type', choices=['CA', 'H'], default=None, help="Atom type to use for distogram: 'CA' or 'H'")
    args = parser.parse_args()
    
    # Get input paths and output directory
    pdb_files = get_pdb_paths(args.directory)
    output_dir = get_output_dir(args.output_dir)
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    if not pdb_files:
        print(f"No PDB files found for {args.directory}")
        sys.exit(1)
    
    print(f"Found {len(pdb_files)} PDB files to process")
    print("[DEBUG] Finished input/output directory setup.")

    # --- Sequence Extraction and Alignment ---
    sequences = []
    seq_objs = []
    from tqdm import tqdm
    for pdb_file in tqdm(pdb_files, desc="Extracting sequences"):
        seq_obj = Sequence.from_pdb_single(pdb_file)
        seq_objs.append(seq_obj)
        seq_str = ''.join([s.sequence for s in seq_obj.protein_sequences + seq_obj.peptide_sequences])
        sequences.append(seq_str)
    print("[DEBUG] Finished sequence extraction.")

    unique_sequences = set(sequences)
    print("\n[DEBUG] Unique sequences found:")
    for us in unique_sequences:
        print(us)

    # Check if all sequences are identical
    if len(unique_sequences) == 1:
        print("[DEBUG] All sequences are identical. Skipping alignment step.")
        # Build sequence_groups directly from original seq_objs
        sequence_groups = []
        for i, seq_obj in enumerate(seq_objs):
            if len(seq_obj.peptide_sequences) == 0:
                sequence_groups.append({
                    'sequence': seq_obj.protein_sequences[0].sequence,
                    'residue_numbers': seq_obj.protein_sequences[0].residue_numbers,
                    'chain_ids': [seq_obj.protein_sequences[0].chain_id]*len(seq_obj.protein_sequences[0].residue_numbers),
                    'files': [pdb_files[i]]
                })                
                continue
            if seq_obj.protein_sequences[0].sequence and seq_obj.peptide_sequences[0].sequence :
                sequence_groups.append({
                    'sequence': seq_obj.protein_sequences[0].sequence + seq_obj.peptide_sequences[0].sequence,
                    'residue_numbers': seq_obj.protein_sequences[0].residue_numbers + seq_obj.peptide_sequences[0].residue_numbers,
                    'chain_ids': [seq_obj.protein_sequences[0].chain_id]*len(seq_obj.protein_sequences[0].residue_numbers) + [seq_obj.peptide_sequences[0].chain_id]*len(seq_obj.peptide_sequences[0].residue_numbers),
                    'files': [pdb_files[i]]
                })
        print("[DEBUG] Built sequence_groups for user selection (no alignment needed).")
        reference_group = sequence_groups[0]
        print(f"[DEBUG] Reference group: {reference_group}")
        print(f"[DEBUG] Reference group sequence: {reference_group['sequence']}")
    else:
        # --- Alignment ---
        concat_subseqs = []
        concat_chain_resnums_list = []
        from tqdm import tqdm
        for seq_obj in tqdm(seq_objs, desc="Building concatenated SubSequence objects"):
            protein_seqs = seq_obj.protein_sequences
            peptide_seqs = seq_obj.peptide_sequences
            protein_concat_seq = ''.join([s.sequence for s in protein_seqs])
            protein_chain_resnums = []
            for s in protein_seqs:
                if s.residue_numbers:
                    protein_chain_resnums.extend([(s.chain_id, rn) for rn in s.residue_numbers])
            peptide_concat_seq = ''.join([s.sequence for s in peptide_seqs])
            peptide_chain_resnums = []
            for s in peptide_seqs:
                if s.residue_numbers:
                    peptide_chain_resnums.extend([(s.chain_id, rn) for rn in s.residue_numbers])
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
            concat_subseqs.append((protein_subseq, peptide_subseq))
            concat_chain_resnums = protein_chain_resnums + peptide_chain_resnums
            concat_chain_resnums_list.append(concat_chain_resnums)
        print("[DEBUG] Finished building concatenated SubSequence objects.")

        protein_subseqs = [pair[0] for pair in concat_subseqs]
        peptide_subseqs = [pair[1] for pair in concat_subseqs]

        def align_chain_resnums(aligned_seq1, aligned_seq2, orig_chain_resnums):
            aligned = []
            idx = 0
            for c1, c2 in zip(aligned_seq1, aligned_seq2):
                if c1 == '_' or c2 == '_':
                    continue
                aligned.append(orig_chain_resnums[idx])
                idx += 1
            return aligned

        def create_aligned_subseq(orig_subseq, aligned_seq1, aligned_seq2, aligned_chain_resnums):
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
            print(f"\nAligning {seq_type} sequences:")
            aligned_subseqs = [sequences[0]]
            ref_seq = sequences[0]
            for i in tqdm(range(1, len(sequences)), desc=f"Aligning {seq_type} sequences"):
                target_seq = sequences[i]
                aligned_ref, aligned_target, _, _, score = Sequence.align_subsequences(ref_seq, target_seq)
                if i == 1:
                    aligned_ref_chain_resnums = align_chain_resnums(aligned_ref, aligned_ref, ref_seq.residue_numbers or [])
                    aligned_subseqs[0] = create_aligned_subseq(ref_seq, aligned_ref, aligned_ref, aligned_ref_chain_resnums)
                aligned_target_chain_resnums = align_chain_resnums(aligned_ref, aligned_target, target_seq.residue_numbers or [])
                aligned_subseqs.append(create_aligned_subseq(target_seq, aligned_ref, aligned_target, aligned_target_chain_resnums))
            return aligned_subseqs

        aligned_protein_subseqs = align_sequences(protein_subseqs, "protein")
        aligned_peptide_subseqs = align_sequences(peptide_subseqs, "peptide")
        print("[DEBUG] Finished sequence alignment.")

        # Build sequence_groups for user selection
        sequence_groups = []
        for i, (protein_subseq, peptide_subseq) in enumerate(zip(aligned_protein_subseqs, aligned_peptide_subseqs)):
            # Always add protein subseq
            sequence_groups.append({
                'sequence': protein_subseq.sequence,
                'residue_numbers': [cr[1] for cr in protein_subseq.residue_numbers],
                'chain_ids': [cr[0] for cr in protein_subseq.residue_numbers],
                'files': [pdb_files[i]]
            })
            # Add peptide subseq if available and has length > 0
            if peptide_subseq and len(peptide_subseq.sequence) > 0:
                sequence_groups.append({
                    'sequence': peptide_subseq.sequence,
                    'residue_numbers': [cr[1] for cr in peptide_subseq.residue_numbers],
                    'chain_ids': [cr[0] for cr in peptide_subseq.residue_numbers],
                    'files': [pdb_files[i]]
                })
        print("[DEBUG] Built sequence_groups for user selection.")
        reference_group = sequence_groups[0]
        print(f"[DEBUG] Reference group: {reference_group}")
        print(f"[DEBUG] Reference group sequence: {reference_group['sequence']}")

    # Prompt user for residue ranges
    selected_ranges = prompt_residue_ranges(sequence_groups, reference_group)
    print(f"[DEBUG] User selected ranges: {selected_ranges}")
    if not selected_ranges:
        print("No ranges selected. Exiting.")
        sys.exit(1)

    def get_selected_chain_resnums(group, selected_ranges):
        chain_ids = group['chain_ids']
        residue_numbers = group['residue_numbers']
        selected = []
        for chain, resnum in zip(chain_ids, residue_numbers):
            for sel_chain, sel_start, sel_end in selected_ranges:
                if chain == sel_chain and sel_start <= resnum <= sel_end:
                    selected.append((chain, resnum))
                    break
        return selected

    # Prompt user for atom type if not provided
    atom_type = args.atom_type
    if atom_type is None:
        while True:
            atom_type = input("Select atom type for distogram ('CA' for alpha carbons, 'H' for all hydrogens): ").strip().upper()
            if atom_type in ['CA', 'H']:
                break
            print("Invalid atom type. Please enter 'CA' or 'H'.")
    print(f"[DEBUG] Using atom type: {atom_type}")

    print("[DEBUG] Starting distogram processing loop.")
    for i, group in enumerate(sequence_groups):
        keep_chain_resnums = get_selected_chain_resnums(group, selected_ranges)
        print(f"[DEBUG] Processing file {group['files'][0]} with {len(keep_chain_resnums)} selected residues.")
        process_pdb_file(group['files'][0], output_dir, keep_chain_resnums=keep_chain_resnums, atom_type=atom_type)

def find_shared_residues(sequence_groups, reference_group):
    """Find residues that are shared among all sequence groups.
    Args:
        sequence_groups: List of sequence group dicts
        reference_group: Reference sequence group dict
    Returns:
        List of indices in the reference sequence that are shared among all groups
    """
    shared_positions = []
    ref_positions = [(c, n) for c, n in zip(reference_group['chain_ids'], reference_group['residue_numbers'])
                    if c != ' ' and n != -1]
    for i, (ref_chain, ref_res) in enumerate(ref_positions):
        is_shared = True
        for group in sequence_groups:
            if group is reference_group:
                continue
            found = False
            for chain, res_num in zip(group['chain_ids'], group['residue_numbers']):
                if chain == ref_chain and res_num == ref_res:
                    found = True
                    break
            if not found:
                is_shared = False
                break
        if is_shared:
            ref_idx = next(i for i, (c, n) in enumerate(zip(reference_group['chain_ids'], reference_group['residue_numbers']))
                         if c == ref_chain and n == ref_res)
            shared_positions.append(ref_idx)
    return sorted(shared_positions)

class RangeSelector(tk.Tk):
    """GUI window for selecting residue ranges."""
    def __init__(self, sequence, residue_numbers, chain_ids, shared_positions):
        super().__init__()
        
        self.title("Select Residue Ranges for Alignment")
        self.sequence = sequence
        self.residue_numbers = residue_numbers
        self.chain_ids = chain_ids
        self.shared_positions = shared_positions
        self.ranges = []  # List of (chain, start, end) tuples
        self.current_selection = None
        
        # Calculate layout parameters
        self.residues_per_row = 100
        self.char_width = 15
        self.row_height = 120
        self.num_rows = (len(sequence) + self.residues_per_row - 1) // self.residues_per_row
        
        # Create main frame
        main_frame = ttk.Frame(self)
        main_frame.pack(padx=10, pady=10, fill=tk.BOTH, expand=True)
        
        # Instructions
        ttk.Label(main_frame, text="Click and drag to select ranges. Only residues shared among all structures are shown.").pack()
        
        # Create scrollable frame for sequence display
        canvas_frame = ttk.Frame(main_frame)
        canvas_frame.pack(fill=tk.BOTH, expand=True)
        
        # Create canvas with scrollbar
        self.canvas = tk.Canvas(canvas_frame, bg='white')
        scrollbar = ttk.Scrollbar(canvas_frame, orient="vertical", command=self.canvas.yview)
        self.canvas.configure(yscrollcommand=scrollbar.set)
        
        # Pack canvas and scrollbar
        self.canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        # Create buttons
        button_frame = ttk.Frame(main_frame)
        button_frame.pack(fill=tk.X, pady=10)
        
        ttk.Button(button_frame, text="Clear Selections", command=self.clear_selections).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Done", command=self.finish_selection).pack(side=tk.RIGHT, padx=5)
        
        # Display sequence
        self.draw_sequence()
        
        # Bind mouse events
        self.canvas.bind('<ButtonPress-1>', self.start_selection)
        self.canvas.bind('<B1-Motion>', self.update_selection)
        self.canvas.bind('<ButtonRelease-1>', self.end_selection)
        
        # Center window and set appropriate size
        self.update_idletasks()
        window_width = min(800, self.residues_per_row * self.char_width + 40)
        window_height = min(600, self.num_rows * self.row_height + 150)
        self.geometry(f'{window_width}x{window_height}')
        
        x = (self.winfo_screenwidth() // 2) - (window_width // 2)
        y = (self.winfo_screenheight() // 2) - (window_height // 2)
        self.geometry(f'+{x}+{y}')
    
    def draw_sequence(self):
        """Draw the sequence on the canvas in multiple rows."""
        self.canvas.delete('all')
        
        # Calculate canvas dimensions
        canvas_width = self.residues_per_row * self.char_width + 20
        canvas_height = self.num_rows * self.row_height
        self.canvas.configure(width=canvas_width, height=canvas_height, scrollregion=(0, 0, canvas_width, canvas_height))
        
        # Draw each row
        for row in range(self.num_rows):
            start_idx = row * self.residues_per_row
            end_idx = min(start_idx + self.residues_per_row, len(self.sequence))
            row_sequence = self.sequence[start_idx:end_idx]
            row_residue_numbers = self.residue_numbers[start_idx:end_idx]
            row_chain_ids = self.chain_ids[start_idx:end_idx]
            row_shared_positions = [i - start_idx for i in self.shared_positions if start_idx <= i < end_idx]
            
            self.draw_sequence_row(row, row_sequence, row_residue_numbers, row_chain_ids, row_shared_positions, start_idx)
        
        # Draw existing selections
        self.draw_selections()
    
    def draw_sequence_row(self, row, sequence, residue_numbers, chain_ids, shared_positions, global_start_idx):
        """Draw a single row of the sequence."""
        y_offset = row * self.row_height
        
        # Draw chain separators and background
        current_x = 10
        current_chain = chain_ids[0] if chain_ids else ' '
        chain_start_x = current_x
        
        # First pass: Draw backgrounds for chains
        for i, (char, chain) in enumerate(zip(sequence, chain_ids)):
            if chain != current_chain:
                # Draw previous chain background
                self.canvas.create_rectangle(
                    chain_start_x, y_offset + 20,
                    current_x, y_offset + 100,
                    fill='lightgray' if current_chain != ' ' else 'white',
                    outline='gray'
                )
                chain_start_x = current_x
                current_chain = chain
            current_x += self.char_width
        
        # Draw last chain background
        self.canvas.create_rectangle(
            chain_start_x, y_offset + 20,
            current_x, y_offset + 100,
            fill='lightgray' if current_chain != ' ' else 'white',
            outline='gray'
        )
        
        # Second pass: Draw text
        current_x = 10
        for i, (char, chain) in enumerate(zip(sequence, chain_ids)):
            x_center = current_x + self.char_width/2
            global_idx = global_start_idx + i
            
            # Only draw residues that are shared
            if global_idx in self.shared_positions:
                # Draw chain ID at top
                self.canvas.create_text(
                    x_center, y_offset + 30,
                    text=chain,
                    font=('Courier', 10),
                    fill='blue'
                )
                
                # Draw residue number in middle
                self.canvas.create_text(
                    x_center, y_offset + 50,
                    text=str(residue_numbers[i]),
                    font=('Courier', 8)
                )
                
                # Draw amino acid at bottom
                self.canvas.create_text(
                    x_center, y_offset + 80,
                    text=char,
                    font=('Courier', 14, 'bold')
                )
            elif char == ':':
                # Draw chain break marker
                self.canvas.create_line(
                    current_x, y_offset + 20,
                    current_x, y_offset + 100,
                    fill='red',
                    width=2
                )
            else:
                # Draw placeholder for non-shared residues
                self.canvas.create_text(
                    x_center, y_offset + 80,
                    text='·',
                    font=('Courier', 14),
                    fill='gray'
                )
            
            current_x += self.char_width
        
        # Draw row separator (except for last row)
        if row < self.num_rows - 1:
            self.canvas.create_line(
                0, y_offset + self.row_height,
                self.residues_per_row * self.char_width + 20, y_offset + self.row_height,
                fill='lightblue',
                width=1
            )
    
    def draw_selections(self):
        """Draw all selected ranges."""
        for chain, start, end in self.ranges:
            # Find positions in sequence
            start_idx = next(i for i, (c, n) in enumerate(zip(self.chain_ids, self.residue_numbers))
                           if c == chain and n == start)
            end_idx = next(i for i, (c, n) in enumerate(zip(self.chain_ids, self.residue_numbers))
                         if c == chain and n == end)
            
            # Calculate row and position for start and end
            start_row = start_idx // self.residues_per_row
            end_row = end_idx // self.residues_per_row
            
            if start_row == end_row:
                # Selection within same row
                start_x = 10 + (start_idx % self.residues_per_row) * self.char_width
                end_x = 10 + (end_idx % self.residues_per_row + 1) * self.char_width
                y_offset = start_row * self.row_height
                
                self.canvas.create_rectangle(
                    start_x, y_offset + 20,
                    end_x, y_offset + 100,
                    fill='yellow',
                    stipple='gray50',
                    outline='orange',
                    width=2
                )
            else:
                # Selection spans multiple rows
                # Draw first row selection
                start_x = 10 + (start_idx % self.residues_per_row) * self.char_width
                y_offset = start_row * self.row_height
                self.canvas.create_rectangle(
                    start_x, y_offset + 20,
                    self.residues_per_row * self.char_width + 10, y_offset + 100,
                    fill='yellow',
                    stipple='gray50',
                    outline='orange',
                    width=2
                )
                
                # Draw middle rows (if any)
                for row in range(start_row + 1, end_row):
                    y_offset = row * self.row_height
                    self.canvas.create_rectangle(
                        10, y_offset + 20,
                        self.residues_per_row * self.char_width + 10, y_offset + 100,
                        fill='yellow',
                        stipple='gray50',
                        outline='orange',
                        width=2
                    )
                
                # Draw last row selection
                end_x = 10 + (end_idx % self.residues_per_row + 1) * self.char_width
                y_offset = end_row * self.row_height
                self.canvas.create_rectangle(
                    10, y_offset + 20,
                    end_x, y_offset + 100,
                    fill='yellow',
                    stipple='gray50',
                    outline='orange',
                    width=2
                )
    
    def get_sequence_index_from_coordinates(self, x, y):
        """Convert mouse coordinates to sequence index."""
        # Adjust for margin
        x = x - 10
        
        # Calculate row from y coordinate
        row = int(y // self.row_height)
        if row < 0 or row >= self.num_rows:
            return None
        
        # Calculate position within row
        col = int(x // self.char_width)
        if col < 0 or col >= self.residues_per_row:
            return None
        
        # Calculate global index
        global_idx = row * self.residues_per_row + col
        
        if global_idx >= len(self.sequence):
            return None
        
        return global_idx
    
    def start_selection(self, event):
        """Handle mouse button press."""
        idx = self.get_sequence_index_from_coordinates(event.x, event.y)
        
        if idx is not None and idx in self.shared_positions:
            self.current_selection = {
                'start_idx': idx,
                'chain': self.chain_ids[idx],
                'start_res': self.residue_numbers[idx]
            }
    
    def update_selection(self, event):
        """Handle mouse drag."""
        if self.current_selection:
            idx = self.get_sequence_index_from_coordinates(event.x, event.y)
            
            if idx is not None and idx in self.shared_positions and self.chain_ids[idx] == self.current_selection['chain']:
                # Draw temporary selection
                self.draw_sequence()  # Redraw to clear previous temp selection
                
                # Calculate selection rectangle
                start_idx = self.current_selection['start_idx']
                start_row = start_idx // self.residues_per_row
                end_row = idx // self.residues_per_row
                
                if start_row == end_row:
                    # Selection within same row
                    start_x = 10 + (start_idx % self.residues_per_row) * self.char_width
                    end_x = 10 + (idx % self.residues_per_row + 1) * self.char_width
                    y_offset = start_row * self.row_height
                    
                    self.canvas.create_rectangle(
                        start_x, y_offset + 20,
                        end_x, y_offset + 100,
                        fill='yellow', stipple='gray50'
                    )
                else:
                    # Multi-row selection - draw full rows
                    for row in range(min(start_row, end_row), max(start_row, end_row) + 1):
                        y_offset = row * self.row_height
                        self.canvas.create_rectangle(
                            10, y_offset + 20,
                            self.residues_per_row * self.char_width + 10, y_offset + 100,
                            fill='yellow', stipple='gray50'
                        )
    
    def end_selection(self, event):
        """Handle mouse button release."""
        if self.current_selection:
            idx = self.get_sequence_index_from_coordinates(event.x, event.y)
            
            if idx is not None and idx in self.shared_positions and self.chain_ids[idx] == self.current_selection['chain']:
                end_res = self.residue_numbers[idx]
                start_res = self.current_selection['start_res']
                chain = self.current_selection['chain']
                
                # Add range (ensure start < end)
                if start_res > end_res:
                    start_res, end_res = end_res, start_res
                self.ranges.append((chain, start_res, end_res))
                
                # Redraw everything
                self.draw_sequence()
            
            self.current_selection = None
    
    def clear_selections(self):
        """Clear all selected ranges."""
        self.ranges = []
        self.draw_sequence()
    
    def finish_selection(self):
        """Complete the selection process."""
        self.quit()

def merge_ranges(ranges):
    """Merge overlapping residue ranges."""
    if not ranges:
        return []
    chain_ranges = {}
    print(f"[DEBUG] Merging ranges: {ranges}")
    for chain, start, end in ranges:
        if chain not in chain_ranges:
            chain_ranges[chain] = []
        chain_ranges[chain].append((start, end))
    merged = []

    for chain, ranges in chain_ranges.items():
        ranges.sort()
        current = list(ranges[0])
        chain_merged = []
        for start, end in ranges[1:]:
            if start <= current[1] + 1:
                current[1] = max(current[1], end)
            else:
                chain_merged.append(tuple(current))
                current = [start, end]
        chain_merged.append(tuple(current))
        merged.extend((chain, start, end) for start, end in chain_merged)
    return merged

def prompt_residue_ranges(sequence_groups, reference_group):
    """Prompt user for residue ranges to use in structural alignment using GUI."""
    shared_positions = find_shared_residues(sequence_groups, reference_group)
    if not shared_positions:
        print("\nError: No residues are shared among all sequence groups.")
        return []
    while True:
        print("\nPlease select residue ranges in the GUI window...")
        print("Note: Only residues shared among all structures are shown and selectable.")
        selector = RangeSelector(reference_group['sequence'],
                               reference_group['residue_numbers'],
                               reference_group['chain_ids'],
                               shared_positions)
        selector.mainloop()
        ranges = selector.ranges
        selector.destroy()
        if not ranges:
            print("\nNo residue ranges selected.")
            response = input("Would you like to try selecting ranges again? (y/n): ").lower()
            if response in ['n', 'no']:
                return []
            continue
        merged_ranges = merge_ranges(ranges)
        print("\nSelected ranges (after merging overlapping regions):")
        for chain, start, end in merged_ranges:
            print(f"Chain {chain}: {start}-{end}")
        while True:
            response = input("\nAre you happy with these ranges? (y/n): ").lower()
            if response in ['y', 'yes']:
                return merged_ranges
            elif response in ['n', 'no']:
                response = input("Would you like to try selecting ranges again? (y/n): ").lower()
                if response in ['n', 'no']:
                    return []
                break
            else:
                print("Please enter 'y' or 'n'")

if __name__ == "__main__":
    main() 