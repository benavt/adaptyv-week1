#!/usr/bin/env python3
"""
Script to create t-SNE visualization from distogram .npy files.
This script reads distogram files for a given PDB ID and creates an interactive 2D t-SNE plot
to visualize the structural relationships between proteins.
"""
import os
import numpy as np
from sklearn.manifold import TSNE
import plotly.express as px
import plotly.graph_objects as go
import argparse
from pathlib import Path
from Bio import PDB
from Bio.PDB import *
from Bio import Align
from Bio.Align import substitution_matrices
import webbrowser
import pandas as pd
import sys
import json
from http.server import HTTPServer, BaseHTTPRequestHandler
import threading
import subprocess
import socket
import time
import pymol
from pymol import cmd
import tkinter as tk
from tkinter import ttk
import math
import glob
from collections import defaultdict
from tqdm import tqdm


# Add the parent directory to the path so we can import the CSPAnalyzer
sys.path.append(str(Path(__file__).parent.parent))
from util import process_model_metrics_file
from structure.csp import CSPAnalyzer
from Sequence import Sequence, SubSequence
from paths import real_CSList_dir, apo_NMR_shift_dir, holo_NMR_shift_dir, PDB_FILES


def convert_aa_name(three_letter_code):
    """Convert three-letter amino acid codes to one-letter codes."""
    aa_dict = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 
        'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 
        'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 
        'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 
        'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    return aa_dict.get(three_letter_code, '?')

# Add PDB file paths to the DataFrame - updated path construction with source-based directories
def get_pdb_path(filepath, pdb_id):
    # Extract the base filename (e.g., 'exp_5vf0_prot_pept' from the full path)
    base_name = Path(filepath).stem
    
    candidate_dirs = []

    # Determine the appropriate directory based on the filename
    if 'comp' in base_name.lower():
        # For computational structures
        candidate_dirs.append(Path('PDB_FILES') / 'computational_structures')
    elif 'exp' in base_name.lower():
        # For experimental structures
        candidate_dirs.append(Path('PDB_FILES') / 'experimental_structures')
    else:
        # For other structures, try both possible locations
        candidate_dirs.append(Path('PDB_FILES') / f"{pdb_id.upper()}_aligned")
        candidate_dirs.append(Path('PDB_FILES') / f"{pdb_id}_haddock_min")
    
    # Recursively remove suffixes (separated by '_') until a matching file is found
    candidate_base = base_name
    while candidate_base:
        for candidate_dir in candidate_dirs:
            candidate_pdb = candidate_dir / (candidate_base + ".pdb")
            if candidate_pdb.exists():
                return str(candidate_pdb)
        # Remove the last suffix (e.g. exp_5vf0_prot_pept -> exp_5vf0_prot)
        candidate_base = '_'.join(candidate_base.split('_')[:-1])
    
    # If no file is found, print a warning and return the original constructed path (using candidate_dir and base_name)
    print(f"Warning: PDB file not found for {base_name} in {candidate_dirs}. Tried candidate bases (recursively removing suffixes) but none exist.")
    return str(candidate_dir / (base_name + ".pdb"))

def get_bmrb_id(pdb_id, csv_file='CSPRANK.csv'):
    """
    Look up the BMRB ID for the apo state from the CSPRANK.csv file.
    
    Args:
        pdb_id: PDB ID of the holo state
        csv_file: Path to the CSPRANK.csv file
        
    Returns:
        BMRB ID for the apo state, or None if not found
    """
    try:
        df = pd.read_csv(csv_file)
        # Find the row where holo_pdb matches the pdb_id
        row = df[df['holo_pdb'] == pdb_id.lower()]
        print(pdb_id.lower())
        print(row)
        print("HERE")
        if not row.empty:
            return row['apo_bmrb'].iloc[0]
        return None
    except Exception as e:
        print(f"Error reading CSPRANK.csv: {e}")
        return None

def align_sequences(seq1: SubSequence, seq2: SubSequence):
    """
    Align two protein sequences using the Sequence class alignment method.
    
    Args:
        seq1: First sequence (SubSequence object)
        seq2: Second sequence (SubSequence object)
        
    Returns:
        Tuple containing:
        - Aligned first sequence
        - Aligned second sequence
        - Dictionary mapping positions in seq2 (CSP) to positions in seq1 (PDB)
    """
    # Use Sequence class alignment method
    aligned_seq1, aligned_seq2, aligned_values1, aligned_values2, score = Sequence.align_subsequences(seq1, seq2)
    
    # Create mapping from CSP positions (seq2) to PDB positions (seq1)
    mapping = {}
    pdb_pos = 0
    csp_pos = 0
    
    for i in range(len(aligned_seq1)):
        if aligned_seq1[i] != '_' and aligned_seq2[i] != '_':
            # Both sequences have a residue at this position
            mapping[csp_pos] = pdb_pos
            pdb_pos += 1
            csp_pos += 1
        elif aligned_seq1[i] != '_':
            # Only PDB sequence has a residue
            pdb_pos += 1
        #elif aligned_seq2[i] != '_':
        else:
            # Only CSP sequence has a residue
            csp_pos += 1
    
    print("\nDebug - Alignment mapping:")
    print(f"Aligned PDB sequence:  {aligned_seq1}")
    print(f"Aligned CSP sequence:  {aligned_seq2}")
    print(f"Mapping (CSP -> PDB): {mapping}")
    
    return aligned_seq1, aligned_seq2, mapping

def calculate_csps(pdb_id, bmrb_id, method):
    """
    Calculate CSPs for a given PDB ID and BMRB ID.
    
    Args:
        pdb_id: PDB ID of the holo state
        bmrb_id: BMRB ID of the apo state
        method: CSP calculation method
        
    Returns:
        Tuple of (csps, cutoff, sequence) or (None, None, None) if calculation fails
    """
    # Construct paths to chemical shift files
    apo_file = os.path.join(real_CSList_dir, f"{bmrb_id}.csv")
    holo_file = os.path.join(real_CSList_dir, f"{pdb_id.upper()}.csv")
    
    # Check if files exist
    if not os.path.exists(apo_file):
        print(f"Error: Apo chemical shift file not found: {apo_file}")
        return None, None, None
    
    if not os.path.exists(holo_file):
        print(f"Error: Holo chemical shift file not found: {holo_file}")
        return None, None, None
    
    try:
        # Calculate CSPs
        analyzer = CSPAnalyzer(method=method)
        csps, cutoff, sequence = analyzer.calc_csp(apo_file, holo_file)
        return csps, cutoff, sequence
    except Exception as e:
        print(f"Error calculating CSPs for {pdb_id}: {e}")
        return None, None, None
    
def get_longer_chain(pdb_file):
    """Get the longer chain from a PDB file.
    
    Args:
        pdb_file: Path to the PDB file
        
    Returns:
        Tuple of (protein_chain_id, peptide_chain_id)
    """
    try:
        # Parse the PDB file
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure('protein', pdb_file)
        
        # Get chain lengths
        chain_lengths = {}
        for model in structure:
            for chain in model:
                chain_lengths[chain.id] = len(list(chain.get_residues()))
        
        # Find the longer chain
        if not chain_lengths:
            print("Error: No chains found in PDB file")
            return None, None
        
        # Sort chains by length
        sorted_chains = sorted(chain_lengths.items(), key=lambda x: x[1], reverse=True)
        protein_chain = sorted_chains[0][0]
        peptide_chain = sorted_chains[1][0] if len(sorted_chains) > 1 else None
        
        print(f"Found protein chain {protein_chain} ({chain_lengths[protein_chain]} residues)")
        if peptide_chain:
            print(f"Found peptide chain {peptide_chain} ({chain_lengths[peptide_chain]} residues)")
        
        return protein_chain, peptide_chain
    
    except Exception as e:
        print(f"Error getting chain lengths from PDB file {pdb_file}: {e}")
        return None, None
    
def get_pdb_sequence(pdb_file, chain_id):
    """Get the sequence from a specific chain in a PDB file.
    
    Args:
        pdb_file: Path to the PDB file
        chain_id: Chain ID to get sequence from
        
    Returns:
        SubSequence object containing the sequence
    """
    try:
        # Parse the PDB file
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure('protein', pdb_file)
        
        # Get specified chain
        target_chain = None
        for model in structure:
            for chain in model:
                if chain.id == chain_id:
                    target_chain = chain
                    break
            if target_chain:
                break
        
        if target_chain is None:
            print(f"Warning: Could not find chain {chain_id} in {pdb_file}")
            return None
        
        # Get sequence from chain
        sequence = []
        residue_numbers = []
        for res in target_chain:
            if 'CA' in res:  # Only consider residues with CA atoms
                residue_name = res.get_resname()
                aa_code = convert_aa_name(residue_name)
                sequence.append(aa_code)
                residue_numbers.append(res.get_id()[1])  # Get residue number
        
        # Create SubSequence object
        subsequence = SubSequence(
            sequence=''.join(sequence),
            chain_id=chain_id,
            is_protein=True,
            is_peptide=False,
            residue_numbers=residue_numbers
        )
        
        return subsequence
    
    except Exception as e:
        print(f"Error getting sequence from PDB file {pdb_file}: {e}")
        return None

def get_pdb_sequence_full(pdb_path):
    """Extract sequence and residue numbers from PDB file (from align_files.py)."""
    sequence = []
    residue_numbers = []
    chain_ids = []
    current_chain = ""
    prev_chain = ""
    prev_res_id = None

    with open(pdb_path, "r") as pdb_file:
        for line in pdb_file:
            if line.startswith("ATOM"):
                chain = line[21]
                res_name = line[17:20].strip()
                res_id = int(line[22:26])
                
                if prev_chain != "" and prev_chain != chain:
                    sequence.append(':')
                    residue_numbers.append(-1)
                    chain_ids.append(' ')
                if prev_res_id != res_id:
                    sequence.append(convert_aa_name(res_name))
                    residue_numbers.append(res_id)
                    chain_ids.append(chain)

                prev_chain = chain
                prev_res_id = res_id

    formatted_sequence = "".join(sequence)
    formatted_residue_numbers = format_residue_numbers(residue_numbers, chain_ids)

    return {
        'sequence': formatted_sequence,
        'residue_numbers': residue_numbers,
        'chain_ids': chain_ids,
        'formatted_numbers': formatted_residue_numbers
    }

def format_residue_numbers(residue_numbers, chain_ids):
    """Format residue numbers for display."""
    if not residue_numbers:
        return ""
        
    num_lines = math.ceil(math.log(max(residue_numbers), 10) + 1)
    positions = [[" "] * len(residue_numbers) for i in range(0, num_lines)]
    
    for line in range(0, num_lines):
        if line == 0:
            for col, chain_id in enumerate(chain_ids):
                positions[line][col] = chain_id
        else:
            for col, residue_number in enumerate(residue_numbers):
                if residue_number == -1:
                    continue
                t_res = residue_number - (residue_number // int(math.pow(10, num_lines-line))) * int(math.pow(10, num_lines-line))
                positions[line][col] = str(t_res // int(math.pow(10, num_lines - line - 1)))

    return '\n'.join([''.join(position) for position in positions])

class SequenceGroup:
    """Class to store sequence information and related files."""
    def __init__(self, sequence_info, file_path):
        self.sequence = sequence_info['sequence']
        self.residue_numbers = sequence_info['residue_numbers']
        self.chain_ids = sequence_info['chain_ids']
        self.formatted_numbers = sequence_info['formatted_numbers']
        self.files = [file_path]
        
    def matches(self, other_sequence_info):
        """Check if another sequence matches this group's sequence."""
        return self.sequence == other_sequence_info['sequence']
        
    def add_file(self, file_path):
        """Add a file to this sequence group."""
        if file_path not in self.files:
            self.files.append(file_path)

def find_unique_sequences(directory_path):
    """Find all unique sequences in PDB files in the directory."""
    if not os.path.isdir(directory_path):
        print(f"Directory not found: {directory_path}")
        return []

    pdb_files = glob.glob(os.path.join(directory_path, "*.pdb"))
    if not pdb_files:
        print(f"No PDB files found in {directory_path}")
        return []

    print(f"Found {len(pdb_files)} PDB files to process")
    
    sequence_groups = []
    
    for file_path in tqdm(pdb_files):
        try:
            sequence_info = get_pdb_sequence_full(file_path)
            
            # Check if sequence matches any existing group
            matched = False
            for group in sequence_groups:
                if group.matches(sequence_info):
                    group.add_file(file_path)
                    matched = True
                    break
                    
            # If no match found, create new group
            if not matched:
                sequence_groups.append(SequenceGroup(sequence_info, file_path))
                
        except Exception as e:
            print(f"Error processing {file_path}: {str(e)}")
            continue
    
    return sequence_groups

def find_shared_residues(sequence_groups, reference_group):
    """Find residues that are shared among all sequence groups.
    
    Args:
        sequence_groups: List of SequenceGroup objects
        reference_group: Reference SequenceGroup object
        
    Returns:
        List of indices in the reference sequence that are shared among all groups
    """
    # Start with all residues from reference group
    shared_positions = []
    
    # Create a list of (chain, residue_number) pairs for each position in reference
    ref_positions = [(c, n) for c, n in zip(reference_group.chain_ids, reference_group.residue_numbers)
                    if c != ' ' and n != -1]  # Exclude chain breaks
    
    # Check each position against all other groups
    for i, (ref_chain, ref_res) in enumerate(ref_positions):
        is_shared = True
        for group in sequence_groups:
            if group == reference_group:
                continue
                
            # Find if this residue exists in the other group
            found = False
            for chain, res_num in zip(group.chain_ids, group.residue_numbers):
                if chain == ref_chain and res_num == ref_res:
                    found = True
                    break
            
            if not found:
                is_shared = False
                break
        
        if is_shared:
            # Find the index in the reference sequence
            ref_idx = next(i for i, (c, n) in enumerate(zip(reference_group.chain_ids, reference_group.residue_numbers))
                         if c == ref_chain and n == ref_res)
            shared_positions.append(ref_idx)
    
    return sorted(shared_positions)

def merge_ranges(ranges):
    """Merge overlapping residue ranges.
    
    Args:
        ranges: List of (chain, start, end) tuples
    
    Returns:
        List of merged (chain, start, end) tuples
    """
    if not ranges:
        return []
    
    # Group ranges by chain
    chain_ranges = {}
    for chain, start, end in ranges:
        if chain not in chain_ranges:
            chain_ranges[chain] = []
        chain_ranges[chain].append((start, end))
    
    # Merge ranges for each chain
    merged = []
    for chain, ranges in chain_ranges.items():
        # Sort ranges by start position
        ranges.sort()
        
        # Merge overlapping ranges
        current = list(ranges[0])  # Convert to list for mutability
        chain_merged = []
        
        for start, end in ranges[1:]:
            if start <= current[1] + 1:  # +1 to merge adjacent ranges
                current[1] = max(current[1], end)
            else:
                chain_merged.append(tuple(current))
                current = [start, end]
        chain_merged.append(tuple(current))
        
        # Add merged ranges with chain ID
        merged.extend((chain, start, end) for start, end in chain_merged)
    
    return merged

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

def prompt_residue_ranges(pdb_dir):
    """Prompt user for residue ranges to use in structural alignment using GUI."""
    # Find unique sequences in the PDB directory
    sequence_groups = find_unique_sequences(pdb_dir)
    
    if not sequence_groups:
        print(f"\nError: No sequences found in {pdb_dir}")
        return []
    
    if len(sequence_groups) == 1:
        # Only one sequence group, use it as reference
        reference_group = sequence_groups[0]
    else:
        # Multiple sequence groups, let user choose
        print("\nAvailable sequence groups:")
        for i, group in enumerate(sequence_groups, 1):
            print(f"\nGroup {i} ({len(group.files)} files):")
            print(group.formatted_numbers)
            print(group.sequence)
            print("-" * len(group.sequence))
        
        while True:
            try:
                selection = int(input(f"\nEnter the number of the sequence group to use as reference (1-{len(sequence_groups)}): "))
                if 1 <= selection <= len(sequence_groups):
                    reference_group = sequence_groups[selection - 1]
                    break
                print("Invalid selection. Please try again.")
            except ValueError:
                print("Please enter a valid number.")
    
    # Find shared residues
    shared_positions = find_shared_residues(sequence_groups, reference_group)
    
    if not shared_positions:
        print("\nError: No residues are shared among all sequence groups.")
        return []
    
    # Use GUI for residue selection (now always called from main thread)
    return prompt_residue_ranges_gui(reference_group, shared_positions)

def prompt_residue_ranges_gui(reference_group, shared_positions):
    """Prompt user for residue ranges using GUI (main thread only)."""
    while True:
        print("\nPlease select residue ranges in the GUI window...")
        print("Note: Only residues shared among all structures are shown and selectable.")
        
        # Create and run GUI
        selector = RangeSelector(reference_group.sequence,
                               reference_group.residue_numbers,
                               reference_group.chain_ids,
                               shared_positions)
        selector.mainloop()
        
        # Get selected ranges and destroy window
        ranges = selector.ranges
        selector.destroy()
        
        if not ranges:
            print("\nNo residue ranges selected.")
            response = input("Would you like to try selecting ranges again? (y/n): ").lower()
            if response in ['n', 'no']:
                return []
            continue
        
        # Merge overlapping ranges
        merged_ranges = merge_ranges(ranges)
        
        print("\nSelected ranges (after merging overlapping regions):")
        for chain, start, end in merged_ranges:
            print(f"Chain {chain}: {start}-{end}")
        
        # Ask for confirmation
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

def prompt_residue_ranges_cli(reference_group, shared_positions):
    """Prompt user for residue ranges using command-line interface."""
    print("\nCommand-line residue selection interface")
    print("Available shared residues:")
    
    # Group residues by chain
    chain_residues = {}
    for pos in shared_positions:
        chain = reference_group.chain_ids[pos]
        res_num = reference_group.residue_numbers[pos]
        aa = reference_group.sequence[pos]
        if chain not in chain_residues:
            chain_residues[chain] = []
        chain_residues[chain].append((res_num, aa))
    
    # Display available residues by chain
    for chain in sorted(chain_residues.keys()):
        residues = sorted(chain_residues[chain])
        print(f"\nChain {chain}:")
        for res_num, aa in residues:
            print(f"  {res_num}: {aa}")
    
    ranges = []
    print("\nEnter residue ranges for alignment (format: chain start-end, e.g., A 1-50)")
    print("Enter 'done' when finished, or 'skip' to use default alignment")
    
    while True:
        try:
            user_input = input("\nEnter range (or 'done'/'skip'): ").strip()
            
            if user_input.lower() == 'done':
                break
            elif user_input.lower() == 'skip':
                return []
            elif not user_input:
                continue
            
            # Parse input
            parts = user_input.split()
            if len(parts) != 2:
                print("Invalid format. Use: chain start-end (e.g., A 1-50)")
                continue
            
            chain = parts[0]
            range_part = parts[1]
            
            if '-' not in range_part:
                print("Invalid range format. Use: start-end (e.g., 1-50)")
                continue
            
            try:
                start, end = map(int, range_part.split('-'))
                if start > end:
                    start, end = end, start
            except ValueError:
                print("Invalid range numbers. Use: start-end (e.g., 1-50)")
                continue
            
            # Validate that the range contains shared residues
            valid_residues = []
            for res_num, aa in chain_residues.get(chain, []):
                if start <= res_num <= end:
                    valid_residues.append(res_num)
            
            if not valid_residues:
                print(f"No shared residues found in range {chain} {start}-{end}")
                continue
            
            ranges.append((chain, start, end))
            print(f"Added range: {chain} {start}-{end} (contains {len(valid_residues)} shared residues)")
            
        except KeyboardInterrupt:
            print("\nCancelled by user.")
            return []
        except EOFError:
            print("\nCancelled.")
            return []
    
    if not ranges:
        print("No ranges selected.")
        return []
    
    # Merge overlapping ranges
    merged_ranges = merge_ranges(ranges)
    
    print("\nSelected ranges (after merging overlapping regions):")
    for chain, start, end in merged_ranges:
        print(f"Chain {chain}: {start}-{end}")
    
    # Ask for confirmation
    while True:
        response = input("\nAre you happy with these ranges? (y/n): ").lower()
        if response in ['y', 'yes']:
            return merged_ranges
        elif response in ['n', 'no']:
            response = input("Would you like to try selecting ranges again? (y/n): ").lower()
            if response in ['n', 'no']:
                return []
            # Reset and try again
            ranges = []
            print("\nEnter residue ranges for alignment (format: chain start-end, e.g., A 1-50)")
            print("Enter 'done' when finished, or 'skip' to use default alignment")
            break
        else:
            print("Please enter 'y' or 'n'")

def color_structure_by_csp(obj_name, protein_chain, peptide_chain, csps=None, cutoff=None, sequence=None, pdb_subsequence=None):
    """
    Color a structure in PyMOL based on CSP significance.
    
    Args:
        obj_name: Name of the PyMOL object
        protein_chain: ID of the protein chain
        peptide_chain: ID of the peptide chain (if any)
        csps: CSP values array
        cutoff: CSP significance cutoff
        sequence: Sequence object containing CSP data
        pdb_subsequence: SubSequence object for the PDB structure
    """
    # Basic visualization setup
    cmd.hide('everything', obj_name)
    cmd.show('cartoon', obj_name)
    
    # Color the protein chain gray by default
    cmd.color('gray', f'{obj_name} and chain {protein_chain}')
    
    # Color the peptide chain yellow if it exists
    if peptide_chain:
        cmd.color('yellow', f'{obj_name} and chain {peptide_chain}')
    
    # If CSP data is available, color based on significance
    if csps is not None and sequence is not None and pdb_subsequence is not None:
        # Align the sequences
        aligned_pdb_seq, aligned_csp_seq, mapping = align_sequences(pdb_subsequence, sequence.subsequences[0])
        
        # Calculate z-scores
        csps_below_cutoff = [csp for csp in csps if csp <= cutoff and csp > 0]
        mean = np.mean(csps_below_cutoff)
        std = np.std(csps_below_cutoff)
        
        # Find residues with significant CSPs and their z-scores
        significant_residues = []
        z_scores = {}
        
        # Map CSP positions to PDB positions
        for i, res in enumerate(sequence.subsequences[0].sequence):
            if i in mapping:
                pdb_pos = mapping[i]
                csp = csps[i]
                z_score = (csp - mean) / std if std > 0 else 0
                z_scores[pdb_pos] = z_score
                if z_score > 0:
                    significant_residues.append(pdb_pos)
        
        # Create selections for different z-score ranges
        light_red_residues = []
        red_residues = []
        deep_red_residues = []
        
        for res_num in significant_residues:
            z_score = z_scores[res_num]
            if 0 <= z_score < 1:
                light_red_residues.append(res_num)
            elif 1 <= z_score < 2:
                red_residues.append(res_num)
            else:  # z_score >= 2
                deep_red_residues.append(res_num)
        
        # Color residues based on z-score ranges
        if light_red_residues:
            cmd.color('salmon', f'{obj_name} and chain {protein_chain} and resi {"+".join(map(str, light_red_residues))}')
        if red_residues:
            cmd.color('red', f'{obj_name} and chain {protein_chain} and resi {"+".join(map(str, red_residues))}')
        if deep_red_residues:
            cmd.color('ruby', f'{obj_name} and chain {protein_chain} and resi {"+".join(map(str, deep_red_residues))}')

class PyMOLServer(BaseHTTPRequestHandler):
    # Class variables to track PyMOL state
    pymol_initialized = False
    pymol_window_launched = False
    first_structure_loaded = False
    pymol_process = None  # Track the PyMOL process
    first_obj_name = None
    first_chain = None
    csp_data = None  # Store CSP data for reuse
    peptide_colors = [  # Define a set of distinct colors for peptides
        'yellow', 'orange', 'red', 'pink', 'purple', 'blue', 'cyan', 'green',
        'lime', 'salmon', 'violet', 'magenta', 'marine', 'teal', 'forest', 'olive'
    ]
    current_color_index = 0  # Track which color to use next
    # New variables for residue selection
    selected_residue_ranges = None  # Store selected residue ranges for alignment
    pdb_dir = None  # Store PDB directory for residue selection

    def get_next_peptide_color(self):
        """Get the next color for a peptide chain, cycling through the color list."""
        color = PyMOLServer.peptide_colors[PyMOLServer.current_color_index]
        PyMOLServer.current_color_index = (PyMOLServer.current_color_index + 1) % len(PyMOLServer.peptide_colors)
        print(f"[PyMOL Server] Current color index: {PyMOLServer.current_color_index}, using color: {color}")
        return color

    def do_OPTIONS(self):
        """Handle OPTIONS request for CORS preflight"""
        self.send_response(200)
        self.send_header('Access-Control-Allow-Origin', '*')
        self.send_header('Access-Control-Allow-Methods', 'POST, OPTIONS')
        self.send_header('Access-Control-Allow-Headers', 'Content-Type')
        self.end_headers()
    
    def send_json_response(self, status_code, data):
        """Helper method to send JSON response with proper headers"""
        self.send_response(status_code)
        self.send_header('Access-Control-Allow-Origin', '*')
        self.send_header('Content-type', 'application/json')
        self.end_headers()
        self.wfile.write(json.dumps(data).encode('utf-8'))
    

    def launch_pymol(self, pdb_path):
        # If a PyMOL window is already open, close it first
        if PyMOLServer.pymol_process is not None:
            try:
                PyMOLServer.pymol_process.terminate()
                PyMOLServer.pymol_process.wait(timeout=5)
                print("[PyMOL Server] Closed previous PyMOL window")
            except Exception as e:
                print(f"[PyMOL Server] Error closing previous PyMOL window: {e}")

        # Launch GUI after loading first structure
        try:
            # Save current session to a temporary file
            temp_session = f"{os.path.splitext(pdb_path)[0]}_temp.pse"
            cmd.save(temp_session)
            
            # Launch PyMOL GUI with the session
            import subprocess
            PyMOLServer.pymol_process = subprocess.Popen(['pymol', temp_session])
            
            # Wait a moment for PyMOL to open
            import time
            time.sleep(2)
            
            # Remove the temporary session file
            os.remove(temp_session)
            
            print("[PyMOL Server] PyMOL GUI launched")
            PyMOLServer.pymol_window_launched = True
        except Exception as e:
            error_msg = f"Error launching PyMOL GUI: {str(e)}"
            print(f"[PyMOL Server] {error_msg}")
            self.send_json_response(500, {
                'status': 'error',
                'message': error_msg
            })
            return

    def clear_pymol_cache(self):
        """Clear all loaded structures in PyMOL and reinitialize for fresh start."""
        try:
            print("[PyMOL Server] Clearing PyMOL cache and reinitializing...")
            
            # Close PyMOL process if it exists
            if PyMOLServer.pymol_process is not None:
                try:
                    PyMOLServer.pymol_process.terminate()
                    PyMOLServer.pymol_process.wait(timeout=5)
                    print("[PyMOL Server] Closed PyMOL window")
                except Exception as e:
                    print(f"[PyMOL Server] Error closing PyMOL window: {e}")
                finally:
                    PyMOLServer.pymol_process = None
            # Clear all objects in PyMOL if it's initialized
            if PyMOLServer.pymol_initialized:
                try:
                    cmd.delete('all')
                    print("[PyMOL Server] Cleared all PyMOL objects")
                except Exception as e:
                    print(f"[PyMOL Server] Error clearing PyMOL objects: {e}")

            # Reset PyMOL state
            PyMOLServer.pymol_initialized = False
            PyMOLServer.pymol_window_launched = False
            PyMOLServer.first_structure_loaded = False
            PyMOLServer.first_obj_name = None
            PyMOLServer.first_chain = None
            PyMOLServer.current_color_index = 0
            
            print("[PyMOL Server] PyMOL cache cleared and ready for fresh start")
            return True
            
        except Exception as e:
            error_msg = f"Error clearing PyMOL cache: {str(e)}"
            print(f"[PyMOL Server] {error_msg}")
            return False

    def do_POST(self):
        """Handle POST request with CORS headers"""
        try:
            content_length = int(self.headers['Content-Length'])
            post_data = self.rfile.read(content_length)
            data = json.loads(post_data.decode('utf-8'))
            
            # Handle logging requests
            if self.path == '/log':
                print("\n[DEBUG LOG]", data.get('message', 'No message'))
                if 'data' in data:
                    print("[DEBUG DATA]", json.dumps(data['data'], indent=2))
                self.send_json_response(200, {'status': 'logged'})
                return
            
            # Handle clear cache requests
            if self.path == '/clear_cache':
                success = self.clear_pymol_cache()
                if success:
                    self.send_json_response(200, {'status': 'success', 'message': 'PyMOL cache cleared'})
                else:
                    self.send_json_response(500, {'status': 'error', 'message': 'Failed to clear PyMOL cache'})
                return
            
            # Handle normal PDB loading requests
            requested_path = data.get('pdb_path')
            highlight_color = data.get('highlight_color')  # Get the highlight color from the request
            
            print(f"\n[PyMOL Server] Received request to load: {requested_path} with color: {highlight_color}")
            
            if not requested_path:
                self.send_json_response(400, {
                    'status': 'error',
                    'message': 'No PDB path provided'
                })
                return
            

            # Extract the basename (without extension) from the requested path
            requested_basename = os.path.splitext(os.path.basename(requested_path))[0]
            for i in range(0, len(requested_basename.split('_'))):
                pdb_dir = self.server.pdb_dir
                # Find a .pdb file in pdb_dir with the same basename
                matched_pdb = None
                for fname in os.listdir(pdb_dir):
                    if fname.lower().endswith('.pdb') and os.path.splitext(fname)[0] == requested_basename:
                        matched_pdb = os.path.join(pdb_dir, fname)
                        break
                if matched_pdb:
                    break
                requested_basename = '_'.join(requested_basename.split('_')[:-1])
                
            if not matched_pdb:
                self.send_json_response(404, {
                    'status': 'error',
                    'message': f'PDB file not found for basename {requested_basename} in {pdb_dir}'
                })
                return
            pdb_path = matched_pdb
            print(f"[PyMOL Server] Loading structure in PyMOL: {pdb_path}")
            
            # Initialize PyMOL and launch GUI on first click
            if not PyMOLServer.pymol_initialized:
                print("[PyMOL Server] Initializing PyMOL for the first time")
                pymol.finish_launching(['pymol', '-c'])
                PyMOLServer.pymol_initialized = True
                print("[PyMOL Server] PyMOL initialized in command-line mode")
            
            # Get a unique object name based on the PDB filename
            obj_name = os.path.splitext(os.path.basename(pdb_path))[0]
            
            # Check if object already exists
            if obj_name in cmd.get_names('objects'):
                print(f"[PyMOL Server] Object {obj_name} already exists, reusing it")
                # Just update the view
                cmd.zoom(obj_name)
                self.launch_pymol(pdb_path)
                self.send_json_response(200, {'status': 'success'})
                return
            
            # Load the structure
            cmd.load(pdb_path, obj_name)

            # Color the entire model by the highlight color (from t-SNE)
            cmd.color(highlight_color, obj_name)
            print(f"[PyMOL Server] Colored model {obj_name} with {highlight_color}")
         
            
            if not PyMOLServer.first_structure_loaded:
                # First structure - store reference for alignment
                PyMOLServer.first_obj_name = obj_name
                PyMOLServer.pdb_dir = self.server.pdb_dir  # Store PDB directory for later use
                # Store the pre-selected residue ranges from the server
                PyMOLServer.selected_residue_ranges = self.server.selected_residue_ranges
                cmd.orient()
                cmd.zoom('all', 1.5)
                self.launch_pymol(pdb_path)
                PyMOLServer.first_structure_loaded = True
            else:
                # Subsequent structures - use pre-selected residue ranges for alignment
                if PyMOLServer.selected_residue_ranges:
                    print(f"[PyMOL Server] Using pre-selected residue ranges for alignment:")
                    for chain, start, end in PyMOLServer.selected_residue_ranges:
                        print(f"  Chain {chain}: {start}-{end}")
                    
                    # Create PyMOL selection strings for the selected residues
                    ref_selection_parts = []
                    mobile_selection_parts = []
                    
                    for chain, start, end in PyMOLServer.selected_residue_ranges:
                        ref_selection_parts.append(f"{PyMOLServer.first_obj_name} and chain {chain} and resi {start}-{end}")
                        mobile_selection_parts.append(f"{obj_name} and chain {chain} and resi {start}-{end}")
                    
                    ref_selection = " or ".join(ref_selection_parts)
                    mobile_selection = " or ".join(mobile_selection_parts)
                    
                    print(f"[PyMOL Server] Aligning {mobile_selection} to {ref_selection}")
                    try:
                        cmd.align(mobile_selection, ref_selection)
                    except Exception as e:
                        print(f"[PyMOL Server] Error during alignment: {e}")
                        print("[PyMOL Server] Falling back to default alignment.")
                        cmd.align(obj_name, PyMOLServer.first_obj_name)
                else:
                    print("[PyMOL Server] No residue ranges selected. Using default alignment.")
                    # Use default alignment
                    if PyMOLServer.first_obj_name and obj_name != PyMOLServer.first_obj_name:
                        ref_selection = f"{PyMOLServer.first_obj_name}"
                        mobile_selection = f"{obj_name}"
                        print(f"[PyMOL Server] Aligning {mobile_selection} to {ref_selection}")
                        cmd.align(mobile_selection, ref_selection)
                
                cmd.zoom(obj_name)
                self.launch_pymol(pdb_path)
            
            print("[PyMOL Server] Structure loaded successfully")
            self.send_json_response(200, {'status': 'success'})
            
        except json.JSONDecodeError as e:
            print(f"[PyMOL Server] Invalid JSON in request: {e}")
            self.send_json_response(400, {
                'status': 'error',
                'message': 'Invalid JSON in request'
            })
        except Exception as e:
            error_msg = str(e)
            print(f"[PyMOL Server] Error loading structure: {error_msg}")
            self.send_json_response(500, {
                'status': 'error',
                'message': f"Error loading structure in PyMOL: {error_msg}"
            })

    def get_csp_data(self, pdb_id):
        """Get CSP data for a given PDB ID."""
        try:
            # Get the BMRB ID for the apo state
            bmrb_id = get_bmrb_id(pdb_id, 'CSPRANK.csv')
            if bmrb_id is None:
                print(f"Error: Could not find BMRB ID for PDB ID {pdb_id}")
                return None, None, None

            # Construct paths to chemical shift files
            apo_file = os.path.join(real_CSList_dir, f"{bmrb_id}.csv")
            holo_file = os.path.join(real_CSList_dir, f"{pdb_id.upper()}.csv")

            # Check if files exist
            if not os.path.exists(apo_file):
                print(f"Error: Apo chemical shift file not found: {apo_file}")
                return None, None, None
            if not os.path.exists(holo_file):
                print(f"Error: Holo chemical shift file not found: {holo_file}")
                return None, None, None

            # Calculate CSPs using CSPAnalyzer
            analyzer = CSPAnalyzer(method='MONTE')
            csps, cutoff, sequence = analyzer.calc_csp(apo_file, holo_file)
            return csps, cutoff, sequence
        except Exception as e:
            print(f"Error getting CSP data: {e}")
            return None, None, None

class PyMOLHTTPServer(HTTPServer):
    """Custom HTTPServer class that can store the PDB ID, PDB directory, and selected residue ranges."""
    def __init__(self, server_address, RequestHandlerClass, pdb_id, pdb_dir, selected_residue_ranges=None):
        super().__init__(server_address, RequestHandlerClass)
        self.pdb_id = pdb_id
        self.pdb_dir = pdb_dir
        self.selected_residue_ranges = selected_residue_ranges

def find_free_port():
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(('', 0))
        return s.getsockname()[1]

def start_pymol_server(pdb_dir, selected_residue_ranges=None):
    port = find_free_port()
    pdb_id = pdb_dir.split('/')[-1].upper()
    server = PyMOLHTTPServer(('localhost', port), PyMOLServer, pdb_id, pdb_dir, selected_residue_ranges)
    server_thread = threading.Thread(target=server.serve_forever)
    server_thread.daemon = False
    server_thread.start()
    print(f"[PyMOL Server] Server started on port {port}")
    return port, server, server_thread

def find_distogram_dirs(pdb_id):
    """
    Find all directories matching the pattern ./{pdb_id.upper()}/*_distograms
    
    Args:
        pdb_id: PDB ID to process
        
    Returns:
        List of Path objects for matching directories
    """
    base_dir = Path(f"./{pdb_id.upper()}")
    if not base_dir.exists():
        raise FileNotFoundError(f"Base directory not found: {base_dir}")
        
    distogram_dirs = list(base_dir.glob("*_distograms"))
    return distogram_dirs

def select_distogram_dir(pdb_id, distogram_dirs):
    """
    Prompt user to select a distogram directory from the list of available directories.
    
    Args:
        pdb_id: PDB ID being processed
        distogram_dirs: List of Path objects for available distogram directories
        
    Returns:
        Selected Path object
    """
    if not distogram_dirs:
        raise FileNotFoundError(f"No distogram directories found in ./{pdb_id.upper()}/")
        
    print("\nAvailable distogram directories:")
    for i, dir_path in enumerate(distogram_dirs, 1):
        print(f"{i}. {dir_path.name}")
    
    while True:
        try:
            choice = input("\nSelect a directory number (or 'q' to quit): ")
            if choice.lower() == 'q':
                sys.exit(0)
                
            idx = int(choice) - 1
            if 0 <= idx < len(distogram_dirs):
                return distogram_dirs[idx]
            else:
                print(f"Please enter a number between 1 and {len(distogram_dirs)}")
        except ValueError:
            print("Please enter a valid number")

def load_distograms(input_dir):
    """
    Load all distogram .npy files from the input directory.
    
    Args:
        input_dir: Path to the directory containing distogram .npy files
        
    Returns:
        Tuple of (distograms, filenames, filepaths, distogram_dir_name) where distograms is a list of numpy arrays,
        filenames is a list of corresponding filenames, filepaths is a list of full paths, and distogram_dir_name is the name of the directory used.
    """
    input_dir = Path(input_dir)
    if not input_dir.exists():
        raise FileNotFoundError(f"Directory not found: {input_dir}")
    
    print(f"\nLoading distograms from: {input_dir}")
    
    distograms = []
    filenames = []
    filepaths = []
    
    # Dictionary to store shape counts
    shape_counts = {}
    
    for npy_file in input_dir.glob('*.npy'):
        try:
            distogram = np.load(str(npy_file))
            shape = distogram.shape
            shape_counts[shape] = shape_counts.get(shape, 0) + 1
            
            distograms.append(distogram)
            filenames.append(npy_file.stem)
            filepaths.append(str(npy_file))
        except Exception as e:
            print(f"Error loading {npy_file}: {e}")
    
    # Print shape statistics
    print("\nDistogram shape statistics:")
    print("---------------------------")
    for shape, count in shape_counts.items():
        print(f"Shape {shape}: {count} files")
    print("---------------------------\n")
    
    return distograms, filenames, filepaths, input_dir.name

def load_csp_results(pdb_id):
    """
    Load CSP results from CSV file.
    
    Args:
        pdb_id: PDB ID to process
        
    Returns:
        DataFrame containing CSP results
    """
    results_file = Path(f"./{pdb_id.upper()}/csp_results.csv")
    if not results_file.exists():
        print(f"[DEBUG] CSP results file not found: {results_file}")
        return None
    
    return pd.read_csv(results_file)

def calculate_f1_score(row):
    """
    Calculate F1 score from TP, TN, and FP values.
    F1 = 2 * (precision * recall) / (precision + recall)
    where precision = TP / (TP + FP) and recall = TP / (TP + FN)
    """
    tp = row['TP']
    fp = row['FP']
    tn = row['TN']
    fn = row['FN'] if 'FN' in row else (tn - fp)  # Estimate FN if not present
    
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    
    f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
    return f1

def find_matching_basename(basename, csp_results):
    """
    Find a matching basename in the CSP results by iteratively removing suffixes.
    
    Args:
        basename: The original basename from the distogram file
        csp_results: DataFrame containing CSP results
        
    Returns:
        The matching basename if found, or the original basename if no match found
    """
    # Try the original basename first
    if basename in csp_results['pdb_basename'].values:
        return basename
        
    # Split the basename by '_' and try progressively shorter versions
    parts = basename.split('_')
    for i in range(len(parts)-1, 0, -1):
        test_basename = '_'.join(parts[:i])
        if test_basename in csp_results['pdb_basename'].values:
            # print(f"[DEBUG] Found matching basename: {test_basename} (original: {basename})")
            return test_basename
            
    print(f"[WARNING] No matching basename found for {basename}, using original")
    return basename

def save_tsne_results(input_dir, embeddings, filenames, filepaths, perplexity, distogram_dir_name):
    """Save t-SNE results to a file.
    
    Args:
        input_dir: Input directory being processed
        embeddings: t-SNE embeddings array
        filenames: List of filenames
        filepaths: List of filepaths
        perplexity: t-SNE perplexity parameter used
        distogram_dir_name: Name of the directory used for distograms
    """
    results_dir = Path(f"./{input_dir}/tsne_results")
    results_dir.mkdir(exist_ok=True)
    
    # Sanitize directory name for filename usage
    safe_dir_name = distogram_dir_name.replace('/', '_').replace(' ', '_')
    # Create a unique filename based on perplexity and distogram directory
    results_file = results_dir / f"tsne_{safe_dir_name}_perplexity_{perplexity}.npz"
    
    # Save the results
    np.savez(
        results_file,
        embeddings=embeddings,
        filenames=filenames,
        filepaths=filepaths,
        perplexity=perplexity,
        distogram_dir_name=distogram_dir_name
    )
    print(f"\nSaved t-SNE results to: {results_file}")

def load_tsne_results(input_dir, perplexity, distogram_dir_name):
    """Load precomputed t-SNE results if they exist.
    
    Args:
        perplexity: t-SNE perplexity parameter to look for
        distogram_dir_name: Name of the directory used for distograms
        
    Returns:
        Tuple of (embeddings, filenames, filepaths) if found, None otherwise
    """
    safe_dir_name = distogram_dir_name.replace('/', '_').replace(' ', '_')
    results_file = Path(f"./{input_dir}/tsne_results/tsne_{safe_dir_name}_perplexity_{perplexity}.npz")
    
    if results_file.exists():
        try:
            data = np.load(results_file, allow_pickle=True)
            print(f"\nFound precomputed t-SNE results for perplexity {perplexity} and directory {distogram_dir_name}")
            return (
                data['embeddings'],
                data['filenames'].tolist(),
                data['filepaths'].tolist()
            )
        except Exception as e:
            print(f"Error loading precomputed results: {e}")
            return None
    return None

def determine_special_case(filename):
    """Determine the color and label based on the model source"""
    if 'exp_' in filename.lower():
        return 'green', 'NMR'
    elif 'comp_' in filename.lower():
        return 'cyan', 'Baseline AF2'
    elif 'v3_' in filename and 'dropout' in filename:
        return 'blue', 'AFSample V3'
    elif 'v2_' in filename and 'dropout' in filename:
        return 'pink', 'AFSample V2'
    elif 'dropout' in filename:
        return 'red', 'AFSample'
    elif 'v2_' in filename:
        return 'purple', 'AFSample2 V2'
    elif 'v3_' in filename:
        return 'magenta', 'AFSample2 V3'
    elif 'notemplate' in filename:
        return 'orange', 'AFCluster'
    elif 'multimer' in filename:
        return 'yellow', 'AFSample2 Multimer'
    elif 'cluster' in filename:
        return 'brown', 'Haddock Decoy'
    else:
        return 'gray', 'Other'

def load_ranking_confidences(pdb_id):
    """
    Load ranking confidence values from CSV file.
    
    Args:
        pdb_id: PDB ID to process
        
    Returns:
        DataFrame containing ranking confidence values
    """
    results_file = Path(f"./{pdb_id.upper()}/ranking_confidences.csv")
    if not results_file.exists():
        print(f"Warning: Ranking confidence file not found: {results_file}")
        return None
    
    return pd.read_csv(results_file)

def get_model_basename(pdb_path):
    """
    Extract the model basename from a PDB path by removing beginning and ending tags.
    
    Args:
        pdb_path: Full path to PDB file
        
    Returns:
        Extracted model basename
    """
    # Get the filename without extension
    basename = Path(pdb_path).stem
    
    # Remove common prefix patterns
    prefixes_to_remove = ['min_', 'exp_', 'comp_']
    for prefix in prefixes_to_remove:
        if basename.startswith(prefix):
            basename = basename[len(prefix):]
    
    # Remove common suffix patterns
    suffixes_to_remove = ['_af2', '_aligned', '_haddock_min']
    for suffix in suffixes_to_remove:
        if basename.endswith(suffix):
            basename = basename[:-len(suffix)]
    
    return basename

def find_matching_ranking_confidence(model_basename, ranking_confidences):
    """
    Find the ranking confidence for a model basename using exact matching.
    Uses a dictionary for O(1) lookup.
    
    Args:
        model_basename: The model basename to look up
        ranking_confidences: DataFrame containing model and ranking_confidence columns
        
    Returns:
        The ranking confidence value for the exact model match
        
    Raises:
        ValueError: If no exact match is found for the model basename
    """
    # Create a dictionary mapping model names to ranking confidences for O(1) lookup
    model_to_conf = dict(zip(ranking_confidences['model'].astype(str), ranking_confidences['ranking_confidence']))
    
    # Try exact match
    if model_basename in model_to_conf:
        return model_to_conf[model_basename]
    
    # If no exact match, raise an exception
    available_models = sorted(model_to_conf.keys())
    raise ValueError(f"No exact match found for model '{model_basename}'")

def find_matching_metrics_file(pdb_file, metrics_dir, debug=False):
    """
    Given a pdb_file and a directory, find a .pkl or .json file where the basename (with 'result_' removed)
    matches the pdb_file basename (with 'min_' and '_af2' removed).
    Returns the path to the matching file or None if not found.
    """
    import os
    pdb_basename = os.path.splitext(os.path.basename(str(pdb_file)))[0].split('.')[0]

    pdb_basename = pdb_basename.replace('_af2', '').replace('min_', '').replace('_prot_pept', '')

    if debug:
        print(f"[DEBUG] Searching for matching metrics file for {pdb_basename} in {metrics_dir}")
    
    # List all .pkl and .json files in the directory
    for fname in os.listdir(metrics_dir):
        if fname.endswith('.pkl') or fname.endswith('.json'):
            file_base = os.path.splitext(fname)[0].split('.')[0]
            # Remove 'result_' prefix if present
            if file_base.startswith('result_'):
                file_base = file_base[len('result_'):]
            if file_base == pdb_basename:
                if debug:
                    print(f"[DEBUG] Found matching metrics file: {fname}")
                return os.path.join(metrics_dir, fname)
    if debug:
        print(f"[DEBUG] No matching metrics file found for {pdb_basename}")
    return None

def create_tsne_plot(input_dir, distograms, filenames, filepaths, output_file, perplexity=30, distogram_dir_name="", pdb_dir=None):
    """
    Create interactive t-SNE visualization from distograms using Plotly, colored by model type (determine_special_case).
    """
    # Prompt for residue selection in the main thread before starting PyMOL server
    print("\n" + "="*60)
    print("RESIDUE SELECTION FOR STRUCTURAL ALIGNMENT")
    print("="*60)
    print("Before starting the interactive visualization, please select residue ranges")
    print("that will be used for structural alignment when loading structures in PyMOL.")
    print("These ranges will be applied to all subsequent structures loaded.")
    print()
    
    selected_residue_ranges = prompt_residue_ranges(pdb_dir)
    
    if selected_residue_ranges:
        print(f"\nSelected residue ranges for alignment:")
        for chain, start, end in selected_residue_ranges:
            print(f"  Chain {chain}: {start}-{end}")
    else:
        print("\nNo residue ranges selected. Default alignment will be used.")
    
    print("\n" + "="*60)
    print("STARTING INTERACTIVE VISUALIZATION")
    print("="*60)
    
    # Start PyMOL server and keep references to server and thread
    pymol_port, pymol_server, pymol_thread = start_pymol_server(pdb_dir, selected_residue_ranges)

    # Create custom JavaScript for PyMOL integration (unchanged)
    custom_js = f"""
    <script>
    var pymolPort = {pymol_port};
    var allPlots = [];
    var highlightCounter = 0;
    var peptideColors = [
        'yellow', 'orange', 'red', 'pink', 'purple', 'blue', 'cyan', 'green',
        'lime', 'salmon', 'violet', 'magenta', 'marine', 'teal', 'forest', 'olive'
    ];
    var currentColorIndex = 0;
    function getNextHighlightColor() {{
        var color = peptideColors[currentColorIndex];
        currentColorIndex = (currentColorIndex + 1) % peptideColors.length;
        return color;
    }}
    function logToServer(message, data) {{
        if (message.toLowerCase().includes('error')) {{
            fetch('http://localhost:' + pymolPort + '/log', {{
                method: 'POST',
                headers: {{ 'Content-Type': 'application/json', }},
                body: JSON.stringify({{ message: message, data: data }})
            }}).catch(error => console.error('Error logging to server:', error));
        }}
    }}
    function loadInPyMOL(pdbPath) {{
        var highlightColor = getNextHighlightColor();
        fetch('http://localhost:' + pymolPort, {{
            method: 'POST',
            headers: {{ 'Content-Type': 'application/json', }},
            body: JSON.stringify({{ pdb_path: pdbPath, highlight_color: highlightColor }})
        }})
        .then(response => response.json())
        .then(data => {{
            if (data.status === 'error') {{
                console.error('Error loading in PyMOL:', data.message);
                alert('Error loading structure in PyMOL: ' + data.message);
            }}
        }})
        .catch(error => {{
            console.error('Error:', error);
            alert('Error communicating with PyMOL server');
        }});
        return highlightColor;
    }}
    function highlightPoint(plotDiv, point, highlightColor) {{
        var highlightTrace = {{
            x: [point.x],
            y: [point.y],
            mode: 'markers',
            marker: {{
                size: 15,
                color: 'rgba(0,0,0,0)',
                line: {{ width: 4, color: highlightColor }},
                symbol: 'circle'
            }},
            name: 'highlight_' + highlightCounter++,
            showlegend: false,
            hoverinfo: 'skip'
        }};
        var currentData = plotDiv.data;
        var newData = [...currentData, highlightTrace];
        Plotly.react(plotDiv, newData, plotDiv.layout)
            .catch(error => console.error('Error updating plot:', error));
    }}
    function clearAllHighlights() {{
        allPlots.forEach(function(plotDiv) {{
            var currentData = plotDiv.data;
            var filteredData = currentData.filter(function(trace) {{
                return !trace.name || !trace.name.startsWith('highlight_');
            }});
            Plotly.react(plotDiv, filteredData, plotDiv.layout)
                .catch(error => console.error('Error clearing highlights:', error));
        }});
        highlightCounter = 0;
        console.log('All highlights cleared');
    }}
    
    function clearPyMOLCache() {{
        fetch('http://localhost:' + pymolPort + '/clear_cache', {{
            method: 'POST',
            headers: {{ 'Content-Type': 'application/json', }},
            body: JSON.stringify({{}})
        }})
        .then(response => response.json())
        .then(data => {{
            if (data.status === 'success') {{
                console.log('PyMOL cache cleared successfully');
            }} else {{
                console.error('Error clearing PyMOL cache:', data.message);
            }}
        }})
        .catch(error => {{
            console.error('Error communicating with PyMOL server:', error);
        }});
    }}
    
    function clearAllSelections() {{
        // Clear plot highlights
        clearAllHighlights();
        // Clear PyMOL cache
        clearPyMOLCache();
    }}
    
    function findMatchingPoint(clickedPoint, plotDiv) {{
        function normalizeString(str) {{
            if (!str) return '';
            return str.toString().trim().toLowerCase();
        }}
        function comparePaths(path1, path2) {{
            if (!path1 || !path2) return false;
            var p1 = normalizeString(path1);
            var p2 = normalizeString(path2);
            return p1 === p2 || p1.endsWith(p2) || p2.endsWith(p1);
        }}
        var plotData = plotDiv.data;
        var clickedFilename = clickedPoint.text;
        var clickedPdbPath = clickedPoint.customdata ? clickedPoint.customdata[0] : null;
        var clickedX = clickedPoint.customdata ? clickedPoint.customdata[1] : null;
        var clickedY = clickedPoint.customdata ? clickedPoint.customdata[2] : null;
        for (var i = 0; i < plotData.length; i++) {{
            var trace = plotData[i];
            if (trace.name && trace.name.startsWith('highlight_')) continue;
            if (trace.customdata && Array.isArray(trace.customdata)) {{
                for (var j = 0; j < trace.customdata.length; j++) {{
                    var customData = trace.customdata[j];
                    if (customData && customData[0]) {{
                        var compared_path = customData[0];
                        if (comparePaths(clickedPdbPath, compared_path)) {{
                            var x = customData[1];
                            var y = customData[2];
                            if (typeof x === 'number' && !isNaN(x) && typeof y === 'number' && !isNaN(y)) {{
                                return {{ x: x, y: y, text: clickedFilename }};
                            }}
                        }}
                        if (trace.text && trace.text[j] && normalizeString(clickedFilename) === normalizeString(trace.text[j])) {{
                            var x = customData[1];
                            var y = customData[2];
                            if (typeof x === 'number' && !isNaN(x) && typeof y === 'number' && !isNaN(y)) {{
                                return {{ x: x, y: y, text: clickedFilename }};
                            }}
                        }}
                    }}
                }}
            }}
        }}
        return null;
    }}
    document.addEventListener('DOMContentLoaded', function() {{
        // Create Clear Selections button
        var buttonContainer = document.createElement('div');
        buttonContainer.style.position = 'fixed';
        buttonContainer.style.top = '20px';
        buttonContainer.style.right = '20px';
        buttonContainer.style.zIndex = '1000';
        
        var clearButton = document.createElement('button');
        clearButton.textContent = 'Clear Selections';
        clearButton.style.padding = '10px 15px';
        clearButton.style.backgroundColor = '#ff6b6b';
        clearButton.style.color = 'white';
        clearButton.style.border = 'none';
        clearButton.style.borderRadius = '5px';
        clearButton.style.cursor = 'pointer';
        clearButton.style.fontSize = '14px';
        clearButton.style.fontWeight = 'bold';
        clearButton.style.boxShadow = '0 2px 4px rgba(0,0,0,0.2)';
        
        clearButton.addEventListener('mouseover', function() {{
            this.style.backgroundColor = '#ff5252';
        }});
        
        clearButton.addEventListener('mouseout', function() {{
            this.style.backgroundColor = '#ff6b6b';
        }});
        
        clearButton.addEventListener('click', function() {{
            clearAllSelections();
            this.textContent = 'Cleared!';
            this.style.backgroundColor = '#4caf50';
            setTimeout(() => {{
                this.textContent = 'Clear Selections';
                this.style.backgroundColor = '#ff6b6b';
            }}, 1000);
        }});
        
        buttonContainer.appendChild(clearButton);
        document.body.appendChild(buttonContainer);
        
        // Set up plot event listeners
        var plotContainers = document.querySelectorAll('.plotly-graph-div');
        plotContainers.forEach(function(container) {{
            allPlots.push(container);
            container.on('plotly_click', function(data) {{
                var point = data.points[0];
                var pdbPath = point.customdata[0];
                if (pdbPath) {{
                    var highlightColor = loadInPyMOL(pdbPath);
                    allPlots.forEach(function(targetPlot) {{
                        var matchingPoint = findMatchingPoint(point, targetPlot);
                        if (matchingPoint) {{
                            highlightPoint(targetPlot, matchingPoint, highlightColor);
                        }}
                    }});
                }}
            }});
        }});
    }});
    </script>
    """
    output_file = Path(output_file)
    # t-SNE computation (unchanged)
    precomputed_results = load_tsne_results(input_dir, perplexity, distogram_dir_name)
    if precomputed_results is not None:
        embeddings, precomputed_filenames, precomputed_filepaths = precomputed_results
        if (set(precomputed_filenames) == set(filenames) and set(precomputed_filepaths) == set(filepaths)):
            print("\nPrecomputed results match current files exactly.")
            use_precomputed = True
        else:
            print("\nPrecomputed results exist but don't match current files exactly.")
            print("Differences found in filenames or filepaths.")
            response = input("Would you like to recompute t-SNE? (y/n): ").lower()
            use_precomputed = response != 'y'
    else:
        use_precomputed = False
    if not use_precomputed:
        print("\nFirst 5 distogram shapes:")
        print("------------------------")
        for i, d in enumerate(distograms[:5]):
            print(f"Distogram {i} ({filenames[i]}): {d.shape}")
        print("------------------------\n")
        try:
            flattened_distograms = np.array([d.flatten() for d in distograms])
        except ValueError as e:
            print("\nError during flattening:")
            print("------------------------")
            print(f"Error message: {str(e)}")
            print("\nChecking individual distogram sizes:")
            for i, d in enumerate(distograms):
                print(f"Distogram {i} ({filenames[i]}): {d.shape}, flattened size: {d.flatten().shape}")
            raise
        print(f"\nComputing t-SNE with perplexity {perplexity}...")
        tsne = TSNE(n_components=2, perplexity=perplexity, random_state=42)
        embeddings = tsne.fit_transform(flattened_distograms)
        save_tsne_results(input_dir, embeddings, filenames, filepaths, perplexity, distogram_dir_name)
    else:
        print("Using precomputed t-SNE results.")
    # Only keep minimal DataFrame for model type coloring
    df = pd.DataFrame({
        'x': embeddings[:, 0],
        'y': embeddings[:, 1],
        'filename': filenames,
        'filepath': filepaths
    })
    # Customdata for click/hover: [filepath, x, y]
    df['customdata'] = df.apply(lambda row: [row['filepath'], row['x'], row['y']], axis=1)
    hover_template = 'filename: %{text}<br>filepath: %{customdata[0]}<br>x: %{customdata[1]:.3f}<br>y: %{customdata[2]:.3f}<br><br>Click to load in PyMOL'
    # Group points by category using determine_special_case
    categories = {}
    for idx, row in df.iterrows():
        color, label = determine_special_case(row['filename'])
        if label not in categories:
            categories[label] = {
                'x': [], 'y': [], 'text': [], 'customdata': [], 'color': color
            }
        categories[label]['x'].append(row['x'])
        categories[label]['y'].append(row['y'])
        categories[label]['text'].append(row['filename'])
        categories[label]['customdata'].append(row['customdata'])
    # Create the plot
    fig = go.Figure()
    for label, data in categories.items():
        fig.add_trace(go.Scatter(
            x=data['x'],
            y=data['y'],
            mode='markers',
            marker=dict(color=data['color'], size=10),
            name=label,
            text=data['text'],
            hovertemplate=hover_template,
            customdata=data['customdata']
        ))
    fig.update_layout(
        xaxis_title='t-SNE dimension 1',
        yaxis_title='t-SNE dimension 2',
        hovermode='closest',
        title=f'Interactive t-SNE Visualization of Protein Structures from {distogram_dir_name} (Click points to load in PyMOL)',
        legend=dict(
            yanchor="top",
            y=1,
            xanchor="left",
            x=1.02,
            bgcolor="rgba(255, 255, 255, 0.8)"
        ),
        margin=dict(r=150)
    )
    # Add custom JavaScript
    fig_html = fig.to_html(include_plotlyjs='cdn')
    fig_html = fig_html.replace('</body>', f'{custom_js}</body>')
    with open(output_file, 'w') as f:
        f.write(fig_html)
    print(f"\nPyMOL server is running on port {pymol_port}")
    print(f"Interactive t-SNE plot saved to: {output_file}")
    print("Click on any point in the plot to load the corresponding structure in PyMOL")
    print("The PyMOL server will continue running after this script exits.")
    print("To stop the server, you'll need to close the terminal or press Ctrl+C")
    webbrowser.open(f'file://{output_file.absolute()}')

    # --- NEW: Ranking Confidence Plot ---
    if isinstance(input_dir, Path):
         pdb_id_upper = input_dir.name.upper()
    else:
         pdb_id_upper = str(input_dir).split(os.sep)[-1].upper()
    ranking_conf_file = Path(f"./{pdb_id_upper}/{pdb_id_upper}_ranking_confidences.csv")
    if ranking_conf_file.exists():
         print(f"Found ranking confidence file: {ranking_conf_file}")
         ranking_df = pd.read_csv(ranking_conf_file)
         ranking_confidences = []
         for fname in filenames:
             print(f"DEBUG: fname: {fname}")
             model_basename = fname.replace('_h_atoms', '')
             print(f"DEBUG: model_basename: {model_basename}")
             rc = find_matching_ranking_confidence(model_basename, ranking_df)
             print(f"DEBUG: rc: {rc}")
             ranking_confidences.append(rc)
         df_rank = pd.DataFrame({
             'x': embeddings[:, 0],
             'y': embeddings[:, 1],
             'filename': filenames,
             'ranking_confidence': ranking_confidences
         })
         fig_rank = px.scatter(
             df_rank,
             x='x', y='y',
             color='ranking_confidence',
             color_continuous_scale='Viridis',
             hover_data=['filename', 'ranking_confidence'],
             title=f't-SNE colored by ranking_confidence ({pdb_id_upper})',
             labels={'ranking_confidence': 'Ranking Confidence'}
         )
         fig_rank.update_traces(marker=dict(size=10))
         fig_rank_html = fig_rank.to_html(include_plotlyjs='cdn')
         fig_rank_html = fig_rank_html.replace('</body>', f'{custom_js}</body>')
         output_file_rank = output_file.parent / f"tsne_plot_ranking_confidence.html"
         with open(output_file_rank, 'w') as f:
             f.write(fig_rank_html)
         webbrowser.open(f'file://{output_file_rank.absolute()}')
         print(f"Ranking confidence t-SNE plot saved to: {output_file_rank}")
    else:
         print(f"No ranking_confidence file found at {ranking_conf_file}, skipping ranking confidence plot.")
    # Keep the main thread alive to maintain the server
    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        print("\nShutting down PyMOL server...")
        pymol_server.shutdown()
        pymol_server.server_close()
        print("PyMOL server stopped.")


def update_csp_results_with_metrics(pdb_id, csp_results_df, metrics_dir):
    """
    Update the CSP results DataFrame with metrics from corresponding .pkl or .json files, but only if needed.
    If all metrics columns are present and fully populated, just return the DataFrame as-is.
    """
    metrics_cols = ['ptm', 'iptm', 'ranking_confidence', 'weighted_score']
    # Check if all metrics columns exist and are fully populated (no NaN)
    if all(col in csp_results_df.columns for col in metrics_cols):
        # If all values are present (not null), skip processing
        if csp_results_df[metrics_cols].notnull().all().all():
            print("All metrics columns are present and fully populated in CSP results. Skipping metrics file processing.")
            return csp_results_df
        else:
            print("Some metrics columns have missing values. Will update only missing entries.")
    else:
        # Add missing columns as None
        for col in metrics_cols:
            if col not in csp_results_df.columns:
                csp_results_df[col] = None
        print("Added missing metrics columns to CSP results.")

    print("\nUpdating CSP results with metrics data (only for missing values)...")
    # Create a mapping of PDB basenames to their metrics files, only for missing or invalid ranking_confidence
    metrics_mapping = {}
    from tqdm import tqdm
    for idx, row in tqdm(csp_results_df.iterrows(), desc="Looking for systems with missing confidence values.", unit="row"):
        pdb_file = row['pdb_file']
        rc = row.get('ranking_confidence', None)
        # Only update if rc is missing or not a float in [0, 1]
        needs_update = False
        try:
            rc_float = float(rc)
            if not (0.0 <= rc_float <= 1.0):
                needs_update = True
        except (TypeError, ValueError):
            needs_update = True
        if needs_update:
            metrics_file = find_matching_metrics_file(pdb_file, metrics_dir)
            if metrics_file:
                metrics_mapping[pdb_file] = metrics_file
    print(f"[DEBUG] Found {len(metrics_mapping)} systems with missing confidence values.")
    # Only update rows where ranking_confidence is missing
    from tqdm import tqdm
    for pdb_file, metrics_file in tqdm(metrics_mapping.items(), desc="Updating metrics", unit="file"):
        mask = csp_results_df['pdb_file'] == pdb_file
        needs_update = mask & (
            csp_results_df.loc[mask, 'ranking_confidence'].isnull() |
            csp_results_df.loc[mask, 'ptm'].isnull() |
            csp_results_df.loc[mask, 'iptm'].isnull() |
            csp_results_df.loc[mask, 'weighted_score'].isnull()
        )
        if needs_update.any():
            metrics = process_model_metrics_file(metrics_file)
            if metrics is not None:
                csp_results_df.loc[needs_update, 'ptm'] = metrics['ptm']
                csp_results_df.loc[needs_update, 'iptm'] = metrics['iptm']
                csp_results_df.loc[needs_update, 'ranking_confidence'] = metrics['ranking_confidence']
                csp_results_df.loc[needs_update, 'weighted_score'] = metrics['weighted_score']
                
    # Save the updated DataFrame back to the CSV file
    output_file = Path(f"./{pdb_id.upper()}/csp_results.csv")
    csp_results_df.to_csv(output_file, index=False)
    print(f"Updated CSP results saved to: {output_file}")
    # Print statistics about the update
    total_rows = len(csp_results_df)
    updated_rows = csp_results_df['ranking_confidence'].notna().sum()
    print(f"\nMetrics update statistics:")
    print(f"Total structures: {total_rows}")
    print(f"Structures with metrics: {updated_rows}")
    print(f"Structures without metrics: {total_rows - updated_rows}")
    return csp_results_df

def find_pdb_dirs(pdb_id):
    """
    Find all directories matching the pattern ./PDB_FILES/{pdb_id.upper()}*
    
    Args:
        pdb_id: PDB ID to process
        
    Returns:
        List of Path objects for matching directories
    """
    base_dir = Path("./PDB_FILES")
    if not base_dir.exists():
        raise FileNotFoundError(f"Base directory not found: {base_dir}")
        
    pdb_dirs = list(base_dir.glob(f"{pdb_id.upper()}*"))
    return pdb_dirs

def select_pdb_dir(pdb_id, pdb_dirs):
    """
    Prompt user to select a PDB directory from the list of available directories.
    If only one directory is found, return it without prompting.
    
    Args:
        pdb_id: PDB ID being processed
        pdb_dirs: List of Path objects for available PDB directories
        
    Returns:
        Selected Path object
    """
    if not pdb_dirs:
        raise FileNotFoundError(f"No PDB directories found in ./PDB_FILES/ matching {pdb_id.upper()}")
    
    if len(pdb_dirs) == 1:
        print(f"\nFound single PDB directory: {pdb_dirs[0].name}")
        return pdb_dirs[0]
        
    print("\nAvailable PDB directories:")
    for i, dir_path in enumerate(pdb_dirs, 1):
        print(f"{i}. {dir_path.name}")
    
    while True:
        try:
            choice = input("\nSelect a directory number (or 'q' to quit): ")
            if choice.lower() == 'q':
                sys.exit(0)
                
            idx = int(choice) - 1
            if 0 <= idx < len(pdb_dirs):
                return pdb_dirs[idx]
            else:
                print(f"Please enter a number between 1 and {len(pdb_dirs)}")
        except ValueError:
            print("Please enter a valid number")

def main():
    parser = argparse.ArgumentParser(description='Create interactive t-SNE visualization from distogram files')
    parser.add_argument('pdb_id', help='PDB ID to process (e.g., 5vf0)')
    parser.add_argument('--perplexity', type=float, default=30.0,
                      help='t-SNE perplexity parameter (default: 30.0)')
    parser.add_argument('--force-recompute', action='store_true',
                      help='Force recomputation of t-SNE even if precomputed results exist')
    args = parser.parse_args()
    
    pdb_id = args.pdb_id.upper()
    
    try:
        # Find and select PDB directory
        pdb_dirs = find_pdb_dirs(pdb_id)
        pdb_dir = select_pdb_dir(pdb_id, pdb_dirs)
        print(f"\nSelected PDB directory: {pdb_dir}")
        
        # Find and select distogram directory
        distogram_dirs = find_distogram_dirs(pdb_id)
        if len(distogram_dirs) == 1:
            print(f"\nFound single distogram directory: {distogram_dirs[0].name}")
            distogram_dir = distogram_dirs[0]
        else:
            distogram_dir = select_distogram_dir(pdb_id, distogram_dirs)
        print(f"\nSelected distogram directory: {distogram_dir}")
        
        # Set up output file path
        output_file = Path(f"./{pdb_id}/tsne_plot.html")
        
        # Load and process distograms
        distograms, filenames, filepaths, distogram_dir_name = load_distograms(distogram_dir)
        if not distograms:
            print("No distogram files found!")
            return
        print(f"Loaded {len(distograms)} distogram files")
        
        # Create t-SNE plot
        create_tsne_plot(pdb_id, distograms, filenames, filepaths, output_file, args.perplexity, distogram_dir_name, str(pdb_dir))
        
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return

if __name__ == "__main__":
    main() 