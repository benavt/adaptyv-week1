import pandas as pd
import glob
import os
import sys
import tempfile
from os.path import isdir, exists
import math
import argparse
from os import listdir
from collections import defaultdict
from tqdm import tqdm
import numpy as np
from Bio import PDB
import tkinter as tk
from tkinter import ttk
import readline  # Add at top of file with other imports

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

def get_pdb_sequence(pdb_path):
    """Extract sequence and residue numbers from PDB file."""
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
    if not isdir(directory_path):
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
            sequence_info = get_pdb_sequence(file_path)
            
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

def select_reference_group(sequence_groups):
    """Prompt user to select a reference sequence group."""
    if not sequence_groups:
        print("No sequences available for selection.")
        return None
        
    print("\nAvailable sequence groups:")
    for i, group in enumerate(sequence_groups, 1):
        print(f"\nGroup {i} ({len(group.files)} files):")
        print(group.formatted_numbers)
        print(group.sequence)
        print("-" * len(group.sequence))
    
    while True:
        try:
            selection = int(input("\nEnter the number of the sequence group to use as reference (1-{}): ".format(len(sequence_groups))))
            if 1 <= selection <= len(sequence_groups):
                return sequence_groups[selection - 1]
            print("Invalid selection. Please try again.")
        except ValueError:
            print("Please enter a valid number.")

def align_sequences(ref_seq, target_seq):
    """Perform sequence alignment between reference and target sequences, handling subsequences."""
    # Split sequences at chain breaks
    ref_chains = ref_seq.split(':')
    target_chains = target_seq.split(':')
    
    aligned_indices = []
    current_offset = 0
    
    # Align each chain separately
    for ref_chain, target_chain in zip(ref_chains, target_chains):
        if ref_chain == target_chain:
            # If chains are identical, maintain current indices
            chain_length = len(ref_chain)
            aligned_indices.extend(range(current_offset, current_offset + chain_length))
            current_offset += chain_length
        else:
            # Look for the target chain within the reference chain
            found_pos = ref_chain.find(target_chain)
            if found_pos >= 0:
                # Target is a substring of reference - map indices accordingly
                aligned_indices.extend(range(current_offset + found_pos, 
                                          current_offset + found_pos + len(target_chain)))
            else:
                # For cases where there's no exact substring match
                # Try to find best alignment position by matching characters
                best_match_pos = -1
                best_match_score = -1
                
                # Slide the shorter sequence along the longer one to find best match
                for i in range(len(ref_chain) - len(target_chain) + 1):
                    score = sum(1 for a, b in zip(ref_chain[i:i+len(target_chain)], target_chain) if a == b)
                    if score > best_match_score:
                        best_match_score = score
                        best_match_pos = i
                
                if best_match_score > len(target_chain) * 0.5:  # At least 50% match
                    aligned_indices.extend(range(current_offset + best_match_pos,
                                              current_offset + best_match_pos + len(target_chain)))
                else:
                    # No good match found
                    aligned_indices.extend([-1] * len(target_chain))
        
        current_offset += len(ref_chain)
        
        # Add -1 for the chain break character
        if ref_chain != ref_chains[-1]:
            aligned_indices.append(-1)
            current_offset += 1
    
    return aligned_indices

def adjust_residue_numbers(sequence_groups, reference_group):
    """Adjust residue numbers of all groups to match the reference group."""
    for group in sequence_groups:
        if group == reference_group:
            continue
            
        aligned_indices = align_sequences(reference_group.sequence, group.sequence)
        
        # Create new residue numbers based on alignment
        new_residue_numbers = []
        new_chain_ids = []
        
        for i, aligned_idx in enumerate(aligned_indices):
            if aligned_idx == -1:
                new_residue_numbers.append(-1)
                new_chain_ids.append(' ')
            else:
                new_residue_numbers.append(reference_group.residue_numbers[aligned_idx])
                new_chain_ids.append(reference_group.chain_ids[aligned_idx])
        
        group.residue_numbers = new_residue_numbers
        group.chain_ids = new_chain_ids
        group.formatted_numbers = format_residue_numbers(new_residue_numbers, new_chain_ids)

def print_sequence_alignment(sequence_groups):
    """Print alignment of all unique sequences found."""
    if not sequence_groups:
        print("No sequences found.")
        return
        
    print(f"\nFound {len(sequence_groups)} unique sequences:")
    
    for i, group in enumerate(sequence_groups, 1):
        print(f"\nSequence Group {i} ({len(group.files)} files):")
        print("\nSequence:")
        print(group.formatted_numbers)
        print(group.sequence)
        print("-" * len(group.sequence))

def validate_alignment(sequence_groups, reference_group):
    """Validate that all sequences are properly aligned with the reference sequence."""
    validation_results = []
    
    for i, group in enumerate(sequence_groups):
        if group == reference_group:
            continue
            
        # Check each position in the sequence
        mismatches = []
        ref_seq = reference_group.sequence
        group_seq = group.sequence
        
        # Create sets of (residue, chain_id, residue_number) tuples for each sequence
        ref_residues = set()
        group_residues = set()
        
        # Build reference residues set
        for char, chain_id, res_num in zip(ref_seq, reference_group.chain_ids, reference_group.residue_numbers):
            if char not in [' ', ':']:  # Skip spaces and chain breaks
                ref_residues.add((char, chain_id, res_num))
        
        # Build group residues set
        for char, chain_id, res_num in zip(group_seq, group.chain_ids, group.residue_numbers):
            if char not in [' ', ':']:  # Skip spaces and chain breaks
                group_residues.add((char, chain_id, res_num))
        
        # Find residues that appear in both sequences (matching residue and chain)
        shared_residues = set((c,i) for r, c, i in ref_residues) & set((c,i) for r, c, i in group_residues)
        
        # For each shared residue, check if the residue numbers match
        for chain, index in shared_residues:
            ref_char = next(char for char, c, i in ref_residues if c == chain and i == index)
            group_char = next(char for char, c, i in group_residues if c == chain and i == index)
            
            if ref_char != group_char:
                mismatches.append({
                    'index': index,
                    'chain': chain,
                    'ref_char': ref_char,
                    'group_char': group_char
                })
        
        if mismatches:
            validation_results.append({
                'group': i + 1,
                'mismatches': mismatches
            })
    
    return validation_results

def print_validation_results(validation_results):
    """Print the results of alignment validation."""
    if not validation_results:
        print("\nAlignment Validation: SUCCESS - All sequences are properly aligned with the reference.")
    else:
        print("\nAlignment Validation: FAILED - Found the following issues:")
        for result in validation_results:
            print(f"\nGroup {result['group']} has {len(result['mismatches'])} mismatches:")
            for mismatch in result['mismatches']:
                print(f"  Residue {mismatch['residue']}{mismatch['chain']}: Reference number {mismatch['ref_number']} != Group number {mismatch['group_number']}")

def create_aligned_pdb(input_pdb_path, output_pdb_path, new_residue_numbers, new_chain_ids):
    """Create a new PDB file with updated residue numbers and chain IDs."""
    # Track both old and new residue numbers
    current_old_res = None
    current_new_res = None
    current_new_chain = None
    new_residue_numbers = new_residue_numbers.copy()  # Work with a copy to not modify original
    new_chain_ids = new_chain_ids.copy()
    
    with open(input_pdb_path, 'r') as infile, open(output_pdb_path, 'w') as outfile:
        for line in infile:
            if line.startswith(('ATOM', 'HETATM', 'TER')):
                old_res_num = int(line[22:26])
                
                # If we hit a new residue in the original file, get the next new residue number
                if old_res_num != current_old_res:
                    current_old_res = old_res_num
                    if new_residue_numbers:  # Check if we have more residue numbers
                        current_new_res = new_residue_numbers.pop(0)
                        current_new_chain = new_chain_ids.pop(0)
                        if current_new_res == -1 and current_new_chain == ' ': 
                            # this happens when we have a chain break
                            current_new_res = new_residue_numbers.pop(0)
                            current_new_chain = new_chain_ids.pop(0)
                    else:
                        # If we run out of new residue numbers, something is wrong
                        print(f"Warning: Ran out of new residue numbers at old residue {old_res_num}")
                        current_new_res = -1
                        current_new_chain = ' '
                
                # Only modify lines if we have a valid new residue number
                if current_new_res != -1:
                    # Update residue number and chain ID
                    # Format: columns 1-21, chain(21), residue number(22-25), insertion code(26), rest
                    new_line = (
                        line[:21] +                         # Keep everything before chain ID
                        current_new_chain +                 # New chain ID
                        str(current_new_res).rjust(4) +    # New residue number (right justified in 4 spaces)
                        line[26:]                          # Keep everything after residue number
                    )
                    outfile.write(new_line)
                else:
                    outfile.write(line)
            else:
                outfile.write(line)

def save_aligned_structures(sequence_groups, reference_group, input_directory):
    """Save all structures with aligned residue numbers to a new directory."""
    # Create output directory
    input_dir_name = os.path.basename(input_directory)
    output_directory = os.path.join(os.path.dirname(input_directory), f"{input_dir_name}_sequence_aligned")
    
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
        print(f"\nCreated output directory: {output_directory}")
    
    # Process each sequence group
    total_files = sum(len(group.files) for group in sequence_groups)
    processed_files = 0
    
    print("\nSaving aligned structures:")
    for group in sequence_groups:
        if group == reference_group:
            # For reference group, just copy the files
            for file_path in tqdm(group.files, desc="Processing reference group"):
                output_path = os.path.join(output_directory, os.path.basename(file_path))
                with open(file_path, 'rb') as src, open(output_path, 'wb') as dst:
                    dst.write(src.read())
                processed_files += 1
        else:
            # For other groups, update residue numbers
            for file_path in tqdm(group.files, desc=f"Processing group with {len(group.files)} files"):
                output_path = os.path.join(output_directory, os.path.basename(file_path))
                create_aligned_pdb(
                    file_path,
                    output_path,
                    group.residue_numbers.copy(),  # Pass copies to avoid modifying originals
                    group.chain_ids.copy()
                )
                processed_files += 1
    
    print(f"\nSuccessfully processed {processed_files} files")
    print(f"Aligned structures saved to: {output_directory}")
    return output_directory

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
        
        # Create main frame
        main_frame = ttk.Frame(self)
        main_frame.pack(padx=10, pady=10, fill=tk.BOTH, expand=True)
        
        # Instructions
        ttk.Label(main_frame, text="Click and drag to select ranges. Only residues shared among all structures are shown.").pack()
        
        # Create sequence display
        self.canvas = tk.Canvas(main_frame, height=120, bg='white')
        self.canvas.pack(fill=tk.X, expand=True)
        
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
        
        # Center window
        self.update_idletasks()
        width = self.winfo_width()
        height = self.winfo_height()
        x = (self.winfo_screenwidth() // 2) - (width // 2)
        y = (self.winfo_screenheight() // 2) - (height // 2)
        self.geometry(f'+{x}+{y}')
    
    def draw_sequence(self):
        """Draw the sequence on the canvas."""
        self.canvas.delete('all')
        
        # Calculate dimensions
        char_width = 15
        total_width = len(self.sequence) * char_width
        self.canvas.configure(width=total_width + 20, height=120)
        
        # Draw chain separators and background
        current_x = 10
        current_chain = self.chain_ids[0]
        chain_start_x = current_x
        
        # First pass: Draw backgrounds for chains
        for i, (char, chain) in enumerate(zip(self.sequence, self.chain_ids)):
            if chain != current_chain:
                # Draw previous chain background
                self.canvas.create_rectangle(
                    chain_start_x, 20,
                    current_x, 100,
                    fill='lightgray' if current_chain != ' ' else 'white',
                    outline='gray'
                )
                chain_start_x = current_x
                current_chain = chain
            current_x += char_width
        
        # Draw last chain background
        self.canvas.create_rectangle(
            chain_start_x, 20,
            current_x, 100,
            fill='lightgray' if current_chain != ' ' else 'white',
            outline='gray'
        )
        
        # Second pass: Draw text
        current_x = 10
        for i, (char, chain) in enumerate(zip(self.sequence, self.chain_ids)):
            x_center = current_x + char_width/2
            
            # Only draw residues that are shared
            if i in self.shared_positions:
                # Draw chain ID at top
                self.canvas.create_text(
                    x_center, 30,
                    text=chain,
                    font=('Courier', 10),
                    fill='blue'
                )
                
                # Draw residue number in middle
                self.canvas.create_text(
                    x_center, 50,
                    text=str(self.residue_numbers[i]),
                    font=('Courier', 8)
                )
                
                # Draw amino acid at bottom
                self.canvas.create_text(
                    x_center, 80,
                    text=char,
                    font=('Courier', 14, 'bold')
                )
            elif char == ':':
                # Draw chain break marker
                self.canvas.create_line(
                    current_x, 20,
                    current_x, 100,
                    fill='red',
                    width=2
                )
            else:
                # Draw placeholder for non-shared residues
                self.canvas.create_text(
                    x_center, 80,
                    text='·',
                    font=('Courier', 14),
                    fill='gray'
                )
            
            current_x += char_width
        
        # Draw existing selections
        self.draw_selections()
    
    def draw_selections(self):
        """Draw all selected ranges."""
        char_width = 15
        for chain, start, end in self.ranges:
            # Find positions in sequence
            start_idx = next(i for i, (c, n) in enumerate(zip(self.chain_ids, self.residue_numbers))
                           if c == chain and n == start)
            end_idx = next(i for i, (c, n) in enumerate(zip(self.chain_ids, self.residue_numbers))
                         if c == chain and n == end)
            
            # Draw highlight
            x1 = 10 + start_idx * char_width
            x2 = 10 + (end_idx + 1) * char_width
            self.canvas.create_rectangle(
                x1, 20,
                x2, 100,  # Increased height
                fill='yellow',
                stipple='gray50',
                outline='orange',
                width=2
            )
    
    def start_selection(self, event):
        """Handle mouse button press."""
        char_width = 15
        x = event.x - 10  # Adjust for margin
        idx = x // char_width
        
        if 0 <= idx < len(self.sequence) and idx in self.shared_positions:
            self.current_selection = {
                'start_idx': idx,
                'chain': self.chain_ids[idx],
                'start_res': self.residue_numbers[idx]
            }
    
    def update_selection(self, event):
        """Handle mouse drag."""
        if self.current_selection:
            char_width = 15
            x = event.x - 10  # Adjust for margin
            idx = min(max(0, x // char_width), len(self.sequence) - 1)
            
            if idx in self.shared_positions and self.chain_ids[idx] == self.current_selection['chain']:
                # Draw temporary selection
                self.draw_sequence()  # Redraw to clear previous temp selection
                x1 = 10 + self.current_selection['start_idx'] * char_width
                x2 = 10 + (idx + 1) * char_width
                self.canvas.create_rectangle(x1, 20, x2, 100,
                                          fill='yellow', stipple='gray50')
    
    def end_selection(self, event):
        """Handle mouse button release."""
        if self.current_selection:
            char_width = 15
            x = event.x - 10  # Adjust for margin
            idx = min(max(0, x // char_width), len(self.sequence) - 1)
            
            if idx in self.shared_positions and self.chain_ids[idx] == self.current_selection['chain']:
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

def complete_filename(text, state):
    """Tab completion function for filenames."""
    # text is the current input text
    # state is the index of the completion to return (0 for first match, 1 for second, etc.)
    if not hasattr(complete_filename, "matches"):
        # Get the list of matching files
        if not text:
            complete_filename.matches = complete_filename.available_files[:]
        else:
            complete_filename.matches = [f for f in complete_filename.available_files
                                      if f.lower().startswith(text.lower())]
    
    # Return a match if we have one, or None if we've run out
    try:
        return complete_filename.matches[state]
    except IndexError:
        return None

def prompt_residue_ranges(sequence_groups, reference_group):
    """Prompt user for residue ranges to use in structural alignment using GUI."""
    # Find shared residues
    shared_positions = find_shared_residues(sequence_groups, reference_group)
    
    if not shared_positions:
        print("\nError: No residues are shared among all sequence groups.")
        return []
    
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
                # Ask about reference structure
                while True:
                    response = input("\nWould you like to specify a reference structure for alignment? (y/n): ").lower()
                    if response in ['y', 'yes']:
                        # Get list of available PDB files
                        available_files = []
                        for group in sequence_groups:
                            available_files.extend([os.path.basename(f) for f in group.files])
                        available_files = sorted(set(available_files))  # Remove duplicates
                        
                        # Set up tab completion
                        complete_filename.available_files = available_files
                        readline.set_completer(complete_filename)
                        readline.set_completer_delims(' \t\n;')
                        readline.parse_and_bind('tab: complete')
                        
                        print("\nEnter the name of the reference PDB file.")
                        print("(Type part of the name and press TAB for suggestions)")
                        print("(Enter 'cancel' to use default)")
                        
                        while True:
                            try:
                                ref_file = input("\nReference file: ").strip()
                                if ref_file.lower() == 'cancel':
                                    # Clean up readline settings
                                    readline.set_completer(None)
                                    return merged_ranges
                                if ref_file in available_files:
                                    # Clean up readline settings
                                    readline.set_completer(None)
                                    return merged_ranges, ref_file
                                print("Invalid file name. Please try again.")
                            except EOFError:  # Handle Ctrl+D
                                print("\nCancelled. Using default reference.")
                                readline.set_completer(None)
                                return merged_ranges
                            except KeyboardInterrupt:  # Handle Ctrl+C
                                print("\nCancelled. Using default reference.")
                                readline.set_completer(None)
                                return merged_ranges
                    elif response in ['n', 'no']:
                        return merged_ranges
                    else:
                        print("Please enter 'y' or 'n'")
            elif response in ['n', 'no']:
                response = input("Would you like to try selecting ranges again? (y/n): ").lower()
                if response in ['n', 'no']:
                    return []
                break
            else:
                print("Please enter 'y' or 'n'")

def standardize_atom_name(atom_name):
    """Standardize atom name format."""
    if atom_name[0].isdigit():
        return atom_name[1:] + atom_name[0]
    return atom_name

def read_pdb_atoms(file_path, residue_ranges):
    """Read atoms from PDB file within specified residue ranges.
    
    Args:
        file_path: Path to PDB file
        residue_ranges: List of tuples (chain_id, min_res, max_res)
    """
    atoms = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("ATOM"):
                chain = line[21]
                res_num = int(line[22:26])
                
                # Check if this residue is in any of our ranges
                for chain_id, min_res, max_res in residue_ranges:
                    if chain == chain_id and min_res <= res_num <= max_res:
                        atom_info = standardize_atom_name(line[12:16].strip()) + line[17:26].strip()
                        atoms.append(atom_info)
                        break
    return atoms

def find_exclusive_atoms(list1, list2, file1=None, file2=None):
    """Find atoms exclusive to each list, with their source files.
    
    Args:
        list1: First list of atoms
        list2: Second list of atoms
        file1: Source file path for list1
        file2: Source file path for list2
        
    Returns:
        tuple: (exclusive_to_list1, exclusive_to_list2) where each element is a list of (atom, file) tuples
    """
    exclusive_to_list1 = [(atom, file1) for atom in set(list1) - set(list2)]
    exclusive_to_list2 = [(atom, file2) for atom in set(list2) - set(list1)]
    return exclusive_to_list1, exclusive_to_list2

def remove_atoms_from_pdb(file_path, new_file_path, exclusive_atoms):
    """Remove specified atoms from PDB file."""
    with open(file_path, 'r') as infile, open(new_file_path, 'w') as outfile:
        for line in infile:
            if line.startswith("ATOM"):
                atom_info = standardize_atom_name(line[12:16].strip()) + line[17:26].strip()
                if atom_info not in exclusive_atoms:
                    outfile.write(line)
            else:
                outfile.write(line)

def calculate_transformation(reference_structure, target_structure, residue_ranges):
    """Calculate transformation matrix for structural alignment."""
    ref_atoms = []
    target_atoms = []

    for chain_id, min_res, max_res in residue_ranges:
        ref_atoms.extend([atom.get_coord() for residue in reference_structure[0][chain_id]
                         for atom in residue if min_res <= residue.get_id()[1] <= max_res])
        target_atoms.extend([atom.get_coord() for residue in target_structure[0][chain_id]
                           for atom in residue if min_res <= residue.get_id()[1] <= max_res])

    ref_atoms = np.array(ref_atoms)
    target_atoms = np.array(target_atoms)

    ref_center = np.mean(ref_atoms, axis=0)
    target_center = np.mean(target_atoms, axis=0)

    ref_atoms -= ref_center
    target_atoms -= target_center

    correlation_matrix = np.dot(np.transpose(target_atoms), ref_atoms)
    u, s, vh = np.linalg.svd(correlation_matrix)
    rotation = np.dot(u, vh)
    translation = ref_center - np.dot(target_center, rotation)

    return rotation, translation

def apply_transformation(structure, rotation, translation):
    """Apply transformation to structure."""
    for atom in structure.get_atoms():
        atom.coord = np.dot(atom.coord, rotation) + translation
    return structure

def perform_structural_alignment(sequence_aligned_dir, residue_ranges, reference_file=None):
    """Perform structural alignment on sequence-aligned structures.
    
    Args:
        sequence_aligned_dir: Directory containing sequence-aligned PDB files
        residue_ranges: List of (chain, start, end) tuples for alignment
        reference_file: Optional filename to use as reference structure
    """
    # Create output directory
    structure_aligned_dir = sequence_aligned_dir.replace('_sequence_aligned', '_structure_aligned')
    if not os.path.exists(structure_aligned_dir):
        os.makedirs(structure_aligned_dir)
        print(f"\nCreated output directory: {structure_aligned_dir}")

    # Get all PDB files
    pdb_files = sorted(glob.glob(os.path.join(sequence_aligned_dir, "*.pdb")))
    if not pdb_files:
        print("No PDB files found in sequence aligned directory")
        return

    # Set reference file
    if reference_file:
        ref_path = os.path.join(sequence_aligned_dir, reference_file)
        if os.path.exists(ref_path):
            # Move reference file to front of list
            pdb_files.remove(ref_path)
            pdb_files.insert(0, ref_path)
        else:
            print(f"Warning: Specified reference file {reference_file} not found. Using default.")
            reference_file = os.path.basename(pdb_files[0])
    else:
        reference_file = os.path.basename(pdb_files[0])
    
    print(f"\nUsing {reference_file} as reference structure")

    # Find atoms to remove (atoms not present in reference or any other structure)
    atoms_remove = {}  # Dictionary to store atoms and their source files
    ref_atoms = read_pdb_atoms(pdb_files[0], residue_ranges)
    
    for file_path in tqdm(pdb_files[1:], desc="Finding non-shared atoms"):
        target_atoms = read_pdb_atoms(file_path, residue_ranges)
        exclusive_ref, exclusive_target = find_exclusive_atoms(ref_atoms, target_atoms, 
                                                            pdb_files[0], file_path)
        
        # Update atoms_remove with source file information
        for atom, source_file in exclusive_ref + exclusive_target:
            if atom not in atoms_remove:
                atoms_remove[atom] = set()
            atoms_remove[atom].add(source_file)

    if atoms_remove:
        print("\nThe following atoms will be removed from all structures:")
        for atom in sorted(atoms_remove.keys()):
            files = sorted(atoms_remove[atom])
            num_files = len(files)
            print(f"  {atom} (found in: {num_files} files including {files[0]}")
        
        while True:
            response = input("\nDo you want to continue? (y/n): ").lower()
            if response in ['n', 'no']:
                print("Aborting structural alignment.")
                return
            elif response in ['y', 'yes']:
                break
            print("Please enter 'y' or 'n'")

    # Remove non-shared atoms and perform alignment
    parser = PDB.PDBParser(QUIET=True)
    
    # Process reference structure
    ref_output_path = os.path.join(structure_aligned_dir, os.path.basename(pdb_files[0]))
    remove_atoms_from_pdb(pdb_files[0], ref_output_path, set(atoms_remove.keys()))
    reference_structure = parser.get_structure('reference', ref_output_path)

    # Process and align other structures
    for file_path in tqdm(pdb_files[1:], desc="Aligning structures"):
        output_path = os.path.join(structure_aligned_dir, os.path.basename(file_path))
        
        # Remove non-shared atoms
        remove_atoms_from_pdb(file_path, output_path, set(atoms_remove.keys()))
        
        # Perform structural alignment
        structure = parser.get_structure('target', output_path)
        rotation, translation = calculate_transformation(reference_structure, structure, residue_ranges)
        structure = apply_transformation(structure, rotation, translation)
        
        # Save aligned structure
        io = PDB.PDBIO()
        io.set_structure(structure)
        io.save(output_path)

    print(f"\nStructural alignment complete. Results saved in: {structure_aligned_dir}")
    return structure_aligned_dir

def main():
    parser = argparse.ArgumentParser(description='Find unique sequences in PDB files.')
    parser.add_argument('directory', help='Directory containing PDB files to process')
    
    args = parser.parse_args()
    directory_path = os.path.abspath(args.directory)
    
    print(f"Processing directory: {directory_path}")
    
    sequence_groups = find_unique_sequences(directory_path)
    
    # Select reference group and perform alignment
    reference_group = select_reference_group(sequence_groups)
    if reference_group:
        adjust_residue_numbers(sequence_groups, reference_group)
        print("\nAfter alignment with reference sequence:")
        print_sequence_alignment(sequence_groups)
        
        # Validate the alignment
        validation_results = validate_alignment(sequence_groups, reference_group)
        print_validation_results(validation_results)
        
        if not validation_results:  # Only proceed if validation passed
            # Save sequence-aligned structures
            sequence_aligned_dir = save_aligned_structures(sequence_groups, reference_group, directory_path)
            
            # Get residue ranges for structural alignment
            print("\nNow we'll perform structural alignment.")
            print("Please specify residue ranges to use for superposition.")
            result = prompt_residue_ranges(sequence_groups, reference_group)
            
            if result:
                # Check if reference structure was specified
                if isinstance(result, tuple):
                    residue_ranges, ref_file = result
                else:
                    residue_ranges, ref_file = result, None
                
                # Perform structural alignment
                perform_structural_alignment(sequence_aligned_dir, residue_ranges, ref_file)
            else:
                print("\nNo residue ranges specified. Skipping structural alignment.")
        else:
            print("\nSkipping structure alignment due to validation failures.")

if __name__ == "__main__":
    main()
