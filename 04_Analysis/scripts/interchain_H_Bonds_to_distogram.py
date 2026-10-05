#!/usr/bin/env python3
"""
Script to convert PDB files to distograms based on hydrogen atoms involved in inter-chain hydrogen bonds.
A distogram is an NxN matrix where each element (i,j) represents the distance between hydrogen atoms i and j.
Processes PDB files based on a directory containing PDB files.
"""
import os
import sys
import numpy as np
from Bio import PDB
from Bio.PDB.PDBParser import PDBParser
import argparse
from pathlib import Path
from tqdm import tqdm
from collections import defaultdict
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Sequence import Sequence

def find_interchain_hydrogen_bonds(pdb_file, h_bond_distance_threshold=3.5, h_bond_angle_threshold=30.0):
    """
    Find hydrogen atoms involved in inter-chain hydrogen bonds.
    
    Args:
        pdb_file: Path to the PDB file
        h_bond_distance_threshold: Maximum distance for H-bond (in Angstroms)
        h_bond_angle_threshold: Minimum angle for H-bond (in degrees)
        
    Returns:
        Set of hydrogen atom identifiers in format "chain:residue:atom_id"
    """
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure('protein', pdb_file)
    except Exception as e:
        print(f"Error parsing PDB file {pdb_file}: {e}")
        return set()
    
    # Collect all hydrogen atoms and their parent atoms
    h_atoms = []
    parent_atoms = []
    chains = []
    residues = []
    
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.element == 'H':
                        h_atoms.append(atom)
                        chains.append(chain.id)
                        residues.append(residue.get_id()[1])  # Residue number
                        # Find the parent atom (usually N, O, or S)
                        parent = None
                        min_distance = float('inf')
                        for other_atom in residue:
                            if other_atom.element in ['N', 'O', 'S']:
                                # Check if this atom is covalently bonded to the H atom
                                distance = np.linalg.norm(atom.get_coord() - other_atom.get_coord())
                                if distance < min_distance and distance < 1.5:
                                    parent = other_atom
                                    min_distance = distance
                        parent_atoms.append(parent)
    
    # Collect all electronegative atoms (potential acceptors)
    acceptor_atoms = []
    acceptor_chains = []
    acceptor_residues = []
    
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.element in ['N', 'O', 'S']:
                        acceptor_atoms.append(atom)
                        acceptor_chains.append(chain.id)
                        acceptor_residues.append(residue.get_id()[1])  # Residue number
    
    # Find inter-chain hydrogen bonds
    h_bond_h_atoms = set()
    
    for i, h_atom in enumerate(h_atoms):
        if parent_atoms[i] is None:
            continue
            
        h_chain = chains[i]
        h_residue = residues[i]
        
        for j, acceptor_atom in enumerate(acceptor_atoms):
            acceptor_chain = acceptor_chains[j]
            acceptor_residue = acceptor_residues[j]
            
            # Skip if acceptor is the same as the parent atom
            if acceptor_atom == parent_atoms[i]:
                continue
            
            # Check if this is an inter-chain interaction
            if h_chain != acceptor_chain:
                # Calculate distance between H atom and acceptor
                distance = np.linalg.norm(h_atom.get_coord() - acceptor_atom.get_coord())
                
                if distance <= h_bond_distance_threshold:
                    # Check angle between H-parent and H-acceptor
                    h_coord = h_atom.get_coord()
                    parent_coord = parent_atoms[i].get_coord()
                    acceptor_coord = acceptor_atom.get_coord()
                    
                    # Vector from parent to H
                    parent_h_vector = h_coord - parent_coord
                    # Vector from H to acceptor
                    h_acceptor_vector = acceptor_coord - h_coord
                    
                    # Calculate angle
                    dot_product = np.dot(parent_h_vector, h_acceptor_vector)
                    norm_product = np.linalg.norm(parent_h_vector) * np.linalg.norm(h_acceptor_vector)
                    
                    if norm_product > 0:
                        cos_angle = dot_product / norm_product
                        cos_angle = np.clip(cos_angle, -1.0, 1.0)
                        angle = np.arccos(cos_angle) * 180 / np.pi
                        
                        if angle >= h_bond_angle_threshold:
                            # This is a valid hydrogen bond
                            # Create unique identifier: chain:residue:atom_id
                            h_atom_id = f"{h_chain}:{h_residue}:{h_atom.get_id()}"
                            h_bond_h_atoms.add(h_atom_id)
    
    return h_bond_h_atoms

def get_h_atom_coordinates(pdb_file, h_atom_ids):
    """
    Extract coordinates of specified hydrogen atoms from a PDB file.
    
    Args:
        pdb_file: Path to the PDB file
        h_atom_ids: Set of hydrogen atom identifiers in format "chain:residue:atom_id"
        
    Returns:
        List of (atom_identifier, coordinates) tuples for the specified hydrogen atoms
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    
    coordinates = []
    
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.element == 'H':
                        # Create the same identifier format as in find_interchain_hydrogen_bonds
                        atom_identifier = f"{chain.id}:{residue.get_id()[1]}:{atom.get_id()}"
                        if atom_identifier in h_atom_ids:
                            coordinates.append((atom_identifier, atom.get_coord()))
    
    return coordinates

def collect_all_h_bond_h_atoms(pdb_files):
    """
    Collect all hydrogen atom IDs involved in inter-chain hydrogen bonds across all PDB files.
    
    Args:
        pdb_files: List of PDB file paths
        
    Returns:
        Set of all hydrogen atom identifiers in format "chain:residue:atom_id" involved in inter-chain H-bonds in at least one file
    """
    all_h_bond_h_atoms = set()
    
    print("Finding hydrogen atoms involved in inter-chain hydrogen bonds...")
    for pdb_file in tqdm(pdb_files, desc="Analyzing PDB files"):
        try:
            h_bond_h_atoms = find_interchain_hydrogen_bonds(pdb_file)
            all_h_bond_h_atoms.update(h_bond_h_atoms)
            print(f"Found {len(h_bond_h_atoms)} H-bond hydrogen atoms in {os.path.basename(pdb_file)}")
            if len(h_bond_h_atoms) > 0:
                print(f"  H-atom IDs: {sorted(h_bond_h_atoms)}")
        except Exception as e:
            print(f"Error processing {pdb_file}: {e}")
    
    print(f"\nTotal unique hydrogen atoms involved in inter-chain H-bonds: {len(all_h_bond_h_atoms)}")
    if len(all_h_bond_h_atoms) > 0:
        print(f"All H-atom IDs: {sorted(all_h_bond_h_atoms)}")
    return all_h_bond_h_atoms

def calculate_h_atom_distogram(h_atom_coordinates):
    """
    Calculate the distogram from hydrogen atom coordinates.
    
    Args:
        h_atom_coordinates: List of (atom_id, coordinates) tuples
        
    Returns:
        NxN matrix of distances between hydrogen atoms
    """
    n = len(h_atom_coordinates)
    if n == 0:
        return np.array([])
    
    distogram = np.zeros((n, n))
    
    for i in range(n):
        for j in range(n):
            distogram[i, j] = np.linalg.norm(h_atom_coordinates[i][1] - h_atom_coordinates[j][1])
    
    return distogram

def process_pdb_file_for_h_atoms(pdb_file, all_h_atom_ids, output_dir):
    """
    Process a single PDB file and save its hydrogen atom distogram.
    
    Args:
        pdb_file: Path to the PDB file
        all_h_atom_ids: Set of all hydrogen atom identifiers in format "chain:residue:atom_id" to consider
        output_dir: Directory to save the distogram
    """
    try:
        # Get coordinates of hydrogen atoms present in this file
        h_atom_coords = get_h_atom_coordinates(pdb_file, all_h_atom_ids)
        
        if len(h_atom_coords) == 0:
            print(f"No hydrogen atoms found in {pdb_file}")
            return
        
        # Calculate distogram
        distogram = calculate_h_atom_distogram(h_atom_coords)
        
        if distogram.size == 0:
            print(f"Empty distogram for {pdb_file}")
            return
        
        # Get basename of PDB file
        basename = os.path.splitext(os.path.basename(pdb_file))[0]
        
        # Save distogram
        output_file = os.path.join(output_dir, f"{basename}_h_atoms.npy")
        np.save(output_file, distogram)
        
        # Save atom IDs for reference
        atom_ids_file = os.path.join(output_dir, f"{basename}_h_atom_ids.txt")
        with open(atom_ids_file, 'w') as f:
            f.write(f"Hydrogen atoms found in {basename}:\n")
            f.write("=" * 50 + "\n")
            f.write("Format: chain:residue:atom_id -> coordinates\n")
            f.write("-" * 50 + "\n")
            for atom_identifier, coords in h_atom_coords:
                f.write(f"{atom_identifier} -> {coords}\n")
            f.write(f"\nTotal: {len(h_atom_coords)} hydrogen atoms\n")
        
        print(f"Distogram saved to {output_file} with {len(h_atom_coords)} hydrogen atoms")
        
    except Exception as e:
        print(f"Error processing {pdb_file}: {e}")
        import traceback
        traceback.print_exc()

def get_pdb_paths(directory):
    """
    Get paths to PDB files for a given directory.
    
    Args:
        directory: The directory to search for PDB files
        
    Returns:
        List of paths to PDB files
    """
    pdb_files = []
    base_dir = Path(directory)
    print(f"[DEBUG] Searching for PDB files in {base_dir}")

    return [str(f) for f in base_dir.glob('*.pdb')]

def get_output_dir(output_dir_name):
    """
    Get the output directory for the results.
    
    Args:
        output_dir_name: The output directory name
        
    Returns:
        Path to the output directory
    """
    output_dir = Path(f'./{output_dir_name}_h_bond_distograms')
    return str(output_dir)

def filter_h_atoms_by_availability(pdb_files, all_h_atom_ids):
    """
    Filter hydrogen atom identifiers to remove those that cannot be found in any PDB file.
    
    Args:
        pdb_files: List of PDB file paths
        all_h_atom_ids: Set of all hydrogen atom identifiers in format "chain:residue:atom_id"
        
    Returns:
        Set of hydrogen atom identifiers that can be found in at least one PDB file
    """
    # Track which atoms are found in each file
    atoms_found_per_file = {}
    
    print("Checking availability of hydrogen atoms across all PDB files...")
    for pdb_file in tqdm(pdb_files, desc="Checking H-atom availability"):
        try:
            # Get coordinates for all hydrogen atoms in this file
            h_atom_coords = get_h_atom_coordinates(pdb_file, all_h_atom_ids)
            
            # Store which atoms were found in this file
            atoms_found_per_file[pdb_file] = {atom_identifier for atom_identifier, _ in h_atom_coords}
                
        except Exception as e:
            print(f"Error checking {pdb_file}: {e}")
            atoms_found_per_file[pdb_file] = set()
    
    # Find atoms that are not found in ANY file
    atoms_to_remove = set()
    
    for atom_id in all_h_atom_ids:
        found_in_any_file = False
        for pdb_file in pdb_files:
            if atom_id in atoms_found_per_file[pdb_file]:
                found_in_any_file = True
                break
        
        if not found_in_any_file:
            atoms_to_remove.add(atom_id)
            print(f"Atom {atom_id} not found in any PDB file")
    
    # Remove atoms that cannot be found in any file
    filtered_h_atoms = all_h_atom_ids - atoms_to_remove
    
    print(f"\nFiltering Results:")
    print(f"Original H-atom set size: {len(all_h_atom_ids)}")
    print(f"Atoms to remove: {len(atoms_to_remove)}")
    print(f"Filtered H-atom set size: {len(filtered_h_atoms)}")
    
    # Show which atoms were removed
    if atoms_to_remove:
        print(f"\nRemoved {len(atoms_to_remove)} H-atoms that cannot be found in any PDB file:")
        for atom_id in sorted(atoms_to_remove):
            print(f"  - {atom_id}")
    else:
        print("\nNo atoms were removed - all atoms are present in at least one PDB file")
    
    return filtered_h_atoms

def main():
    parser = argparse.ArgumentParser(description='Convert PDB file(s) to distogram(s) based on inter-chain hydrogen bonds')
    parser.add_argument('directory', help='Directory containing PDB files to process')
    parser.add_argument('output_dir', help='Output directory name')
    parser.add_argument('--h-bond-distance', type=float, default=3.5, 
                       help='Maximum distance for hydrogen bond detection (default: 3.5 Angstroms)')
    parser.add_argument('--h-bond-angle', type=float, default=30.0,
                       help='Minimum angle for hydrogen bond detection (default: 30.0 degrees)')
    args = parser.parse_args()
    
    # Get input paths and output directory
    pdb_files = get_pdb_paths(args.directory)
    output_dir = get_output_dir(args.output_dir)
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    if not pdb_files:
        print(f"No PDB files found in {args.directory}")
        sys.exit(1)
    
    print(f"Found {len(pdb_files)} PDB files to process")
    
    # Step 1: Collect all hydrogen atoms involved in inter-chain hydrogen bonds
    all_h_bond_h_atoms = collect_all_h_bond_h_atoms(pdb_files)
    
    if len(all_h_bond_h_atoms) == 0:
        print("No hydrogen atoms involved in inter-chain hydrogen bonds found!")
        sys.exit(1)
    
    # Step 2: Filter hydrogen atoms by availability across all PDB files
    filtered_h_atoms = filter_h_atoms_by_availability(pdb_files, all_h_bond_h_atoms)
    
    if len(filtered_h_atoms) == 0:
        print("No hydrogen atoms left after filtering!")
        sys.exit(1)
    
    # Step 3: Process each PDB file to create distograms
    print("\nCreating distograms for each PDB file...")
    for pdb_file in tqdm(pdb_files, desc="Processing PDB files"):
        process_pdb_file_for_h_atoms(pdb_file, filtered_h_atoms, output_dir)
    
    # Step 4: Create a summary file with all hydrogen atom IDs
    summary_file = os.path.join(output_dir, "all_h_bond_h_atom_ids.txt")
    with open(summary_file, 'w') as f:
        f.write("Hydrogen atoms involved in inter-chain hydrogen bonds:\n")
        f.write("=" * 60 + "\n")
        f.write("Format: chain:residue:atom_id\n")
        f.write("-" * 60 + "\n")
        for atom_id in sorted(filtered_h_atoms):
            f.write(f"{atom_id}\n")
        f.write(f"\nTotal: {len(filtered_h_atoms)} hydrogen atoms\n")
        f.write("\nNote: Each identifier includes chain ID, residue number, and atom ID\n")
        f.write("Example: A:15:H means chain A, residue 15, atom H\n")
        f.write("\nNote: This list contains only H-atoms that are present in at least one PDB file\n")
    
    print(f"\nProcessing complete! Results saved to {output_dir}")
    print(f"Summary file: {summary_file}")

if __name__ == "__main__":
    main() 