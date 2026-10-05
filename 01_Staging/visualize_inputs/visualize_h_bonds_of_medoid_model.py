import pymol
from pymol import cmd
import os
import tempfile
from pathlib import Path

# --- Script Configuration ---
PDB_FILE = str(
    Path(__file__).resolve().parent.parent
    / "known_binders"
    / "selected_models"
    / "LANA_ET_selected_models_medoid.pdb"
)
OUTPUT_PSE = 'temp_protein_view.pse'

# --- H-bond Detection Function ---
def add_hbond_labels():
    """
    Detects backbone hydrogen bonds between chain A and chain B using explicit hydrogens,
    creates distance objects for each H-bond (H...O), and labels them.
    """
    # Select donors (H attached to backbone N) and acceptors (backbone O)
    cmd.select('donors_A', 'chain A and elem H and neighbor (chain A and name N)')
    
    cmd.select('acceptors_B', 'chain B and name O')
    
    cmd.select('acceptors_A', 'chain A and name O')
    
    cmd.select('donors_B', 'chain B and elem H and neighbor (chain B and name N)')

    # Debug: Print counts and details of donors
    print(f"Number of donor H atoms in chain A: {cmd.count_atoms('donors_A')}")
    cmd.iterate('donors_A', 'print(f"Donor H (chain A): chain {chain} resi {resi} name {name} id {ID} index {index}")')
    
    print(f"Number of donor H atoms in chain B: {cmd.count_atoms('donors_B')}")
    cmd.iterate('donors_B', 'print(f"Donor H (chain B): chain {chain} resi {resi} name {name} id {ID} index {index}")')

    # Find potential pairs (H ... O) within cutoff
    hbond_pairs_AB = cmd.find_pairs('(donors_A)', '(acceptors_B)', cutoff=3.0, mode=0)
    
    hbond_pairs_BA = cmd.find_pairs('(donors_B)', '(acceptors_A)', cutoff=3.0, mode=0)
    
    # Combine pairs
    hbond_pairs = hbond_pairs_AB + hbond_pairs_BA

    print(f"hbond_pairs (indices): {hbond_pairs}")

    # Initialize selection for H-bond atoms
    cmd.select('hbond_atoms', 'none')

    # Filter by angle and create distance objects
    for i, pair in enumerate(hbond_pairs):
        # Atom selections
        obj1 = pair[0][0]
        index1 = pair[0][1]
        obj2 = pair[1][0]
        index2 = pair[1][1]
        
        h_sel = f"{obj1} and index {index1}"
        a_sel = f"{obj2} and index {index2}"
        
        # Select donor heavy atom (N) using distance-based selection
        cmd.select('d_temp', f"{obj1} and elem N within 2.0 of ({h_sel})")
        
        # Check if selection is valid
        if cmd.count_atoms('d_temp') == 0:
            print(f"Warning: No donor N found for H atom index {index1}")
            continue
        
        # Get angle N-H...O
        angle = cmd.get_angle('d_temp', h_sel, a_sel)
        
        # Get residue info for debugging
        stored = {'d_chain': '', 'd_resi': '', 'd_name': '', 'd_id': 0, 'd_index': 0,
                  'h_chain': '', 'h_resi': '', 'h_name': '', 'h_id': 0, 'h_index': 0,
                  'a_chain': '', 'a_resi': '', 'a_name': '', 'a_id': 0, 'a_index': 0}
        cmd.iterate('d_temp', 'stored.d_chain = chain; stored.d_resi = resi; stored.d_name = name; stored.d_id = ID; stored.d_index = index')
        cmd.iterate(h_sel, 'stored.h_chain = chain; stored.h_resi = resi; stored.h_name = name; stored.h_id = ID; stored.h_index = index')
        cmd.iterate(a_sel, 'stored.a_chain = chain; stored.a_resi = resi; stored.a_name = name; stored.a_id = ID; stored.a_index = index')
        
        print(f"Potential H-bond: Donor N (chain {stored['d_chain']} resi {stored['d_resi']} name {stored['d_name']} id {stored['d_id']} index {stored['d_index']}) - "
              f"H (chain {stored['h_chain']} resi {stored['h_resi']} name {stored['h_name']} id {stored['h_id']} index {stored['h_index']}) ... "
              f"Acceptor O (chain {stored['a_chain']} resi {stored['a_resi']} name {stored['a_name']} id {stored['a_id']} index {stored['a_index']}), angle={angle:.2f}")

        # Clean up temp selection
        cmd.delete('d_temp')
        
        # If angle > 120 degrees, create distance
        if angle > 120:
            # Add to hbond_atoms selection
            cmd.select('hbond_atoms', f"hbond_atoms or ({h_sel}) or ({a_sel})")
            cmd.select('d_temp', f"{obj1} and elem N within 2.0 of ({h_sel})")  # Re-select N for adding
            cmd.select('hbond_atoms', 'hbond_atoms or d_temp')
            cmd.delete('d_temp')
            
            dist_name = f"hbond_{i+1}"
            
            cmd.distance(dist_name, h_sel, a_sel)
            
            # Customize
            cmd.set('dash_color', 'yellow', dist_name)
            
            cmd.set('dash_radius', 0.05, dist_name)
            
            cmd.set('label_color', 'yellow', dist_name)
            
            cmd.set('label_size', 12)
            
            distance_value = cmd.get_distance(h_sel, a_sel)
            cmd.label(dist_name, f"{distance_value:.2f}")

    # Apply coloring and ensure visibility for H-bond atoms
    cmd.show('sticks', 'hbond_atoms')
    cmd.color('atomic', 'hbond_atoms')

    # Clean up selections
    cmd.delete('donors_A')
    cmd.delete('acceptors_B')
    cmd.delete('donors_B')
    cmd.delete('acceptors_A')
    # Note: hbond_atoms is kept for the session, but can be deleted if not needed
    # cmd.delete('hbond_atoms')

# --- Main Visualization Function ---
def generate_view(pdb_path, output_path):
    """
    Loads a PDB file into PyMOL, applies visualization and selection commands,
    adds H-bond labels between chain A and B, and saves as .pse file.
    """
    # Check if the file exists
    if not os.path.exists(pdb_path):
        print(f"Error: PDB file not found at '{pdb_path}'")
        return

    # Launch PyMOL
    pymol.finish_launching(['pymol', '-cqp'])

    # Set background to black
    cmd.bg_color('black')

    # Load the PDB
    cmd.load(pdb_path, 'target')
    
    # Force re-add hydrogens to standardize across chains
    cmd.remove('target and elem H')
    cmd.h_add('target')

    # Hide specific residue
    cmd.hide('everything', 'chain A and resi -26')

    # Show chain B as sticks
    cmd.show('sticks', 'chain B')
    
    cmd.hide('cartoon', 'chain B')

    # Select near chain B
    cmd.select('near_chainB_exclusive', '(all and not chain B) within 5 of (chain B)')

    # Show sticks
    cmd.show('sticks', 'near_chainB_exclusive')

    # Hide cartoons
    cmd.hide('cartoon', '*')

    # Show cartoon for A
    cmd.show('cartoon', 'chain A')
    
    cmd.hide('everything', 'chain A and resi -26')
    # For 68-75
    cmd.hide('cartoon', 'chain A and resi 68-75')
    
    cmd.show('sticks', 'chain A and resi 68-75')

    # Colors
    cmd.color('marine', 'chain A')
    
    cmd.color('salmon', 'chain B')
    
    cmd.color('orange', 'near_chainB_exclusive')

    # Hide sidechain atoms in sticks representation
    cmd.hide('sticks', 'sidechain')
    
    # Add H-bond labels
    add_hbond_labels()
    
    cmd.zoom('all', buffer=1.1)
    
    cmd.orient()
    
    # Save
    cmd.save(output_path)
    
    # Quit
    cmd.quit()
    
    print(f"✅ Script finished. Session saved to '{output_path}'")
    print(f"Opening PyMOL with the session file...")
    
    os.system(f'pymol {output_path}')

# Run
if __name__ == "__main__":
    generate_view(PDB_FILE, OUTPUT_PSE)