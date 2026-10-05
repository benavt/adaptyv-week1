import pymol
from pymol import cmd
import os
from pathlib import Path

# --- Script Configuration ---
PDB_FILE = str(
    Path(__file__).resolve().parent.parent
    / "known_binders"
    / "selected_models"
    / "LANA_ET_selected_models_medoid.pdb"
)
OUTPUT_PSE = 'temp_pocket_view.pse'

# --- Main Visualization Function ---
def generate_view(pdb_path, output_path):
    """
    Loads a PDB file into PyMOL, shows stick representation of the model, visualizes interior surfaces of pockets 
    big enough for at least a water molecule (using small probe radius), with transparency to see interiors,
    hides all H atoms, and saves as .pse file.
    """
    if not os.path.exists(pdb_path):
        print(f"Error: PDB file not found at '{pdb_path}'")
        return

    # Launch PyMOL
    pymol.finish_launching(['pymol', '-cqp'])

    # Set background to black
    cmd.bg_color('black')

    # Load the PDB
    cmd.load(pdb_path, 'target')
    
    # Remove all H atoms
    cmd.remove('target and elem H')

    # Hide everything initially
    cmd.hide('everything', '*')

    # Show stick representation for the entire model
    cmd.show('sticks', 'target')

    # Color chains
    cmd.color('marine', 'chain A')
    cmd.color('salmon', 'chain B')

    # Hide specific residue
    cmd.hide('everything', 'chain A and resi -26')

    # Create separate objects for chains to apply surfaces
    cmd.create('chainA', 'target and chain A')
    cmd.create('chainB', 'target and chain B')

    # Set cavity visualization for chainA (pockets and cavities)
    cmd.set('surface_cavity_mode', 1, 'chainA')  # 1: Cavities & Pockets
    cmd.set('surface_cavity_radius', 1.4, 'chainA')  # Probe radius for water ~1.4 Å
    cmd.set('surface_cavity_cutoff', 2.0, 'chainA')  # Adjust cutoff for small pockets

    # Show surface for chainA pockets, with transparency to see interior
    cmd.show('surface', 'chainA')
    cmd.set('transparency', 0.5, 'chainA')
    cmd.set('surface_color', 'gray', 'chainA')

    # Set cavity visualization for chainB
    cmd.set('surface_cavity_mode', 1, 'chainB')
    cmd.set('surface_cavity_radius', 1.4, 'chainB')
    cmd.set('surface_cavity_cutoff', 2.0, 'chainB')

    # Show surface for chainB pockets, with transparency
    cmd.show('surface', 'chainB')
    cmd.set('transparency', 0.5, 'chainB')
    cmd.set('surface_color', 'gray', 'chainB')

    # Zoom and orient
    cmd.zoom('all', buffer=1.1)
    
    cmd.orient()
    
    # Save the session
    cmd.save(output_path)
    
    # Quit PyMOL
    cmd.quit()
    
    print(f"✅ Script finished. Session saved to '{output_path}'")
    print(f"Opening PyMOL with the session file...")
    
    os.system(f'pymol {output_path}')

# Run the script
if __name__ == "__main__":
    generate_view(PDB_FILE, OUTPUT_PSE)