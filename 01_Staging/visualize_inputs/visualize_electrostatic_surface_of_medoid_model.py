import pymol
from pymol import cmd, util
import os
from pathlib import Path

# --- Script Configuration ---
PDB_FILE = str(
    Path(__file__).resolve().parent.parent
    / "known_binders"
    / "selected_models"
    / "LANA_ET_selected_models_medoid.pdb"
)
OUTPUT_PSE = 'temp_electrostatic_view.pse'

# --- Main Visualization Function ---
def generate_view(pdb_path, output_path):
    """
    Loads a PDB file into PyMOL, creates separate objects for chain A and chain B with their electrostatic surfaces,
    allowing toggling between them in the PyMOL GUI, hides all H atoms, and saves as .pse file.
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

    # Create separate object for chain A
    cmd.create('chainA', 'target and chain A')
    # Hide specific residue in chainA
    cmd.hide('everything', 'chainA and resi -26')

    # Create separate object for chain B
    cmd.create('chainB', 'target and chain B')

    # Hide the original target
    cmd.hide('everything', 'target')

    # Generate and apply vacuum electrostatics to chainA
    util.protein_vacuum_esp('chainA')
    cmd.show('surface', 'chainA')

    # Generate and apply vacuum electrostatics to chainB
    util.protein_vacuum_esp('chainB')
    cmd.show('surface', 'chainB')

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