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
    Loads a PDB file into PyMOL, shows the model as a cartoon colored by residue charge with prolines and glycines in green,
    shows sticks for sidechain heavy atoms of charged residues at the A/B interface,
    colors those stick atoms by atom type,
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

    # Create separate object for chain A
    cmd.create('chainA', 'target and chain A')
    # Hide specific residue in chainA
    cmd.hide('everything', 'chainA and resi -26')

    # Create separate object for chain B
    cmd.create('chainB', 'target and chain B')

    # Hide the original target
    cmd.hide('everything', 'target')

    # Color all residues white (neutral) by default
    cmd.color('white', 'chainA + chainB')

    # Color negative residues red
    cmd.color('red', '(chainA + chainB) and resn asp+glu')

    # Color positive residues blue
    cmd.color('blue', '(chainA + chainB) and resn lys+arg+his')

    # Color prolines and glycines green
    cmd.color('green', '(chainA + chainB) and resn pro+gly')

    # Show cartoon representation
    cmd.show('cartoon', 'chainA')
    cmd.show('cartoon', 'chainB')

    # Define interface charged selections
    interface_charged_A = 'byres ((chainA within 4.0 of chainB) and (resn asp+glu+lys+arg+his)) and sidechain'
    interface_charged_B = 'byres ((chainB within 4.0 of chainA) and (resn asp+glu+lys+arg+his)) and sidechain'

    # Show sticks for sidechain heavy atoms of charged interface residues
    # cmd.show('sticks', interface_charged_A)
    # cmd.show('sticks', interface_charged_B)

    # Color sticks by atom type
    # cmd.color('atomic', interface_charged_A)
    # cmd.color('atomic', interface_charged_B)

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