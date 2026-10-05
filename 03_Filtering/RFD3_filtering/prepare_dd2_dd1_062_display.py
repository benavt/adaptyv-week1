#!/usr/bin/env python3
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "PyMol_top20"))
from molecular_colors import color_commands
import gemmi

ROOT = Path(__file__).resolve().parents[2]
src = next((ROOT / '03_Filtering/Refolding/Putative_DD1_EGFR_D3').glob('DD1_062_Human_*/*model_0.cif'))
run = ROOT / '02_Design/foundry/runs/DD2_DD1_062_Human_EGFR_hotspot_display'
display = run / 'run_display'
display.mkdir(parents=True, exist_ok=True)
st = gemmi.read_structure(str(src))
for model in st:
    for chain in model:
        if chain.name == 'A': chain.name = 'B'
        elif chain.name == 'B': chain.name = 'A'
out = run / 'DD2_DD1_062_Human_EGFR_complex.pdb'
st.write_pdb(str(out))
# Separate RFD3 input: retain only the designed binder chain.
binder_only = gemmi.read_structure(str(out))
for model in binder_only:
    for chain in list(model):
        if chain.name != 'B': model.remove_chain(chain.name)
rfd3_input = run / 'DD2_DD1_062_rfd3_input_binder_only.pdb'
binder_only.write_pdb(str(rfd3_input))
pml = display / 'DD2_DD1_062_Human_EGFR_hotspot_display.pml'
pml.write_text(f'''reinitialize
load "{out}", complex
hide everything, complex
show cartoon, complex
color lightorange, complex and chain A
color cyan, complex and chain B
select dd1_his, complex and chain B and resn HIS
show sticks, dd1_his
color purple, dd1_his
select dd1_acidic, complex and chain B and resn ASP+GLU
show sticks, dd1_acidic
color blue, dd1_acidic
select dd1_interface, (complex and chain B) within 5 of (complex and chain A)
show sticks, dd1_interface
color cyan, dd1_interface
color purple, dd1_his
color blue, dd1_acidic

# Edit this line to choose DD1 hotspot residues for the next design round.
# Example: select dd1_hotspots, complex and chain B and resi 12+18+31
select dd1_hotspots, complex and chain B and resi 3+30+32+60+62
show sticks, dd1_hotspots
color tv_orange, dd1_hotspots

# Original Human EGFR hotspots from DD1_Human_EGFR_hotspots_09_29.
select original_human_egfr_hotspots, complex and chain A and resi 35+37+67+96+100+122+125
show sticks, original_human_egfr_hotspots
color green, original_human_egfr_hotspots

{color_commands("Human", dd1_chain="B", egfr_chain="A")}
set_view (0.231838644, 0.897433579, -0.375316143, -0.923013210, 0.081153467, -0.376114815, -0.307080328, 0.433619887, 0.847159147, -0.000001472, 0.000013672, -154.528198242, -2.875174522, 5.522553921, 0.838042200, 116.729858398, 192.326477051, 20.000000000)
set orthoscopic, on
bg_color white
set ray_opaque_background, off
zoom complex, 4
''')
(display / 'README.md').write_text(f'''# DD2 DD1_062 Human EGFR hotspot display

Source design: `DD1_062` from the Human EGFR species-intersection subset.

PDB: `{out.name}`

Chain convention: chain A = Human EGFR D3; chain B = DD1_062 binder. The source is the Human EGFR Boltz2 model-0 CIF, with chain IDs remapped for the Foundry display convention.

Open with `pymol {pml.name}` from this directory. Edit the `select dd1_hotspots` line to choose DD1 residues for the next design round. The orange selection is the editable hotspot set; `dd1_interface` shows binder residues within 5 Å of EGFR.

Original Human EGFR hotspots from `DD1_Human_EGFR_hotspots_09_29`: residues 35, 37, 67, 96, 100, 122, and 125. They are green sticks except where HIS/NEG/POS residue colors take precedence.

Backbone cartoon convention: Human EGFR (chain A) is lightorange; DD1 (chain B) is cyan. DD2, when present, is magenta. HIS, NEG (ASP/GLU), POS (ARG/LYS), and NAG are purple, blue, red, and yellow sticks respectively, with standard non-carbon element colors. Residue colors take precedence over editable tv_orange and original green hotspot highlights.
''')
print(run)
