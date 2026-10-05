#!/usr/bin/env python3
"""Stage ranked DD1/Human EGFR complexes for Foundry hotspot inspection."""
import csv, json
from pathlib import Path
import gemmi

ROOT=Path(__file__).resolve().parents[2]; TOP=ROOT/'03_Filtering/Boltz2_metrics/out/top20_grid_search_candidates.csv'; MET=ROOT/'03_Filtering/Boltz2_metrics/out/boltz2_dd1_metrics_long.csv'; VIEW=ROOT/'03_Filtering/PyMol_top20/view_capture/human_view.json'; RUNS=ROOT/'02_Design/foundry/runs'

def main():
    tops={int(r['rank']):r for r in csv.DictReader(TOP.open())}; metrics={(r['design_id'],r['target']):r for r in csv.DictReader(MET.open())}; view=json.loads(VIEW.read_text())
    for rank in (2,7,11):
        row=tops[rank]; design=row['design_id']; src=next(Path(metrics[(design,'Human')]['boltz_dir']).glob('*_model_0.cif')); run=RUNS/f'DD2_top{rank:02d}_Human_EGFR_hotspot_display'; disp=run/'run_display'; disp.mkdir(parents=True,exist_ok=True)
        st=gemmi.read_structure(str(src)); model=st[0]; by={c.name:c for c in model};
        # Boltz2: chain A is DD1 and chain B is EGFR. Foundry display convention: A=EGFR, B=DD1.
        by['A'].name='B'; by['B'].name='A'; out=run/f'DD2_top{rank:02d}_Human_EGFR_complex.pdb'; st.write_pdb(str(out))
        pml=disp/f'DD2_top{rank:02d}_Human_EGFR_hotspot_display.pml'; pml.write_text(f'''reinitialize
load "{out}", complex
hide everything, complex
show cartoon, complex
color lightorange, complex and chain A
color purple, complex and chain B
select dd1_interface, (complex and chain B) within 5 of (complex and chain A)
show sticks, dd1_interface
color purple, dd1_interface

# Edit this line to choose DD1 hotspot residues for the next design round.
# Example: select dd1_hotspots, complex and chain B and resi 12+18+31
select dd1_hotspots, none
show sticks, dd1_hotspots
color tv_orange, dd1_hotspots

set_view ({', '.join(f'{float(v):.9f}' for v in view)})
set orthoscopic, on
bg_color white
set ray_opaque_background, off
zoom complex, 4
''')
        (disp/'README.md').write_text(f'''# DD2 top-{rank:02d} hotspot display\n\nDesign: `{design}`\n\nPDB: `{out.name}`\n\nChain convention in this staged display PDB: chain A = Human EGFR D3; chain B = DD2 binder. The source is the Human EGFR Boltz2 model-0 CIF.\n\nOpen with `pymol {pml.name}` from this directory. Edit the `select dd1_hotspots` line in the PML to select DD2 residues on chain B. The orange selection is the editable hotspot set; `dd1_interface` shows DD2 residues within 5 Å of EGFR.\n''')
        print(out)
if __name__=='__main__': main()
