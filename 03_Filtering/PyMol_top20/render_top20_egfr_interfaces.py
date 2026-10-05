#!/usr/bin/env python3
"""Render fixed-view Human/Mouse EGFR-D3 interfaces for all DD1 designs."""
import csv, json, os, subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from molecular_colors import color_commands

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'03_Filtering/PyMol_top20'; IMG=OUT/'renders'; PML=OUT/'pml'
VIEW_DIR=ROOT/'03_Filtering/PyMol_top20/view_capture'
REF_DIR=ROOT/'01_Staging/EGFR'
METRICS=ROOT/'03_Filtering/Boltz2_metrics/out/boltz2_dd1_metrics_long.csv'
SEQ=ROOT/'03_Filtering/RFD3_filtering/out/dd1_asp_near_both_his_sequences.csv'

def write_pml(cif, out_png, hotspot_expr, species):
    view=json.loads((VIEW_DIR/f'{species.lower()}_view.json').read_text())
    view_cmd='set_view ('+', '.join(f'{float(v):.9f}' for v in view)+')'
    ref = REF_DIR/('Mouse_EGFR_seg_medoid.pdb' if species == 'Mouse' else 'Human_EGFR_seg_medoid.pdb')
    return f'''load {cif}, complex
load {ref}, reference
align complex and chain B and name CA, reference and chain A and name CA
hide everything, reference
remove solvent
hide everything
show cartoon, complex
show cartoon, complex and chain B
set cartoon_transparency, 0, complex and chain B
select interface_hotspots, {hotspot_expr}
show sticks, interface_hotspots
color green, interface_hotspots
{color_commands(species)}
set stick_radius, 0.18
set cartoon_fancy_helices, 1
orient complex
zoom complex and (chain A or chain B), 4
{view_cmd}
bg_color white
set orthoscopic, on
set ray_opaque_background, off
set antialias, 2
png {out_png}, width=900, height=650, dpi=150, ray=1
quit
'''

def main():
    OUT.mkdir(exist_ok=True); IMG.mkdir(exist_ok=True); PML.mkdir(exist_ok=True)
    long=list(csv.DictReader(METRICS.open())); designs=sorted({r['design_id'] for r in csv.DictReader(SEQ.open())})
    jobs=[]
    for rank, design in enumerate(designs, 1):
        for species in ('Human','Mouse'):
            matches=[r for r in long if r['design_id']==design and r['target']==species]
            if not matches: raise FileNotFoundError(f'{design} {species}')
            cif=next(Path(matches[0]['boltz_dir']).glob('*_model_0.cif'))
            hotspot='complex and chain B and resi 35+37+67+96+100+122+125'
            stem=f'{design}_{species}'
            png=IMG/f'{stem}.png'; pml=PML/f'{stem}.pml'
            pml.write_text(write_pml(cif,png,hotspot,species))
            jobs.append((pml, png))
    def render(job):
        pml, png = job
        subprocess.run(['pymol','-cq',str(pml)],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    with ThreadPoolExecutor(max_workers=int(os.environ.get('PYMOL_WORKERS','8'))) as pool:
        list(pool.map(render, jobs))
    print(f'Rendered {len(designs)*2} interface views to {IMG}')
if __name__=='__main__': main()
