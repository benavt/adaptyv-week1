#!/usr/bin/env python3
"""Open Human and Mouse PyMOL view-capture sessions with F5 camera saving."""
import csv, json, subprocess
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from molecular_colors import color_commands

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'03_Filtering/PyMol_top20/view_capture'; OUT.mkdir(parents=True,exist_ok=True)
LONG=ROOT/'03_Filtering/Boltz2_metrics/out/boltz2_dd1_metrics_long.csv'
TOP=ROOT/'03_Filtering/Boltz2_metrics/out/top20_grid_search_candidates.csv'
VIEW=ROOT/'01_Staging/EGFR/EGFR_seg_medoid_hotspots_view.json'; REF=ROOT/'01_Staging/EGFR'

def make(species, cif, view_path):
    pml=OUT/f'{species.lower()}_capture.pml'; helper=OUT/f'{species.lower()}_capture.py'
    ref=REF/('Mouse_EGFR_seg_medoid.pdb' if species=='Mouse' else 'Human_EGFR_seg_medoid.pdb')
    pml.write_text(f'''load {cif}, complex
load {ref}, reference
hide everything, reference
show cartoon, complex
show sticks, complex and chain B and resn HIS and (resi 37 or resi 100)
{color_commands(species)}
align complex and chain B and name CA, reference and chain A and name CA
hide everything, reference
show cartoon, complex
set_view ({', '.join(f'{float(v):.9f}' for v in json.loads(VIEW.read_text()))})
bg_color white
set orthoscopic, on
''')
    helper.write_text(f'''from pymol import cmd
VIEW_PATH = r"{view_path}"
def save_view():
    import json
    with open(VIEW_PATH, "w") as h: json.dump(list(cmd.get_view()), h)
    print("Saved fixed view to", VIEW_PATH)
    cmd.quit()
cmd.set_key("F5", save_view)
print("{species} EGFR capture session: set the desired perspective, then press F5.")
''')
    return pml,helper

def main():
    top=list(csv.DictReader(TOP.open()))[0]; rows=list(csv.DictReader(LONG.open())); procs=[]
    for species in ('Human','Mouse'):
        r=next(x for x in rows if x['design_id']==top['design_id'] and x['target']==species)
        cif=next(Path(r['boltz_dir']).glob('*_model_0.cif'))
        pml,helper=make(species,cif,OUT/f'{species.lower()}_view.json')
        procs.append(subprocess.Popen(['pymol',str(pml),'-r',str(helper)]))
        print(f'Opened {species} capture session; press F5 after setting the view.')
    print('Human PID',procs[0].pid,'Mouse PID',procs[1].pid)
if __name__=='__main__': main()
