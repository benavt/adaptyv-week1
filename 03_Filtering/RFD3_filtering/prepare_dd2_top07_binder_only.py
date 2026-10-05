#!/usr/bin/env python3
from pathlib import Path
import gemmi
ROOT=Path(__file__).resolve().parents[2]
run=ROOT/'02_Design/foundry/runs/DD2_top07_Human_EGFR_hotspot_display'
st=gemmi.read_structure(str(run/'DD2_top07_Human_EGFR_complex.pdb'))
for model in st:
    for chain in list(model):
        if chain.name != 'B': model.remove_chain(chain.name)
out=run/'DD2_top07_rfd3_input_binder_only.pdb'
st.write_pdb(str(out))
print(out)
