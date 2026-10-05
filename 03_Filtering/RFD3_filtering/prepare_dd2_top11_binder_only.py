#!/usr/bin/env python3
from pathlib import Path
import gemmi
ROOT=Path(__file__).resolve().parents[2]
run=ROOT/'02_Design/foundry/runs/DD2_top11_Human_EGFR_hotspot_display'
src=run/'DD2_top11_Human_EGFR_complex.pdb'
st=gemmi.read_structure(str(src))
for model in st:
    for chain in list(model):
        if chain.name != 'B': model.remove_chain(chain.name)
out=run/'DD2_top11_rfd3_input_binder_only.pdb'
st.write_pdb(str(out))
print(out)
