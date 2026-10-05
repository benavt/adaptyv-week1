#!/usr/bin/env python3
"""Run/save DeltaForge sequence requests for the ranked top-20 DD1 set."""
import csv,json,subprocess
from pathlib import Path
HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[1]; OUT=HERE/'out/top20_human_egfr'; OUT.mkdir(parents=True,exist_ok=True)
TARGET=json.loads((ROOT/'01_Staging/EGFR/Human_EGFR_310_501_afs3.json').read_text())['sequences'][0]['protein']['sequence']
SRC=ROOT/'03_Filtering/Boltz2_metrics/out/top20_grid_search_candidates.csv'
def main():
    rows=list(csv.DictReader(SRC.open()))
    for row in rows:
        stem=OUT/f"{int(row['rank']):02d}_{row['design_id']}"
        request={'binder_sequence':row['binder_sequence'],'target_sequence':TARGET,'scorer':'auto','binder_name':row['design_id'],'target_name':'Human_EGFR_D3'}
        (stem.with_suffix('.request.json')).write_text(json.dumps(request,indent=2)+'\n')
        subprocess.run(['python3',str(HERE/'request_deltaforge.py'),'sequence','--output',str(stem),'--binder-sequence',row['binder_sequence'],'--target-sequence',TARGET,'--binder-name',row['design_id'],'--target-name','Human_EGFR_D3'],check=False)
if __name__=='__main__': main()
