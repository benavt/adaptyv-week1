#!/usr/bin/env python3
import csv,re
from pathlib import Path
HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[1]; OUT=HERE/'top20_DD1_EGFR_D3'; rows=list(csv.DictReader((ROOT/'03_Filtering/Boltz2_metrics/out/top20_grid_search_candidates.csv').open()))
data=[]
for row in rows:
    for species in ('Mouse','Human'):
        d=OUT/f"{int(row['rank']):02d}_{row['design_id']}_{species}/predictions_A_HIS37_HIS100"; vals={}
        for n in (37,100):
            m=re.search(r'PKA=([0-9.]+)',(d/f'HIS{n}.txt').read_text()); vals[n]=float(m.group(1)) if m else None
        data.append({'rank':row['rank'],'design_id':row['design_id'],'species':species,'HIS37_pKa':vals[37],'HIS100_pKa':vals[100],'binder_sequence':row['binder_sequence']})
with (OUT/'top20_hispka_results.csv').open('w',newline='') as h:
    w=csv.DictWriter(h,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
W,H,L,T,R,B=1150,740,300,70,60,100; cw=75; ch=24
def cell(v): return f'#{int(239-180*(v-4)/5):02x}{int(246-110*(v-4)/5):02x}{int(255-20*(v-4)/5):02x}'
svg=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}"><rect width="100%" height="100%" fill="white"/><text x="575" y="30" text-anchor="middle" font-family="Arial" font-size="20">Top-20 DD1 EGFR HISpKa predictions</text>']
for i,row in enumerate(rows):
    y=T+i*ch; svg.append(f'<text x="{L-10}" y="{y+17}" text-anchor="end" font-family="Arial" font-size="11">{row["rank"]}. {row["design_id"]}</text>')
    for j,species in enumerate(('Mouse','Human')):
        rr=next(x for x in data if x['rank']==row['rank'] and x['species']==species)
        for k,n in enumerate((37,100)):
            x=L+(j*2+k)*cw; v=rr[f'HIS{n}_pKa']; svg.append(f'<rect x="{x}" y="{y}" width="{cw}" height="{ch}" fill="{cell(v)}" stroke="white"/><text x="{x+cw/2}" y="{y+17}" text-anchor="middle" font-family="Arial" font-size="11">{v:.2f}</text>')
svg += [f'<text x="{L+cw}" y="55" text-anchor="middle" font-family="Arial" font-size="13">Mouse EGFR</text>',f'<text x="{L+3*cw}" y="55" text-anchor="middle" font-family="Arial" font-size="13">Human EGFR</text>',f'<text x="{L+cw/2}" y="{H-55}" text-anchor="middle" font-family="Arial" font-size="11">HIS37</text>',f'<text x="{L+1.5*cw}" y="{H-55}" text-anchor="middle" font-family="Arial" font-size="11">HIS100</text>',f'<text x="{L+2.5*cw}" y="{H-55}" text-anchor="middle" font-family="Arial" font-size="11">HIS37</text>',f'<text x="{L+3.5*cw}" y="{H-55}" text-anchor="middle" font-family="Arial" font-size="11">HIS100</text>','</svg>']
(OUT/'top20_hispka_heatmap.svg').write_text('\n'.join(svg)); print('wrote',len(data),'rows')
