#!/usr/bin/env python3
"""Plot retained designs with ipSAE on x and TM score on y."""
import csv
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]; OUT=ROOT/'03_Filtering/Boltz2_metrics/out'; rows=list(csv.DictReader((OUT/'top20_grid_search.csv').open()))
xs=sorted({float(r['ipSAE_cutoff']) for r in rows}); ys=sorted({float(r['TM_cutoff']) for r in rows}); vals={(float(r['ipSAE_cutoff']),float(r['TM_cutoff'])):int(r['n_retained']) for r in rows}; sel=min(rows,key=lambda r:float(r['distance_from_20']))
W,H,L,T,R,B=920,700,110,70,80,105; cw=(W-L-R)/len(xs); ch=(H-T-B)/len(ys); vmax=max(vals.values())
def color(v):
 q=min(1,v/max(vmax,1)); return f'#{int(239-150*q):02x}{int(246-110*q):02x}{int(255-20*q):02x}'
svg=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}"><rect width="100%" height="100%" fill="white"/>',f'<text x="{W/2}" y="30" text-anchor="middle" font-family="Arial" font-size="20">DD1 grid search: designs retained</text>']
for yi,yv in enumerate(ys):
 for xi,xv in enumerate(xs):
  v=vals[(xv,yv)]; x=L+xi*cw; y=H-B-(yi+1)*ch; svg.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{cw+0.5:.1f}" height="{ch+0.5:.1f}" fill="{color(v)}" stroke="white" stroke-width="0.5"/>')
  if abs(v-20)<=1: svg.append(f'<text x="{x+cw/2:.1f}" y="{y+ch/2+4:.1f}" text-anchor="middle" font-family="Arial" font-size="9" fill="#102A43">{v}</text>')
for xi,xv in enumerate(xs):
 if xi%5==0 or xv==float(sel['ipSAE_cutoff']): svg.append(f'<text x="{L+(xi+.5)*cw:.1f}" y="{H-B+22}" text-anchor="middle" font-family="Arial" font-size="11">{xv:.2f}</text>')
for yi,yv in enumerate(ys):
 if yi%5==0 or yv==float(sel['TM_cutoff']): svg.append(f'<text x="{L-12}" y="{H-B-(yi+.5)*ch+4:.1f}" text-anchor="end" font-family="Arial" font-size="11">{yv:.2f}</text>')
sx=L+(xs.index(float(sel['ipSAE_cutoff']))+.5)*cw; sy=H-B-(ys.index(float(sel['TM_cutoff']))+.5)*ch
svg += [f'<rect x="{sx-cw/2+1:.1f}" y="{sy-ch/2+1:.1f}" width="{cw-2:.1f}" height="{ch-2:.1f}" fill="none" stroke="#DC2626" stroke-width="3"/>',f'<text x="{W/2}" y="{H-35}" text-anchor="middle" font-family="Arial">ipSAE cutoff</text>',f'<text x="24" y="{H/2}" text-anchor="middle" transform="rotate(-90 24 {H/2})" font-family="Arial">TM-score cutoff to original design</text>',f'<text x="{L}" y="{H-B+48}" font-family="Arial" font-size="13" fill="#DC2626">Red outline: selected pair ipSAE &gt; {float(sel["ipSAE_cutoff"]):.2f}, TM-score &gt; {float(sel["TM_cutoff"]):.2f}, n = {sel["n_retained"]}</text>','</svg>']
(OUT/'grid_search_ipsae_x_tm_y_heatmap.svg').write_text('\n'.join(svg))
print(OUT/'grid_search_ipsae_x_tm_y_heatmap.svg')
