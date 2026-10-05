#!/usr/bin/env python3
"""Create uniquely named SVG scatterplots with heuristic cutoffs highlighted."""
import csv
from pathlib import Path

HERE=Path(__file__).resolve().parent
OUT=HERE/"out"; rows=list(csv.DictReader((OUT/"boltz2_dd1_metrics_paired.csv").open()))

def plot(metric, threshold, label, name, highlight_filter=False):
    pairs=[(float(r[f'{metric}_mouse']),float(r[f'{metric}_human'])) for r in rows if r.get(f'{metric}_mouse') and r.get(f'{metric}_human')]
    w,h,L,T,R,B=760,580,90,55,40,85; lo=min(a for a,b in pairs); hi=max(b for a,b in pairs); pad=(hi-lo)*.08; lo-=pad;hi+=pad
    def x(v): return L+(v-lo)/(hi-lo)*(w-L-R)
    def y(v): return h-B-(v-lo)/(hi-lo)*(h-T-B)
    svg=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}">','<rect width="100%" height="100%" fill="white"/>',f'<line x1="{x(lo)}" y1="{y(lo)}" x2="{x(hi)}" y2="{y(hi)}" stroke="#9ca3af" stroke-dasharray="6,5"/>']
    if threshold is not None: svg += [f'<line x1="{x(threshold)}" y1="{T}" x2="{x(threshold)}" y2="{h-B}" stroke="#dc2626" stroke-width="2" stroke-dasharray="7,4"/>',f'<line x1="{L}" y1="{y(threshold)}" x2="{w-R}" y2="{y(threshold)}" stroke="#dc2626" stroke-width="2" stroke-dasharray="7,4"/>']
    for a,b in pairs:
        passed=((a>.65 and b>.65 and float(next(r for r in rows if float(r[f'{metric}_mouse'])==a and float(r[f'{metric}_human'])==b)['tm_score_mouse'])>.7 and float(next(r for r in rows if float(r[f'{metric}_mouse'])==a and float(r[f'{metric}_human'])==b)['tm_score_human'])>.7) if highlight_filter else (threshold is not None and a>threshold and b>threshold))
        svg.append(f'<circle cx="{x(a):.1f}" cy="{y(b):.1f}" r="4" fill="{"#16a34a" if passed else "#64748b"}" fill-opacity=".72"/>')
    ticks=[]
    for i in range(6):
        v=lo+(hi-lo)*i/5; ticks += [f'<text x="{x(v):.1f}" y="{h-B+25}" text-anchor="middle" font-family="Arial" font-size="12">{v:.2f}</text>',f'<text x="{L-10}" y="{y(v)+4:.1f}" text-anchor="end" font-family="Arial" font-size="12">{v:.2f}</text>']
    svg += [f'<line x1="{L}" y1="{h-B}" x2="{w-R}" y2="{h-B}" stroke="black"/><line x1="{L}" y1="{T}" x2="{L}" y2="{h-B}" stroke="black"/>',*ticks,f'<text x="{w/2}" y="25" text-anchor="middle" font-family="Arial" font-size="18">DD1 {label}: Mouse vs Human</text>',*( [f'<text x="{x(threshold)+6}" y="{T+18}" fill="#b91c1c" font-family="Arial" font-size="14">cutoff {threshold}</text>'] if threshold is not None else []),f'<text x="{w/2}" y="{h-25}" text-anchor="middle" font-family="Arial">{label} — Mouse EGFR</text>',f'<text x="20" y="{h/2}" text-anchor="middle" transform="rotate(-90 20 {h/2})" font-family="Arial">{label} — Human EGFR</text>',('</svg>')]
    (OUT/name).write_text("\n".join(svg)); print(name)

plot('ipSAE',.65,'ipSAE','heuristic_v2_ipSAE.svg')
plot('tm_score',.7,'TM-score to original RFD3','heuristic_v2_tm_score.svg')
plot('ipdae_mean_pae',None,'ipDAE mean PAE','heuristic_v2_ipDAE.svg',highlight_filter=True)
