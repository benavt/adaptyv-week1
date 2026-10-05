#!/usr/bin/env python3
"""Render the top-20 pKa/DeltaForge report as a heatmap table."""
import csv, html
from pathlib import Path

HERE=Path(__file__).resolve().parent
IN=HERE/'out/top20_human_egfr_combined.csv'
OUT=HERE/'out/top20_human_egfr_pka_heatmap_table.svg'

def color(v):
    # clipped red (6.5) -> blue (7.4)
    t=max(0.0,min(1.0,(v-6.5)/0.9))
    r=int(220*(1-t)+35*t); g=int(55*(1-t)+95*t); b=int(55*(1-t)+190*t)
    return f'#{r:02x}{g:02x}{b:02x}'

def text_color(v):
    t=max(0.0,min(1.0,(v-6.5)/0.9))
    return '#ffffff' if t < 0.42 else '#102030'

def main():
    rows=list(csv.DictReader(IN.open()))
    W,H=2050,160+len(rows)*39
    x0=35; widths=[270,780,175,175,190,180]
    xs=[x0]
    for w in widths: xs.append(xs[-1]+w)
    parts=[f'''<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">
<title>Top 20 DD1 pKa heatmap with DeltaForge predictions</title>
<rect width="100%" height="100%" fill="#fbfbfd"/>
<text x="35" y="32" font-family="Arial" font-size="25" font-weight="bold">Top 20 DD1 binders: pKa heatmap and DeltaForge predictions</text>
<text x="35" y="58" font-family="Arial" font-size="14">pKa color scale clipped to 6.5–7.4: red = minimum, blue = maximum. DeltaForge values are reported numerically.</text>
<defs><linearGradient id="scale" x1="0" x2="1"><stop offset="0%" stop-color="#dc3737"/><stop offset="100%" stop-color="#235fbe"/></linearGradient></defs>
<text x="35" y="82" font-family="Arial" font-size="12">Full DD1 sequences are left justified; pKa values outside the color domain are clipped to the endpoint.</text>''']
    header_y=125; row_h=39
    headers=['Design','DD1 amino-acid sequence','HIS37 pKa','HIS100 pKa','DeltaForge ΔG (kcal/mol)','DeltaForge Kd (nM)']
    # Put the heatmap legend directly above the two pKa columns.
    legend_x=xs[2]; legend_w=widths[2]+widths[3]
    parts.append(f'<text x="{legend_x+legend_w/2}" y="79" text-anchor="middle" font-family="Arial" font-size="12">HIS pKa heatmap</text>')
    parts.append(f'<rect x="{legend_x+35}" y="84" width="{legend_w-70}" height="12" fill="url(#scale)"/>')
    parts.append(f'<text x="{legend_x+30}" y="95" text-anchor="end" font-family="Arial" font-size="11">6.5</text>')
    parts.append(f'<text x="{legend_x+legend_w-30}" y="95" font-family="Arial" font-size="11">7.4</text>')
    for i,h in enumerate(headers):
        parts.append(f'<rect x="{xs[i]}" y="{header_y}" width="{widths[i]}" height="{row_h}" fill="#dfe5ec" stroke="#aab3bd"/>')
        parts.append(f'<text x="{xs[i]+widths[i]/2}" y="{header_y+25}" text-anchor="middle" font-family="Arial" font-size="13" font-weight="bold">{html.escape(h)}</text>')
    for j,r in enumerate(rows):
        y=header_y+(j+1)*row_h
        name=f"{int(r['rank']):02d}  DD1-{r['design_id'].split('_')[-2]}"
        vals=[float(r['HIS37_pKa']),float(r['HIS100_pKa'])]
        cells=[name, r['binder_sequence'], f"{vals[0]:.3f}", f"{vals[1]:.3f}", f"{float(r['deltaforge_dG_kcal_mol']):.3f}", f"{float(r['deltaforge_KD_nM']):.3f}"]
        for i,val in enumerate(cells):
            fill=color(vals[i-2]) if i in (2,3) else ('#ffffff' if j%2==0 else '#f4f6f8')
            fg=text_color(vals[i-2]) if i in (2,3) else '#17202a'
            parts.append(f'<rect x="{xs[i]}" y="{y}" width="{widths[i]}" height="{row_h}" fill="{fill}" stroke="#c7cdd4"/>')
            anchor='start' if i in (0,1) else 'middle'; tx=xs[i]+12 if i in (0,1) else xs[i]+widths[i]/2
            font_size = 10 if i == 1 else 13
            parts.append(f'<text x="{tx}" y="{y+25}" text-anchor="{anchor}" font-family="Arial" font-size="{font_size}" fill="{fg}">{html.escape(val)}</text>')
    parts.append('</svg>')
    OUT.write_text(''.join(parts))
    print(OUT)

if __name__=='__main__': main()
