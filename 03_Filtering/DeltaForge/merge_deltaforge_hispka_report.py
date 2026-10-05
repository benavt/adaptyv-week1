#!/usr/bin/env python3
"""Merge remote DeltaForge and JustHISpKa top-20 outputs and draw an SVG report."""
import csv, json, html
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DF = HERE / "out/top20_human_egfr_remote/top20_deltaforge_results.json"
PK = ROOT / "03_Filtering/JustHIpKA/top20_DD1_EGFR_D3/top20_hispka_results.csv"
OUT = HERE / "out/top20_human_egfr_combined.csv"
SVG = HERE / "out/top20_human_egfr_deltaforge_hispka.svg"

def f(x):
    try: return float(x)
    except (TypeError, ValueError): return None

def main():
    scores = json.loads(DF.read_text())
    pka = list(csv.DictReader(PK.open()))
    pka_by = {(r['rank'], r['design_id'], r['species']): r for r in pka}
    rows=[]
    for d in scores:
        r = pka_by[(str(d['rank']), d['design_id'], 'Human')]
        h37, h100 = f(r['HIS37_pKa']), f(r['HIS100_pKa'])
        flags=[]
        if h37 is not None and not 6.5 <= h37 <= 7.4: flags.append('HIS37')
        if h100 is not None and not 6.5 <= h100 <= 7.4: flags.append('HIS100')
        rows.append({'rank':d['rank'],'design_id':d['design_id'],'binder_sequence':r['binder_sequence'],
                     'deltaforge_dG_kcal_mol':d.get('dg'),'deltaforge_KD_nM':d.get('kd_nm'),
                     'HIS37_pKa':h37,'HIS100_pKa':h100,'pKa_flag':';'.join(flags) or 'none',
                     'deltaforge_scorer':d.get('scorer_version')})
    fields=list(rows[0])
    with OUT.open('w', newline='') as h:
        w=csv.DictWriter(h, fieldnames=fields); w.writeheader(); w.writerows(rows)
    draw(rows)
    print(f'Wrote {OUT}\nWrote {SVG}\nFlagged {sum(r["pKa_flag"] != "none" for r in rows)} designs')

def draw(rows):
    W,H=1500,930; left, top=270,105; rowh=38; chartw=1130
    def x(v, lo, hi): return left + (v-lo)/(hi-lo)*chartw
    def esc(v): return html.escape(str(v))
    parts=[f'''<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">
<title>Top 20 Human EGFR DD1 DeltaForge and JustHISpKa predictions</title>
<rect width="100%" height="100%" fill="#fbfbfd"/><text x="40" y="38" font-family="Arial" font-size="25" font-weight="bold">Top 20 DD1 binders: DeltaForge affinity and JustHISpKa</text>
<text x="40" y="67" font-family="Arial" font-size="15">Human EGFR D3; pKa acceptable range 6.5–7.4 (red markers are outside range)</text>''']
    # four aligned panels
    panels=[('DeltaForge ΔG (kcal/mol)', 'deltaforge_dG_kcal_mol', -11, -5, '#2b6cb0'),('DeltaForge predicted Kd (nM)', 'deltaforge_KD_nM', 0, 400, '#805ad5'),('HIS37 pKa', 'HIS37_pKa', 5.8, 8.0, '#319795'),('HIS100 pKa', 'HIS100_pKa', 5.8, 8.0, '#d69e2e')]
    gap=18; pw=(chartw-3*gap)/4
    for pi,(label,key,lo,hi,color) in enumerate(panels):
        x0=left+pi*(pw+gap); parts.append(f'<text x="{x0+pw/2:.1f}" y="91" text-anchor="middle" font-family="Arial" font-size="14" font-weight="bold">{esc(label)}</text>')
        parts.append(f'<rect x="{x0}" y="{top-18}" width="{pw}" height="{rowh*20+24}" fill="none" stroke="#9aa3ad"/>')
        for t in range(0,5):
            val=lo+(hi-lo)*t/4; xx=x0+pw*t/4
            parts.append(f'<line x1="{xx:.1f}" y1="{top-18}" x2="{xx:.1f}" y2="{top+rowh*20+6}" stroke="#e1e5ea"/>')
            parts.append(f'<text x="{xx:.1f}" y="{top+rowh*20+24}" text-anchor="middle" font-family="Arial" font-size="11">{val:g}</text>')
        if key.startswith('HIS'):
            for val in (6.5,7.4):
                xx=x(val,lo,hi); parts.append(f'<line x1="{xx:.1f}" y1="{top-18}" x2="{xx:.1f}" y2="{top+rowh*20+6}" stroke="#d53f3f" stroke-dasharray="4 3"/>')
        for i,r in enumerate(rows):
            y=top+i*rowh+rowh/2; val=f(r[key])
            if val is None: continue
            bad=key.startswith('HIS') and not 6.5 <= val <= 7.4
            fill='#d53f3f' if bad else color; xx=x(val,lo,hi)
            parts.append(f'<circle cx="{xx:.1f}" cy="{y:.1f}" r="6" fill="{fill}" stroke="#fff" stroke-width="1"><title>Rank {r["rank"]}: {label} {val:.3f}</title></circle>')
    for i,r in enumerate(rows):
        y=top+i*rowh+rowh/2+5; parts.append(f'<text x="{left-12}" y="{y:.1f}" text-anchor="end" font-family="Arial" font-size="12">{int(r["rank"]):02d} {esc(r["design_id"].replace("Human_EGFR_", ""))}</text>')
    parts.append(f'<text x="40" y="{H-24}" font-family="Arial" font-size="13">Red pKa markers: outside 6.5–7.4. Dashed red lines: threshold boundaries. DeltaForge: remote v10.2 unified parallel scorer.</text></svg>')
    SVG.write_text(''.join(parts))

if __name__ == '__main__': main()
