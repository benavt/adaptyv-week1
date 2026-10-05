#!/usr/bin/env python3
"""Build a 21-slide DD1 top-20 interface deck."""
import csv, json, subprocess
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from molecular_colors import COLOR_KEY
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'03_Filtering/PyMol_top20'
IMG=OUT/'renders'
TABLE=ROOT/'03_Filtering/DeltaForge/out/top20_human_egfr_pka_heatmap_table.png'
SVG=ROOT/'03_Filtering/DeltaForge/out/top20_human_egfr_pka_heatmap_table.svg'
GRID=ROOT/'03_Filtering/Boltz2_metrics/out/grid_search_ipsae_x_tm_y_heatmap.svg'
GRID_PNG=OUT/'grid_search_retained_count_heatmap.png'
TOP=ROOT/'03_Filtering/Boltz2_metrics/out/top20_grid_search_candidates.csv'
LONG=ROOT/'03_Filtering/Boltz2_metrics/out/boltz2_dd1_metrics_long.csv'
HISP=ROOT/'03_Filtering/JustHIpKA/top20_DD1_EGFR_D3/top20_hispka_results.csv'
DF=ROOT/'03_Filtering/DeltaForge/out/top20_human_egfr_combined.csv'
DFH=ROOT/'03_Filtering/DeltaForge/out/top20_human_egfr_remote_v2/top20_deltaforge_results.json'
DFM=ROOT/'03_Filtering/DeltaForge/out/top20_mouse_egfr_remote_v2/top20_deltaforge_results.json'
XTM=ROOT/'03_Filtering/PyMol_top20/out/mouse_vs_human_tm_scores.csv'
DECK=OUT/'DD1_top20_EGFR_interface_renderings_v1.pptx'

def textbox(slide,text,x,y,w,h,size=18,bold=False,color=(30,40,50),align=PP_ALIGN.LEFT):
    box=slide.shapes.add_textbox(Inches(x),Inches(y),Inches(w),Inches(h)); tf=box.text_frame; tf.clear()
    p=tf.paragraphs[0]; p.alignment=align; run=p.add_run(); run.text=text; run.font.name='Arial'; run.font.size=Pt(size); run.font.bold=bold; run.font.color.rgb=RGBColor(*color)
    return box

def main():
    subprocess.run(['rsvg-convert','-w','1920','-o',str(TABLE),str(SVG)],check=True)
    subprocess.run(['rsvg-convert','-w','1500','-o',str(GRID_PNG),str(GRID)],check=True)
    rows=list(csv.DictReader(TOP.open()))
    long_rows={(r['design_id'],r['target']):r for r in csv.DictReader(LONG.open())}
    his_rows={(r['design_id'],r['species']):r for r in csv.DictReader(HISP.open())}
    df_rows={r['design_id']:r for r in csv.DictReader(DF.open())}
    df_h={r['design_id']:r for r in json.loads(DFH.read_text())}
    df_m={r['design_id']:r for r in json.loads(DFM.read_text())}
    xtm={r['design_id']:r for r in csv.DictReader(XTM.open())}
    prs=Presentation(); prs.slide_width=Inches(13.333); prs.slide_height=Inches(7.5)
    blank=prs.slide_layouts[6]
    # Slide 1: heuristic grid-search heatmap.
    s=prs.slides.add_slide(blank); s.background.fill.solid(); s.background.fill.fore_color.rgb=RGBColor(248,250,252)
    textbox(s,'DD1 top-20 selection: heuristic grid search',0.45,0.18,12.4,0.35,22,True)
    textbox(s,'Retained design count across the available ipSAE and TM-score-to-original-design thresholds. The red outline marks the selected pair used to derive the top-20 candidates.',0.45,0.58,12.4,0.3,12,False,(70,80,90))
    s.shapes.add_picture(str(GRID_PNG),Inches(2.25),Inches(1.05),width=Inches(8.8))
    textbox(s,'x-axis: ipSAE threshold | y-axis: TM score to original Foundry design',0.45,7.16,12.4,0.18,9,False,(90,100,110),PP_ALIGN.CENTER)
    # Slide 2: existing approved summary table.
    s=prs.slides.add_slide(blank); s.background.fill.solid(); s.background.fill.fore_color.rgb=RGBColor(248,250,252)
    textbox(s,'DD1 top-20: pKa and DeltaForge summary',0.45,0.18,12.4,0.35,22,True)
    s.shapes.add_picture(str(TABLE),Inches(0.25),Inches(0.62),width=Inches(12.85))
    textbox(s,'Human EGFR D3 pKa and remote DeltaForge predictions; sequence labels and full DD1 sequences are shown in the table.',0.45,7.16,12.4,0.18,9,False,(90,100,110))
    # Slides 2-21: one design each, Human left and Mouse right.
    for row in rows:
        rank=int(row['rank']); design=row['design_id']; short=f'DD1-{design.split("_")[-2]}'
        s=prs.slides.add_slide(blank); s.background.fill.solid(); s.background.fill.fore_color.rgb=RGBColor(248,250,252)
        textbox(s,f'{rank:02d}  {short} | EGFR D3 interface',0.45,0.18,12.4,0.36,22,True)
        textbox(s,'Human EGFR D3',0.55,0.62,5.9,0.25,15,True,(35,55,75),PP_ALIGN.CENTER)
        textbox(s,'Mouse EGFR D3',6.88,0.62,5.9,0.25,15,True,(35,55,75),PP_ALIGN.CENTER)
        for species,x in [('Human',0.35),('Mouse',6.68)]:
            img=IMG/f'{rank:02d}_{design}_{species}.png'
            s.shapes.add_picture(str(img),Inches(x),Inches(0.92),width=Inches(6.3),height=Inches(4.05))
        textbox(s,'DD1 sequence: '+row['binder_sequence'],0.55,5.04,12.1,0.18,7.2,False,(70,80,90))
        textbox(s,COLOR_KEY,0.45,5.25,12.4,0.16,7.2,False,(90,100,110),PP_ALIGN.CENTER)
        add_score_table(s, design, long_rows, his_rows, df_h, df_m, xtm, 0.55, 5.48)
    prs.save(DECK); print(DECK)

def add_score_table(slide, design, long_rows, his_rows, df_h, df_m, xtm, x, y):
    table=slide.shapes.add_table(8,3,Inches(x),Inches(y),Inches(12.2),Inches(1.32)).table
    table.columns[0].width=Inches(2.1); table.columns[1].width=Inches(5.05); table.columns[2].width=Inches(5.05)
    vals=[('Metric','Human EGFR D3','Mouse EGFR D3'),
          ('ipSAE',fmt(long_rows[(design,'Human')].get('ipSAE')),fmt(long_rows[(design,'Mouse')].get('ipSAE'))),
          ('TM score',fmt(long_rows[(design,'Human')].get('tm_score')),fmt(long_rows[(design,'Mouse')].get('tm_score'))),
          ('DeltaForge ΔG (kcal/mol)',fmt(df_h[design].get('dg')),fmt(df_m[design].get('dg'))),
          ('DeltaForge Kd (nM)',fmt(df_h[design].get('kd_nm')),fmt(df_m[design].get('kd_nm'))),
          ('Mouse:DD1 vs Human:DD1 TM',fmt(xtm[design].get('mouse_vs_human_tm_score')),fmt(xtm[design].get('mouse_vs_human_tm_score'))),
          ('HIS37 pKa',fmt(his_rows[(design,'Human')].get('HIS37_pKa')),fmt(his_rows[(design,'Mouse')].get('HIS37_pKa'))),
          ('HIS100 pKa',fmt(his_rows[(design,'Human')].get('HIS100_pKa')),fmt(his_rows[(design,'Mouse')].get('HIS100_pKa')))]
    for i,row in enumerate(vals):
        for j,val in enumerate(row):
            cell=table.cell(i,j); cell.text=val; cell.vertical_anchor=MSO_ANCHOR.MIDDLE
            cell.fill.solid(); cell.fill.fore_color.rgb=RGBColor(225,232,239) if i==0 else RGBColor(248,250,252)
            for p in cell.text_frame.paragraphs:
                p.alignment=PP_ALIGN.LEFT if j==0 else PP_ALIGN.CENTER
                for run in p.runs:
                    run.font.name='Arial'; run.font.size=Pt(7.5); run.font.bold=(i==0 or j==0 or i==5); run.font.color.rgb=RGBColor(190,0,0) if i==5 else RGBColor(30,40,50)

def fmt(v):
    try: return f'{float(v):.3f}'
    except (TypeError,ValueError): return '—'
if __name__=='__main__': main()
