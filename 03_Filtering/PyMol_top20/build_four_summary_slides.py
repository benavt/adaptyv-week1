#!/usr/bin/env python3
import csv,json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from molecular_colors import COLOR_KEY
from pptx import Presentation
from pptx.util import Inches,Pt
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN,MSO_ANCHOR
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'03_Filtering/PyMol_top20/four_summary_slides'; OUT.mkdir(parents=True,exist_ok=True)
REND=ROOT/'03_Filtering/PyMol_top20/renders'
def fmt(v):
    try:return f'{float(v):.3f}'
    except:return '—'

def main():
    top={int(r['rank']):r for r in csv.DictReader((ROOT/'03_Filtering/Boltz2_metrics/out/top20_grid_search_candidates.csv').open())}
    long={(r['design_id'],r['target']):r for r in csv.DictReader((ROOT/'03_Filtering/Boltz2_metrics/out/boltz2_dd1_metrics_long.csv').open())}
    his={(int(r['rank']),r['species']):r for r in csv.DictReader((ROOT/'03_Filtering/JustHIpKA/top20_DD1_EGFR_D3/top20_hispka_results.csv').open())}
    dfh={r['design_id']:r for r in json.loads((ROOT/'03_Filtering/DeltaForge/out/top20_human_egfr_remote_v2/top20_deltaforge_results.json').read_text())}
    dfm={r['design_id']:r for r in json.loads((ROOT/'03_Filtering/DeltaForge/out/top20_mouse_egfr_remote_v2/top20_deltaforge_results.json').read_text())}
    xtm={int(r['rank']):r for r in csv.DictReader((ROOT/'03_Filtering/PyMol_top20/out/mouse_vs_human_tm_scores.csv').open())}
    wanted=[('DD1_062',2),('DD1_02',2),('DD1_07',7),('DD1_11',11)]
    prs=Presentation();prs.slide_width=Inches(13.333);prs.slide_height=Inches(7.5);blank=prs.slide_layouts[6]
    def tx(s,t,x,y,w,h,size=14,bold=False,color=(30,40,50),align=PP_ALIGN.LEFT):
        b=s.shapes.add_textbox(Inches(x),Inches(y),Inches(w),Inches(h));p=b.text_frame.paragraphs[0];p.alignment=align;r=p.add_run();r.text=t;r.font.name='Arial';r.font.size=Pt(size);r.font.bold=bold;r.font.color.rgb=RGBColor(*color)
    for label,rank in wanted:
        row=top[rank];d=row['design_id'];s=prs.slides.add_slide(blank);s.background.fill.solid();s.background.fill.fore_color.rgb=RGBColor(248,250,252)
        tx(s,f'{label} | top-20 rank {rank:02d} | Human and Mouse EGFR D3',.45,.18,12.4,.35,22,True);tx(s,'Human EGFR D3',.55,.62,5.9,.22,14,True,(35,55,75),PP_ALIGN.CENTER);tx(s,'Mouse EGFR D3',6.88,.62,5.9,.22,14,True,(35,55,75),PP_ALIGN.CENTER)
        for sp,x in [('Human',.35),('Mouse',6.68)]:s.shapes.add_picture(str(REND/f'{rank:02d}_{d}_{sp}.png'),Inches(x),Inches(.88),width=Inches(6.3),height=Inches(4.05))
        table=s.shapes.add_table(8,3,Inches(.55),Inches(5.35),Inches(12.2),Inches(1.45)).table;table.columns[0].width=Inches(3.0);table.columns[1].width=Inches(4.6);table.columns[2].width=Inches(4.6)
        vals=[('Metric','Human EGFR D3','Mouse EGFR D3'),('ipSAE',fmt(long[(d,'Human')]['ipSAE']),fmt(long[(d,'Mouse')]['ipSAE'])),('TM score',fmt(long[(d,'Human')]['tm_score']),fmt(long[(d,'Mouse')]['tm_score'])),('DeltaForge ΔG (kcal/mol)',fmt(dfh[d].get('dg')),fmt(dfm[d].get('dg'))),('DeltaForge Kd (nM)',fmt(dfh[d].get('kd_nm')),fmt(dfm[d].get('kd_nm'))),('HIS37 pKa',fmt(his[(rank,'Human')]['HIS37_pKa']),fmt(his[(rank,'Mouse')]['HIS37_pKa'])),('HIS100 pKa',fmt(his[(rank,'Human')]['HIS100_pKa']),fmt(his[(rank,'Mouse')]['HIS100_pKa'])),('Mouse vs Human TM / RMSD',fmt(xtm[rank].get('mouse_vs_human_tm_score'))+' / '+fmt(xtm[rank].get('mouse_vs_human_rmsd'))+' Å',fmt(xtm[rank].get('mouse_vs_human_tm_score'))+' / '+fmt(xtm[rank].get('mouse_vs_human_rmsd'))+' Å')]
        for i,rr in enumerate(vals):
            for j,v in enumerate(rr):
                c=table.cell(i,j);c.text=v;c.vertical_anchor=MSO_ANCHOR.MIDDLE;c.fill.solid();c.fill.fore_color.rgb=RGBColor(225,232,239) if i==0 else RGBColor(248,250,252)
                for p in c.text_frame.paragraphs:
                    p.alignment=PP_ALIGN.LEFT if j==0 else PP_ALIGN.CENTER
                    for r0 in p.runs:r0.font.name='Arial';r0.font.size=Pt(8);r0.font.bold=(i==0 or j==0);r0.font.color.rgb=RGBColor(30,40,50)
        tx(s,'DD1 sequence: '+row['binder_sequence'],.55,4.98,12.1,.18,7.2,False,(70,80,90));tx(s,COLOR_KEY,.45,5.18,12.4,.16,7.0,False,(80,90,100),PP_ALIGN.CENTER)
    prs.save(OUT/'DD1_selected_four_summary_slides.pptx');print(OUT/'DD1_selected_four_summary_slides.pptx')

if __name__=='__main__':main()
