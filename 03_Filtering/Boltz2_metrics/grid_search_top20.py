#!/usr/bin/env python3
"""Grid-search symmetric Mouse/Human cutoffs and rank a near-20 DD1 set."""
import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
rows = list(csv.DictReader((OUT / "boltz2_dd1_metrics_paired.csv").open()))
seq_path = HERE.parent / "RFD3_filtering/out/dd1_asp_near_both_his_sequences.csv"
seqs = {r["design_id"]: r["binder_sequence"] for r in csv.DictReader(seq_path.open())}
for r in rows: r["binder_sequence"] = seqs.get(r["design_id"], "")

def f(r, k, species): return float(r[f"{k}_{species}"])
def keep(r, ipsae, tm): return f(r,"ipSAE","mouse") > ipsae and f(r,"ipSAE","human") > ipsae and f(r,"tm_score","mouse") > tm and f(r,"tm_score","human") > tm

grid=[]
for i in range(30, 71):
    ipsae=i/100
    for j in range(55, 81):
        tm=j/100
        n=sum(keep(r,ipsae,tm) for r in rows)
        grid.append({"ipSAE_cutoff":ipsae,"TM_cutoff":tm,"n_retained":n,"distance_from_20":abs(n-20)})
grid.sort(key=lambda x:(x["distance_from_20"], -x["ipSAE_cutoff"], -x["TM_cutoff"]))
chosen=grid[0]
selected=[r for r in rows if keep(r,chosen["ipSAE_cutoff"],chosen["TM_cutoff"])]

def norm(values, v, reverse=False):
    lo,hi=min(values),max(values)
    if hi == lo: return 1.0
    score=(v-lo)/(hi-lo)
    return 1-score if reverse else score

for r in selected:
    ips=[f(r,"ipSAE",s) for s in ("mouse","human")]
    tms=[f(r,"tm_score",s) for s in ("mouse","human")]
    daes=[f(r,"ipdae_mean_pae",s) for s in ("mouse","human")]
    r["ipSAE_mean"] = sum(ips)/2; r["ipSAE_min"] = min(ips)
    r["TM_mean"] = sum(tms)/2; r["TM_min"] = min(tms)
    r["ipDAE_mean"] = sum(daes)/2; r["ipDAE_max"] = max(daes)
for key, reverse in [("ipSAE_mean",False),("TM_mean",False),("ipDAE_mean",True),("ipSAE_min",False),("TM_min",False)]:
    vals=[float(r[key]) for r in selected]
    for r in selected: r[key+"_rankscore"] = norm(vals,float(r[key]),reverse)
for r in selected:
    r["composite_score"]=(r["ipSAE_mean_rankscore"]+r["TM_mean_rankscore"]+r["ipDAE_mean_rankscore"]+r["ipSAE_min_rankscore"]+r["TM_min_rankscore"])/5
selected.sort(key=lambda r:r["composite_score"], reverse=True)
for i,r in enumerate(selected,1): r["rank"]=i

with (OUT/"top20_grid_search_candidates.csv").open("w",newline="") as h:
    fields=["rank","design_id","source_run","composite_score","ipSAE_mean","ipSAE_min","TM_mean","TM_min","ipDAE_mean","ipDAE_max","binder_sequence"]
    w=csv.DictWriter(h,fieldnames=fields,extrasaction="ignore");w.writeheader();w.writerows(selected[:20])
with (OUT/"top20_grid_search.csv").open("w",newline="") as h:
    w=csv.DictWriter(h,fieldnames=list(grid[0]));w.writeheader();w.writerows(grid)
print(f"Chosen cutoffs: ipSAE > {chosen['ipSAE_cutoff']:.2f} and TM-score > {chosen['TM_cutoff']:.2f}; retained {len(selected)} designs; reported top {min(20,len(selected))}.")
print("Top design:",selected[0]["design_id"] if selected else "none")
