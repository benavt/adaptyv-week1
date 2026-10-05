#!/usr/bin/env python3
"""Collect Boltz2 confidence/interface metrics and Mouse-vs-Human plots."""
from __future__ import annotations
import argparse, csv, gzip, json, re, subprocess, sys
from pathlib import Path
import numpy as np
from Bio.PDB import MMCIFParser

ROOT = Path(__file__).resolve().parents[2]
REFOLDS = ROOT / "03_Filtering/Refolding/Putative_DD1_EGFR_D3"
SEQ_REPORT = ROOT / "03_Filtering/RFD3_filtering/out/dd1_asp_near_both_his_sequences.csv"
ORIGINAL_ROOT = ROOT / "02_Design/foundry/runs"


def conf_file(d): return next(d.glob("confidence*.json"), None)
def pae_file(d): return next(d.glob("pae*.npz"), None)
def plddt_file(d): return next(d.glob("plddt*.npz"), None)
def score_file(d): return next(d.glob("*_15_15.txt"), None)
def model_file(d): return next(d.glob("*_model_0.cif"), None)
AA = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}

def cif_sequences(model):
    s = MMCIFParser(QUIET=True).get_structure("index", model)
    return tuple("".join(AA.get(r.resname, "X") for r in s[0][c] if r.id[0] == " ") for c in ("A", "B"))

def index_models():
    index = {}
    for d in REFOLDS.iterdir():
        if not d.is_dir(): continue
        m = model_file(d)
        if m:
            try: index[cif_sequences(m)] = d
            except Exception: pass
    return index


def parse_ipsae(path):
    if not path: return {}
    lines = [x.split() for x in path.read_text().splitlines() if x.strip() and not x.startswith("#")]
    if len(lines) < 2: return {}
    header, values = lines[0], lines[1:]
    out = {}
    for row in values:
        if len(row) < len(header): continue
        rec = dict(zip(header, row))
        if rec.get("Chn1") == "A" and rec.get("Chn2") == "B" and rec.get("Type") == "max":
            for key in ("ipSAE", "ipSAE_d0chn", "ipSAE_d0dom", "ipTM_af", "ipTM_d0chn", "pDockQ", "pDockQ2", "LIS"):
                try: out[key.lower()] = float(rec[key])
                except (KeyError, ValueError): pass
            return out
    return out


def chain_ca(structure, chain_id):
    chain = structure[0][chain_id]
    return {r.id[1]: r["CA"].coord for r in chain if r.id[0] == " " and "CA" in r}


def rmsd_to_original(predicted, original):
    p = MMCIFParser(QUIET=True).get_structure("pred", predicted)
    oh = gzip.open(original, "rt") if str(original).endswith(".gz") else original
    o = MMCIFParser(QUIET=True).get_structure("orig", oh)
    pxyz, oxyz = [], []
    for chain in ("A", "B"):
        if chain not in p[0] or chain not in o[0]: continue
        pc, oc = chain_ca(p, chain), chain_ca(o, chain)
        ids = sorted(set(pc) & set(oc))
        pxyz.extend(pc[i] for i in ids); oxyz.extend(oc[i] for i in ids)
    if len(pxyz) < 3: return (None, 0)
    x, y = np.asarray(pxyz, float), np.asarray(oxyz, float)
    x -= x.mean(0); y -= y.mean(0)
    u, _, vt = np.linalg.svd(x.T @ y)
    d = np.linalg.det(u @ vt)
    rot = u @ np.diag([1, 1, d]) @ vt
    aligned = x @ rot
    return (float(np.sqrt(np.mean(np.sum((aligned - y) ** 2, axis=1)))), len(ids))


def interface_pae(d, model):
    pf = pae_file(d)
    if not pf or not model: return None
    try:
        pae = np.load(pf)["pae"]
        s = MMCIFParser(QUIET=True).get_structure("x", model)
        lengths = [len([r for r in c if r.id[0] == " "]) for c in s[0]]
        n = lengths[0]
        return float(np.mean(np.concatenate([pae[:n, n:].ravel(), pae[n:, :n].ravel()])))
    except Exception:
        return None

def tm_to_original(predicted, original):
    """Run Foundry's USalign when a runnable build is available."""
    binary = Path(__import__("os").environ.get("USALIGN", "/tmp/USalign_mac"))
    if not binary.is_file(): return None
    try:
        import tempfile, shutil
        with tempfile.TemporaryDirectory() as td:
            op = Path(td) / "original.cif"
            if str(original).endswith(".gz"):
                with gzip.open(original, "rb") as src, op.open("wb") as dst: shutil.copyfileobj(src, dst)
            else: shutil.copyfile(original, op)
            result = subprocess.run([str(binary), str(op), str(predicted), "-mol", "prot", "-mm", "1", "-ter", "1", "-het", "1"], capture_output=True, text=True, check=False)
        rmsds = re.findall(r"Aligned length=\s*\d+\s*, RMSD=\s*([0-9.]+)", result.stdout)
        tms = re.findall(r"TM-score=\s*([0-9.]+)", result.stdout)
        return (float(tms[0]), float(rmsds[0])) if tms and rmsds else None
    except Exception: return None


def collect(row, target, idx, model_index):
    # Refolding names are DD1_<index>_<source>_<target>_EGFR_<design_id>.
    source = "Mouse" if row["run"] == "mouse" else "Human"
    target_seq = json.loads((ROOT / f"01_Staging/EGFR/{target}_EGFR_310_501_afs3.json").read_text())["sequences"][0]["protein"]["sequence"]
    d = model_index.get((row["binder_sequence"], target_seq))
    if d is None:
        tail = re.sub(r"^[A-Za-z]+_EGFR_", "", row["design_id"])
        d = REFOLDS / f"DD1_{idx:03d}_{source}_{target}_EGFR_{tail}"
    cf, model = conf_file(d), model_file(d)
    if not cf or not model: return {"design_id": row["design_id"], "target": target, "status": "missing_confidence_or_model"}
    c = json.loads(cf.read_text())
    plddt = None
    if (pf := plddt_file(d)):
        plddt = float(np.mean(np.load(pf)["plddt"]) * 100)
    ips = parse_ipsae(score_file(d))
    orig_dir = ORIGINAL_ROOT / f"DD1_{source}_EGFR_hotspots_09_29/rfd3/outputs"
    original = orig_dir / f"{row['design_id']}.cif.gz"
    rmsd, nali = rmsd_to_original(model, original) if original.exists() else (None, 0)
    tm = tm_to_original(model, original) if original.exists() else None
    return {"design_id": row["design_id"], "source_run": row["run"], "target": target, "status": "ok", "boltz_dir": str(d),
            "confidence_score": c.get("confidence_score"), "ptm": c.get("ptm"), "iptm": c.get("iptm"),
            "complex_plddt": c.get("complex_plddt", plddt / 100 if plddt else None) * 100 if c.get("complex_plddt", plddt) is not None else None,
            "complex_iplddt": c.get("complex_iplddt"), "complex_pde": c.get("complex_pde"), "complex_ipde": c.get("complex_ipde"),
            "ipdae_mean_pae": interface_pae(d, model), "plddt_mean": plddt, "rmsd_to_rfd3": rmsd, "rmsd_n_ca": nali,
            "ipSAE": ips.get("ipsae"), "ipSAE_d0chn": ips.get("ipsae_d0chn"), "ipSAE_d0dom": ips.get("ipsae_d0dom"),
            "ipTM_af": ips.get("iptm_af"), "ipTM_d0chn": ips.get("iptm_d0chn"), "pDockQ": ips.get("pdockq"), "pDockQ2": ips.get("pdockq2"), "LIS": ips.get("lis"), "tm_score": tm[0] if tm else None, "tm_rmsd": tm[1] if tm else None}


def write_svg(path, metric, x, y):
    w, h, left, top, right, bottom = 700, 560, 80, 45, 35, 75
    lo, hi = min(x + y), max(x + y)
    if lo == hi: lo, hi = lo - 1, hi + 1
    pad = (hi - lo) * .05; lo -= pad; hi += pad
    def px(v): return left + (v - lo) / (hi - lo) * (w-left-right)
    def py(v): return h-bottom - (v - lo) / (hi-lo) * (h-top-bottom)
    points = " ".join(f"{px(a):.1f},{py(b):.1f}" for a,b in zip(x,y))
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}">', '<rect width="100%" height="100%" fill="white"/>',
           f'<line x1="{px(lo)}" y1="{py(lo)}" x2="{px(hi)}" y2="{py(hi)}" stroke="#888" stroke-dasharray="5,4"/>',
           f'<line x1="{left}" y1="{h-bottom}" x2="{w-right}" y2="{h-bottom}" stroke="black"/><line x1="{left}" y1="{top}" x2="{left}" y2="{h-bottom}" stroke="black"/>',
           f'<polyline points="{points}" fill="none" stroke="#2563eb" stroke-width="0"/>']
    svg += [f'<circle cx="{px(a):.1f}" cy="{py(b):.1f}" r="3" fill="#2563eb" fill-opacity=".55"/>' for a,b in zip(x,y)]
    svg += [f'<text x="{w/2}" y="25" text-anchor="middle" font-family="sans-serif" font-size="16">DD1 Boltz2: {metric} (n={len(x)})</text>',
            f'<text x="{w/2}" y="{h-20}" text-anchor="middle" font-family="sans-serif">{metric} — Mouse EGFR</text>',
            f'<text x="18" y="{h/2}" text-anchor="middle" transform="rotate(-90 18 {h/2})" font-family="sans-serif">{metric} — Human EGFR</text>', '</svg>']
    path.write_text("\n".join(svg))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parent / "out")
    args = ap.parse_args(); args.out.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(SEQ_REPORT.open()))
    model_index = index_models()
    print(f"Indexed {len(model_index)} Boltz2 models by chain sequences.")
    data = []
    for i, row in enumerate(rows, 1):
        for target in ("Mouse", "Human"):
            data.append(collect(row, target, i, model_index))
    fields = list(data[0])
    with (args.out / "boltz2_dd1_metrics_long.csv").open("w", newline="") as h:
        w=csv.DictWriter(h, fieldnames=fields); w.writeheader(); w.writerows(data)
    ok = [r for r in data if r["status"] == "ok"]
    by = {(r["design_id"], r["target"]): r for r in ok}
    paired=[]
    for i, row in enumerate(rows, 1):
        a,b=by.get((row["design_id"],"Mouse")),by.get((row["design_id"],"Human"))
        if not a or not b: continue
        out={"design_id":row["design_id"],"source_run":row["run"]}
        for k in fields:
            if k in {"design_id","source_run","target","status","boltz_dir"}: continue
            out[f"{k}_mouse"]=a.get(k); out[f"{k}_human"]=b.get(k)
        paired.append(out)
    pf=args.out/"boltz2_dd1_metrics_paired.csv"
    with pf.open("w",newline="") as h:
        w=csv.DictWriter(h,fieldnames=list(paired[0])); w.writeheader();w.writerows(paired)
    metrics=["ipSAE","ipTM_af","iptm","confidence_score","complex_plddt","complex_iplddt","complex_ipde","ipdae_mean_pae","plddt_mean","pDockQ","pDockQ2","LIS","rmsd_to_rfd3","tm_score"]
    for metric in metrics:
        x,y=[],[]
        for r in paired:
            try:
                xv=float(r[f"{metric}_mouse"]); yv=float(r[f"{metric}_human"])
                if np.isfinite(xv) and np.isfinite(yv): x.append(xv); y.append(yv)
            except (TypeError,ValueError): pass
        if not x: continue
        write_svg(args.out/f"scatter_{metric}.svg", metric, x, y)
    print(f"Wrote {len(data)} target-specific rows, {len(paired)} paired designs and plots to {args.out}")

if __name__ == "__main__": main()
