#!/usr/bin/env python3
"""Count DD1 ASP residues near EGFR-D3 HIS37 and HIS100 in RFD3 models."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
from pathlib import Path

from Bio.PDB import MMCIFParser

ROOT = Path(__file__).resolve().parents[2]
RUNS = {
    "mouse": ROOT / "02_Design/foundry/runs/DD1_Mouse_EGFR_hotspots_09_29/rfd3/outputs",
    "human": ROOT / "02_Design/foundry/runs/DD1_Human_EGFR_hotspots_09_29/rfd3/outputs",
}
AA = dict(zip("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(), "ARNDCQEGHILKMFPSTWYV"))
AA.update({"HID": "H", "HIE": "H", "HIP": "H", "ASH": "D", "GLH": "E", "CYX": "C", "MSE": "M"})


@lru_cache(maxsize=1)
def egfr_target_sequences():
    """Load the recorded D3 sequences, never infer target species from a filename."""
    targets = {}
    for species in ("Human", "Mouse"):
        path = ROOT / f"01_Staging/EGFR/{species}_EGFR_310_501_afs3.json"
        payload = json.loads(path.read_text())
        proteins = [item["protein"]["sequence"] for item in payload["sequences"] if "protein" in item]
        if len(proteins) != 1:
            raise ValueError(f"Expected one EGFR protein sequence in {path}")
        targets[species.lower()] = proteins[0]
    return targets


def target_species(receptor, targets):
    """Require one exact species match to the actual receptor sequence."""
    sequence = "".join(AA.get(r.resname, "X") for r in receptor if r.id[0] == " ")
    matches = [species for species, expected in targets.items() if sequence == expected]
    if len(matches) != 1:
        raise ValueError(f"EGFR chain {receptor.id} sequence matches {len(matches)} staged targets; cannot assign species")
    return matches[0]


def atoms_within(a, b, cutoff_sq):
    return any(sum((x - y) ** 2 for x, y in zip(xyz_a.coord, xyz_b.coord)) <= cutoff_sq
               for xyz_a in a for xyz_b in b)


def analyze(path: Path, cutoff: float) -> dict:
    structure = MMCIFParser(QUIET=True).get_structure(path.stem, gzip.open(path, "rt"))
    chains = list(structure[0].get_chains())
    # RFD3 output maps the designed chain to A and the fixed EGFR chain to B.
    # Fall back to the longest chain as EGFR for robustness to future output naming.
    receptor = next((c for c in chains if len(list(c.get_residues())) >= 190), max(chains, key=lambda c: len(list(c.get_residues()))))
    binder = next(c for c in chains if c.id != receptor.id)
    his = {n: next((r for r in receptor if r.id[1] == n and r.resname == "HIS"), None) for n in (37, 100)}
    asps = [r for r in binder if r.resname == "ASP"]
    cutoff_sq = cutoff * cutoff
    near = {}
    for n, residue in his.items():
        near[n] = [r.id[1] for r in asps if residue is not None and atoms_within(list(r.get_atoms()), list(residue.get_atoms()), cutoff_sq)]
    return {
        "design_id": path.name.removesuffix(".cif.gz"),
        "structure_path": str(path),
        "receptor_chain": receptor.id,
        "binder_chain": binder.id,
        "binder_asp_count": len(asps),
        "asp_near_his37": ";".join(map(str, near[37])),
        "asp_near_his100": ";".join(map(str, near[100])),
        "has_asp_near_his37": bool(near[37]),
        "has_asp_near_his100": bool(near[100]),
        "has_asp_near_either": bool(near[37] or near[100]),
        "has_asp_near_both": bool(near[37] and near[100]),
        "cutoff_angstrom": cutoff,
    }


def analyze_refolded(path: Path, cutoff: float, targets: dict | None = None) -> dict:
    """Analyze a Boltz2 refolded DD1/EGFR complex (binder=A, EGFR=B)."""
    structure = MMCIFParser(QUIET=True).get_structure(path.stem, str(path))
    model = next(structure.get_models())
    chains = list(model.get_chains())
    receptor = next((c for c in chains if c.id == "B"), max(chains, key=lambda c: len(list(c.get_residues()))))
    binder = next(c for c in chains if c.id != receptor.id)
    species = target_species(receptor, egfr_target_sequences() if targets is None else targets)
    his = {n: next((r for r in receptor if r.id[1] == n and r.resname == "HIS"), None) for n in (37, 100)}
    asps = [r for r in binder if r.resname == "ASP"]
    cutoff_sq = cutoff * cutoff
    near = {n: [r.id[1] for r in asps if residue is not None and atoms_within(list(r.get_atoms()), list(residue.get_atoms()), cutoff_sq)] for n, residue in his.items()}
    name = path.parent.name
    return {
        "species": species, "design_id": name.split("_", 2)[0] + "_" + name.split("_", 2)[1] if name.startswith("DD1_") else name,
        "model_id": name, "structure_path": str(path), "receptor_chain": receptor.id, "binder_chain": binder.id,
        "binder_asp_count": len(asps), "asp_near_his37": ";".join(map(str, near[37])), "asp_near_his100": ";".join(map(str, near[100])),
        "has_asp_near_his37": bool(near[37]), "has_asp_near_his100": bool(near[100]), "has_asp_near_either": bool(near[37] or near[100]),
        "has_asp_near_both": bool(near[37] and near[100]), "cutoff_angstrom": cutoff,
    }


def _analyze_job(job):
    """Pickle-friendly worker wrapper."""
    _, path, cutoff = job
    return analyze(path, cutoff)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cutoff", type=float, default=5.0, help="Heavy-atom distance cutoff in Å (default: 5.0).")
    parser.add_argument("--out", type=Path, default=Path(__file__).resolve().parent / "out/dd1_asp_his_proximity.csv")
    parser.add_argument("--refolded-root", type=Path, help="Scan Boltz2 refolded complexes instead of RFD3 outputs.")
    args = parser.parse_args()
    if args.refolded_root:
        paths = sorted(args.refolded_root.glob("*/*model_0.cif"))
        rows = [analyze_refolded(path, args.cutoff) for path in paths]
        fields = list(rows[0]) if rows else []
        args.out.parent.mkdir(parents=True, exist_ok=True)
        with args.out.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields); writer.writeheader(); writer.writerows(rows)
        print(f"Analyzed {len(rows)} refolded complexes at {args.cutoff:.1f} Å.")
        for species in ("mouse", "human"):
            subset = [r for r in rows if r["species"] == species]
            print(species, "either=", sum(r["has_asp_near_either"] for r in subset), "HIS37=", sum(r["has_asp_near_his37"] for r in subset), "HIS100=", sum(r["has_asp_near_his100"] for r in subset), "both=", sum(r["has_asp_near_both"] for r in subset))
        return
    jobs = [(run, path, args.cutoff) for run, directory in RUNS.items() for path in sorted(directory.glob("*.cif.gz"))]
    # CIF parsing is CPU-heavy; use worker processes so the full 20k set is practical.
    with ProcessPoolExecutor(max_workers=min(12, os.cpu_count() or 1)) as pool:
        rows = []
        for (run, _, _), row in zip(jobs, pool.map(_analyze_job, jobs, chunksize=25)):
            row["run"] = run
            rows.append(row)
    fields = ["run"] + [k for k in rows[0] if k != "run"]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Analyzed {len(rows)} structures at {args.cutoff:.1f} Å.")
    for run in RUNS:
        subset = [r for r in rows if r["run"] == run]
        print(run, "either=", sum(r["has_asp_near_either"] for r in subset), "HIS37=", sum(r["has_asp_near_his37"] for r in subset), "HIS100=", sum(r["has_asp_near_his100"] for r in subset), "both=", sum(r["has_asp_near_both"] for r in subset))


if __name__ == "__main__":
    main()
