#!/usr/bin/env python3
"""Build FASTA and AlphaFold 3 inputs for the 342 DD1 designs."""
import csv
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
REPORT = ROOT / "03_Filtering/RFD3_filtering/out/dd1_asp_near_both_his_sequences.csv"
EGFR_JSON = {
    "Human": ROOT / "01_Staging/EGFR/Human_EGFR_310_501_afs3.json",
    "Mouse": ROOT / "01_Staging/EGFR/Mouse_EGFR_310_501_afs3.json",
}


def target_sequence(path):
    data = json.loads(path.read_text())
    return next(item["protein"]["sequence"] for item in data["sequences"] if "protein" in item)


def safe(value):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def main():
    fasta_dir = HERE / "FASTA"
    json_dir = HERE / "AF3_JSON"
    fasta_dir.mkdir(parents=True, exist_ok=True)
    json_dir.mkdir(parents=True, exist_ok=True)
    targets = {species: target_sequence(path) for species, path in EGFR_JSON.items()}
    rows = list(csv.DictReader(REPORT.open()))
    for index, row in enumerate(rows, start=1):
        binder = row["binder_sequence"]
        design = safe(row["design_id"])
        for species, target in targets.items():
            stem = f"DD1_{index:03d}_{species}_{design}"
            fasta = f">DD1_binder|{row['run']}|{row['design_id']}\n{binder}\n>EGFR_D3|{species}\n{target}\n"
            (fasta_dir / f"{stem}.fasta").write_text(fasta)
            af3 = {
                "name": stem,
                "modelSeeds": [1],
                "sequences": [
                    {"protein": {"id": "A", "sequence": binder, "description": f"DD1 binder {row['design_id']}"}},
                    {"protein": {"id": "B", "sequence": target, "description": f"{species} EGFR D3 residues 310-501"}},
                ],
                "dialect": "alphafold3",
                "version": 4,
            }
            (json_dir / f"{stem}.json").write_text(json.dumps(af3, indent=2) + "\n")
    print(f"Generated {len(rows) * 2} FASTA and {len(rows) * 2} AF3 JSON files ({len(rows) * 4} total).")


if __name__ == "__main__":
    main()
