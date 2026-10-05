#!/usr/bin/env python3
"""Summarize and conservatively filter RFD3 JSON outputs."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RUNS = {
    "mouse": ROOT / "02_Design/foundry/runs/DD1_Mouse_EGFR_hotspots_09_29/rfd3/outputs",
    "human": ROOT / "02_Design/foundry/runs/DD1_Human_EGFR_hotspots_09_29/rfd3/outputs",
}
METRICS = [
    "max_ca_deviation", "n_chainbreaks",
    "n_clashing.interresidue_clashes_w_sidechain",
    "n_clashing.interresidue_clashes_w_backbone", "non_loop_fraction",
    "loop_fraction", "helix_fraction", "sheet_fraction", "num_ss_elements",
    "radius_of_gyration", "alanine_content", "glycine_content", "num_residues",
]


def design_id(path: Path) -> str:
    return path.stem


def read_rows(run: str, output_dir: Path) -> list[dict]:
    rows = []
    for path in sorted(output_dir.glob("*.json")):
        data = json.loads(path.read_text())
        metrics = data.get("metrics", {})
        row = {"run": run, "design_id": design_id(path), "json_path": str(path)}
        row["structure_path"] = str(path.with_suffix(".cif.gz"))
        row.update({key: metrics.get(key) for key in METRICS})
        rows.append(row)
    return rows


def classify(row: dict, args: argparse.Namespace) -> tuple[bool, str]:
    failures = []
    if row["n_chainbreaks"] is None or row["n_chainbreaks"] > args.max_chainbreaks:
        failures.append("chainbreaks")
    if row["n_clashing.interresidue_clashes_w_sidechain"] is None or row["n_clashing.interresidue_clashes_w_sidechain"] > args.max_sidechain_clashes:
        failures.append("sidechain_clashes")
    if row["n_clashing.interresidue_clashes_w_backbone"] is None or row["n_clashing.interresidue_clashes_w_backbone"] > args.max_backbone_clashes:
        failures.append("backbone_clashes")
    if row["max_ca_deviation"] is None or row["max_ca_deviation"] > args.max_ca_deviation:
        failures.append("ca_deviation")
    return not failures, ";".join(failures) or "pass"


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else []
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=Path(__file__).resolve().parent / "out")
    parser.add_argument("--max-chainbreaks", type=int, default=0)
    parser.add_argument("--max-sidechain-clashes", type=int, default=0)
    parser.add_argument("--max-backbone-clashes", type=int, default=0)
    parser.add_argument("--max-ca-deviation", type=float, default=2.0)
    args = parser.parse_args()

    rows = []
    for run, directory in DEFAULT_RUNS.items():
        rows.extend(read_rows(run, directory))
    for row in rows:
        row["passes_sanity_filter"], row["failure_reasons"] = classify(row, args)

    write_csv(args.out / "summary.csv", rows)
    passing = [row for row in rows if row["passes_sanity_filter"]]
    write_csv(args.out / "passing_designs.csv", passing)

    counts = []
    for run in DEFAULT_RUNS:
        subset = [r for r in rows if r["run"] == run]
        counts.append({"run": run, "total": len(subset), "passing": sum(r["passes_sanity_filter"] for r in subset), "failed": sum(not r["passes_sanity_filter"] for r in subset)})
    write_csv(args.out / "filter_counts.csv", counts)
    print(f"Wrote {len(rows)} designs; {len(passing)} pass the sanity filter.")


if __name__ == "__main__":
    main()
