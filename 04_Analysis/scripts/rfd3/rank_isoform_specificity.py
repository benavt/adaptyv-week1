#!/usr/bin/env python3
"""Rank designs by isoform specificity from a wide DeltaForge CSV.

A row is kept only when the --on receptor is tighter than --off (lower Kd and
delta-G) and more confident (higher mean plDDT, ipTM, pTM, and peptide ipSAE).
Kept rows are ordered from most specific to least by the log10 Kd ratio, then
the delta-G gap, then the sum of those four confidence gaps.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
RFD3_OUT = REPO / "analysis" / "outputs" / "rfd3"

LOWER_IS_BETTER = ("kd_nm", "delta_g")
HIGHER_IS_BETTER = ("mean_plddt", "iptm", "ptm", "peptide_ipsae")
RANK_COLUMNS = (
    "specificity_rank",
    "log10_kd_ratio",
    "delta_delta_g",
    "iptm_gap",
    "ptm_gap",
    "mean_plddt_gap",
    "peptide_ipsae_gap",
)
PROTECTED_NAMES = {
    "top_0.5pct_sequences.csv",
    "top_0.5pct_sequences_deltaforge.csv",
}


def prefixes(fieldnames: list[str]) -> set[str]:
    found: set[str] = set()
    suffixes = LOWER_IS_BETTER + HIGHER_IS_BETTER
    for name in fieldnames:
        folded = name.casefold()
        for field in suffixes:
            suffix = "_" + field
            if folded.endswith(suffix):
                found.add(name[: -len(suffix)])
                break
    return found


def resolve_prefix(fieldnames: list[str], requested: str) -> str:
    available = prefixes(fieldnames)
    requested_cf = requested.casefold()
    for prefix in available:
        if prefix.casefold() == requested_cf:
            return prefix
    names = ", ".join(sorted(available)) or "(none)"
    raise SystemExit(f"No DeltaForge columns for {requested}. Prefixes in this file: {names}")


def cell(row: dict[str, str], prefix: str, field: str) -> float | None:
    text = (row.get(f"{prefix}_{field}") or "").strip()
    if not text:
        return None
    try:
        value = float(text)
    except ValueError:
        return None
    if not math.isfinite(value):
        return None
    return value


def fmt(value: float) -> str:
    return f"{value:.6g}"


def evaluate(row: dict[str, str], on: str, off: str) -> dict[str, float] | None:
    values: dict[str, float] = {}
    for field in LOWER_IS_BETTER + HIGHER_IS_BETTER:
        on_value = cell(row, on, field)
        off_value = cell(row, off, field)
        if on_value is None or off_value is None:
            return None
        values[f"on_{field}"] = on_value
        values[f"off_{field}"] = off_value
    if values["on_kd_nm"] <= 0 or values["off_kd_nm"] <= 0:
        return None
    if not (values["on_kd_nm"] < values["off_kd_nm"] and values["on_delta_g"] < values["off_delta_g"]):
        return None
    for field in HIGHER_IS_BETTER:
        if not values[f"on_{field}"] > values[f"off_{field}"]:
            return None
    values["log10_kd_ratio"] = math.log10(values["off_kd_nm"] / values["on_kd_nm"])
    values["delta_delta_g"] = values["off_delta_g"] - values["on_delta_g"]
    values["iptm_gap"] = values["on_iptm"] - values["off_iptm"]
    values["ptm_gap"] = values["on_ptm"] - values["off_ptm"]
    values["mean_plddt_gap"] = values["on_mean_plddt"] - values["off_mean_plddt"]
    values["peptide_ipsae_gap"] = values["on_peptide_ipsae"] - values["off_peptide_ipsae"]
    values["confidence_gap"] = (
        values["iptm_gap"] + values["ptm_gap"] + values["mean_plddt_gap"] + values["peptide_ipsae_gap"]
    )
    return values


def default_output(input_path: Path, on: str, off: str) -> Path:
    return RFD3_OUT / f"{input_path.stem}_{on}_over_{off}.csv"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="Wide DeltaForge CSV.")
    parser.add_argument("--on", default="Brd4ET", help="Preferred receptor prefix (default: Brd4ET).")
    parser.add_argument("--off", default="Brd3ET", help="Comparison receptor prefix (default: Brd3ET).")
    parser.add_argument("--output", type=Path, help="Ranked CSV. Defaults next to the input.")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    input_path = args.input.expanduser().resolve()
    if not input_path.is_file():
        print(f"Input CSV not found: {input_path}", file=sys.stderr)
        return 2
    output_path = (args.output or default_output(input_path, args.on, args.off)).expanduser().resolve()
    if output_path.name in PROTECTED_NAMES or output_path == input_path:
        print(f"Refusing to overwrite {output_path}", file=sys.stderr)
        return 2

    with input_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            print(f"No header in {input_path}", file=sys.stderr)
            return 2
        fieldnames = list(reader.fieldnames)
        rows = list(reader)

    on = resolve_prefix(fieldnames, args.on)
    off = resolve_prefix(fieldnames, args.off)
    if on == off:
        print("--on and --off resolved to the same receptor", file=sys.stderr)
        return 2

    kept: list[tuple[dict[str, str], dict[str, float]]] = []
    for row in rows:
        scores = evaluate(row, on, off)
        if scores is not None:
            kept.append((row, scores))
    kept.sort(
        key=lambda item: (
            -item[1]["log10_kd_ratio"],
            -item[1]["delta_delta_g"],
            -item[1]["confidence_gap"],
        )
    )

    out_fields = list(RANK_COLUMNS) + [name for name in fieldnames if name not in RANK_COLUMNS]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = output_path.with_suffix(output_path.suffix + ".tmp")
    with tmp.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=out_fields, extrasaction="ignore")
        writer.writeheader()
        for rank, (row, scores) in enumerate(kept, start=1):
            out = dict(row)
            out["specificity_rank"] = str(rank)
            for name in RANK_COLUMNS[1:]:
                out[name] = fmt(scores[name])
            writer.writerow(out)
    tmp.replace(output_path)
    print(
        f"Kept {len(kept)} of {len(rows)} rows ({on} over {off}). Wrote {output_path}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
