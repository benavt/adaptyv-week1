#!/usr/bin/env python3
"""Keep the top X percent of each Foundry run on ipTM, pTM, plDDT, TM, and RMSD.

Each metric is ranked inside its run. The published top 0.5% table used this
rule for ipTM, TM, and RMSD: keep ceil(n * percent / 100) models, and keep
ties with that cutoff. This script adds pTM and plDDT the same way. A model
that qualifies on more than one metric is written once per metric.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
RUNS = REPO / "design" / "rfd3" / "foundry" / "runs"
RFD3_OUT = REPO / "analysis" / "outputs" / "rfd3"
PUBLISHED = RFD3_OUT / "top_0.5pct_sequences.csv"

# Display name, scores.csv column, output column, higher-is-better.
METRICS = (
    ("ipTM", "iptm", "iptm", True),
    ("pTM", "ptm", "ptm", True),
    ("plDDT", "overall_plddt", "plddt", True),
    ("TM", "tm_score_1", "tm", True),
    ("RMSD", "rmsd", "rmsd", False),
)
OUTPUT_COLUMNS = (
    "run",
    "n",
    "metric",
    "rank",
    "model_id",
    "sequence",
    "iptm",
    "ptm",
    "plddt",
    "tm",
    "rmsd",
    "lists",
    "cutoff",
)


def percent_label(percent: float) -> str:
    return f"{percent:.6f}".rstrip("0").rstrip(".")


def default_output(percent: float) -> Path:
    label = percent_label(percent)
    if label == "0.5":
        return RFD3_OUT / "top_0.5pct_iptm_ptm_plddt_tm_rmsd.csv"
    path = RFD3_OUT / f"top_{label}pct_sequences.csv"
    if path.name == PUBLISHED.name:
        return RFD3_OUT / "top_0.5pct_iptm_ptm_plddt_tm_rmsd.csv"
    return path


def model_id_from_rf3(path: str) -> str:
    name = Path(path).name
    if name.endswith("_model.cif"):
        return name[: -len("_model.cif")]
    return Path(name).stem


def load_sequences(run_dir: Path) -> dict[str, str]:
    """Last FASTA record for each model id. Duplicate headers keep the last one."""
    sequences: dict[str, str] = {}
    mpnn = run_dir / "mpnn" / "outputs"
    if not mpnn.is_dir():
        return sequences
    for fasta in sorted(mpnn.glob("*.fa")):
        header: str | None = None
        chunks: list[str] = []

        def commit() -> None:
            if header is None:
                return
            model_id = header.split(",", 1)[0].strip()
            sequence = "".join(chunks).replace(" ", "").upper()
            if model_id and sequence:
                sequences[model_id] = sequence

        for raw in fasta.read_text(encoding="utf-8").splitlines():
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                commit()
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)
        commit()
    return sequences


def load_scores(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"No header in {path}")
        missing = [column for _name, column, _out, _higher in METRICS if column not in reader.fieldnames]
        if "rf3_model" not in reader.fieldnames:
            missing.append("rf3_model")
        if missing:
            raise ValueError(f"{path} is missing columns: {', '.join(missing)}")
        return list(reader)


def _float(text: str) -> float | None:
    try:
        value = float(text)
    except (TypeError, ValueError):
        return None
    if math.isnan(value) or math.isinf(value):
        return None
    return value


def fmt(value: float) -> str:
    return f"{value:.4f}"


def select_metric(
    rows: list[dict[str, str]],
    column: str,
    higher: bool,
    k: int,
) -> tuple[list[tuple[int, dict[str, str], float]], float | None]:
    """Return (row index, row, value) at least as good as the k-th value."""
    parsed: list[tuple[int, float]] = []
    for index, row in enumerate(rows):
        value = _float(row.get(column) or "")
        if value is None:
            continue
        parsed.append((index, value))
    if not parsed:
        return [], None
    parsed.sort(key=lambda item: item[1], reverse=higher)
    cutoff = parsed[min(k, len(parsed)) - 1][1]
    chosen: list[tuple[int, dict[str, str], float]] = []
    for index, value in parsed:
        if higher and value < cutoff:
            break
        if not higher and value > cutoff:
            break
        chosen.append((index, rows[index], value))
    return chosen, cutoff


def competition_rank(value: float, values: list[float], higher: bool) -> int:
    if higher:
        better = sum(other > value for other in values)
    else:
        better = sum(other < value for other in values)
    return better + 1


def metric_values(row: dict[str, str]) -> dict[str, str]:
    formatted: dict[str, str] = {}
    for _display, source, output, _higher in METRICS:
        value = _float(row.get(source) or "")
        formatted[output] = fmt(value) if value is not None else ""
    return formatted


def filter_run(run_dir: Path, percent: float) -> list[dict[str, str]]:
    scores_path = run_dir / "scores.csv"
    rows = load_scores(scores_path)
    n = len(rows)
    if n == 0:
        return []
    k = max(1, math.ceil(n * percent / 100.0))
    sequences = load_sequences(run_dir)
    lists_by_model: dict[str, list[str]] = {}
    pending: list[tuple[str, str, dict[str, str], float, list[float], float, bool]] = []
    for display, column, _output, higher in METRICS:
        chosen, cutoff = select_metric(rows, column, higher, k)
        if cutoff is None:
            continue
        all_values = [
            parsed
            for row in rows
            if (parsed := _float(row.get(column) or "")) is not None
        ]
        for _index, row, value in chosen:
            model_id = model_id_from_rf3(row["rf3_model"])
            names = lists_by_model.setdefault(model_id, [])
            if display not in names:
                names.append(display)
            pending.append((display, model_id, row, value, all_values, cutoff, higher))
    missing_sequence = 0
    selected: list[dict[str, str]] = []
    for display, model_id, row, value, all_values, cutoff, higher in pending:
        sequence = sequences.get(model_id, "")
        if not sequence:
            missing_sequence += 1
        selected.append(
            {
                "run": run_dir.name,
                "n": str(n),
                "metric": display,
                "rank": str(competition_rank(value, all_values, higher)),
                "model_id": model_id,
                "sequence": sequence,
                "lists": ";".join(lists_by_model[model_id]),
                "cutoff": fmt(cutoff),
                **metric_values(row),
            }
        )
    if missing_sequence:
        print(
            f"  {run_dir.name}: {missing_sequence} selected rows have no FASTA sequence",
            file=sys.stderr,
        )
    return selected


def write_table(path: Path, rows: list[dict[str, str]]) -> None:
    if path.resolve() == PUBLISHED.resolve():
        raise SystemExit(f"Refusing to overwrite {PUBLISHED}")
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    tmp.replace(path)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--percent",
        type=float,
        default=0.5,
        help="Percent of each run to keep per metric (default: 0.5).",
    )
    parser.add_argument(
        "--runs",
        type=Path,
        default=RUNS,
        help="Foundry runs directory.",
    )
    parser.add_argument(
        "--run",
        action="append",
        default=[],
        help="Limit to this run directory name. Repeat to select several.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Output CSV. Defaults under analysis/outputs/rfd3 and does not replace the published top 0.5%% table.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if not math.isfinite(args.percent) or args.percent <= 0 or args.percent > 100:
        print("--percent must be between 0 and 100", file=sys.stderr)
        return 2
    runs_dir = args.runs.expanduser().resolve()
    if not runs_dir.is_dir():
        print(f"Runs directory not found: {runs_dir}", file=sys.stderr)
        return 2
    output = (args.output or default_output(args.percent)).expanduser().resolve()
    wanted = set(args.run)
    selected_rows: list[dict[str, str]] = []
    seen_runs = 0
    for run_dir in sorted(path for path in runs_dir.iterdir() if path.is_dir()):
        if wanted and run_dir.name not in wanted:
            continue
        scores = run_dir / "scores.csv"
        if not scores.is_file():
            continue
        seen_runs += 1
        rows = filter_run(run_dir, args.percent)
        print(f"{run_dir.name}: {len(rows)} rows", flush=True)
        selected_rows.extend(rows)
    if wanted:
        scanned = {
            path.name
            for path in runs_dir.iterdir()
            if path.is_dir() and (path / "scores.csv").is_file()
        }
        absent = sorted(wanted - scanned)
        if absent:
            print(f"No scores.csv for: {', '.join(absent)}", file=sys.stderr)
            return 2
    if seen_runs == 0:
        print(f"No scores.csv under {runs_dir}", file=sys.stderr)
        return 2
    write_table(output, selected_rows)
    print(f"Wrote {len(selected_rows)} rows from {seen_runs} runs to {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
