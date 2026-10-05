#!/usr/bin/env python3
"""DeltaForge-score top_0.5pct sequences against Brd4ET and Brd3ET.

Fusion rows that end with the BRD4 ET sequence are scored as the binder only
(the target suffix is stripped). Other rows are scored as complete binders.
Each unique sequence is folded and scored once per target. Results checkpoint
after every job so a re-run skips successes. The companion CSV has one row per
input row.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
LIGAND_AI = REPO / "design" / "ligand_ai"
SEQUENCES = REPO / "staging" / "sequences"
RFD3_OUT = REPO / "analysis" / "outputs" / "rfd3"
CHECKPOINT_NAME = "top_0.5pct_sequences_deltaforge_checkpoint.json"
for path in (LIGAND_AI, HERE):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from score_outputs_deltaforge import (  # noqa: E402
    _enrich_job_error,
    load_target_sequence,
    score_with_retries,
)
from src.ligandai_local.client import LigandAIClient  # noqa: E402

SOURCE_COLUMNS = (
    "run",
    "n",
    "metric",
    "rank",
    "model_id",
    "sequence",
    "iptm",
    "tm",
    "rmsd",
    "lists",
    "cutoff",
)
METRIC_FIELDS = (
    "delta_g",
    "kd_nm",
    "classification",
    "predicted_affinity_tier",
    "predicted_binder",
    "predicted_binder_call",
    "predicted_binder_probability",
    "affinity_quotable",
    "iptm",
    "ptm",
    "ipsae",
    "peptide_ipsae",
    "mean_plddt",
    "scorer",
    "job_id",
    "error",
)
_DROP_RESULT_KEYS = {"pdbData", "pdb", "pdb_data", "pdbContent", "cif", "cifData", "cifContent"}


def binder_from_sequence(sequence: str, brd4_sequence: str) -> str:
    if sequence.endswith(brd4_sequence):
        binder = sequence[: -len(brd4_sequence)]
        if not binder:
            raise ValueError("Empty binder after stripping BRD4 ET suffix")
        return binder
    return sequence


def load_input_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"No header in {path}")
        missing = [name for name in SOURCE_COLUMNS if name not in reader.fieldnames]
        if missing:
            raise ValueError(f"{path} is missing columns: {', '.join(missing)}")
        return list(reader)


def unique_sequences(rows: list[dict[str, str]]) -> list[tuple[str, str]]:
    """Return (sequence, first model_id) in first-seen order."""
    seen: dict[str, str] = {}
    ordered: list[tuple[str, str]] = []
    for row in rows:
        sequence = (row.get("sequence") or "").strip()
        if not sequence or sequence in seen:
            continue
        model_id = (row.get("model_id") or "").strip() or "binder"
        seen[sequence] = model_id
        ordered.append((sequence, model_id))
    return ordered


def model_run_map(rows: list[dict[str, str]]) -> dict[str, list[str]]:
    mapping: dict[str, list[str]] = {}
    for row in rows:
        model_id = (row.get("model_id") or "").strip()
        run = (row.get("run") or "").strip()
        if not model_id or not run:
            continue
        runs = mapping.setdefault(model_id, [])
        if run not in runs:
            runs.append(run)
    return mapping


def load_checkpoint(path: Path) -> list[dict[str, Any]]:
    if path.is_dir():
        rows: list[dict[str, Any]] = []
        for child in sorted(path.glob(f"*/{CHECKPOINT_NAME}")):
            rows.extend(load_checkpoint(child))
        unassigned = path / "unassigned_checkpoint.json"
        if unassigned.is_file():
            rows.extend(load_checkpoint(unassigned))
        return rows
    if not path.exists():
        return []
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, list):
        raise ValueError(f"Expected a JSON array in {path}")
    return [row for row in data if isinstance(row, dict)]


def write_checkpoints(
    directory: Path,
    rows: list[dict[str, Any]],
    runs_by_model: dict[str, list[str]],
) -> None:
    grouped: dict[str, list[dict[str, Any]]] = {}
    unassigned: list[dict[str, Any]] = []
    for row in rows:
        model_id = row.get("model_id")
        run_names = runs_by_model.get(model_id, []) if isinstance(model_id, str) else []
        if not run_names:
            unassigned.append(row)
            continue
        for run in run_names:
            grouped.setdefault(run, []).append(row)
    for run, group in grouped.items():
        write_json(directory / run / CHECKPOINT_NAME, group)
    if unassigned:
        write_json(directory / "unassigned_checkpoint.json", unassigned)


def is_success(row: dict[str, Any] | None) -> bool:
    return bool(row) and row.get("result") is not None and not row.get("error")


def compact_result(result: Any) -> Any:
    if not isinstance(result, dict):
        return result
    compact: dict[str, Any] = {}
    for key, value in result.items():
        if key in _DROP_RESULT_KEYS and isinstance(value, str):
            compact[key] = f"<omitted {len(value)} chars>"
        else:
            compact[key] = value
    return compact


def write_json(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(rows, indent=2, default=str) + "\n", encoding="utf-8")
    tmp.replace(path)


def _cell(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    return str(value)


def extract_metrics(row: dict[str, Any] | None) -> dict[str, str]:
    result = (row or {}).get("result") or {}
    if not isinstance(result, dict):
        result = {}
    deltaforge = result.get("deltaforge") or {}
    if not isinstance(deltaforge, dict):
        deltaforge = {}
    fold = result.get("fold") or {}
    if not isinstance(fold, dict):
        fold = {}
    values: dict[str, Any] = {
        "delta_g": deltaforge.get("delta_g"),
        "kd_nm": deltaforge.get("kd_nm"),
        "classification": deltaforge.get("classification"),
        "predicted_affinity_tier": deltaforge.get("predicted_affinity_tier"),
        "predicted_binder": deltaforge.get("predicted_binder"),
        "predicted_binder_call": deltaforge.get("predicted_binder_call"),
        "predicted_binder_probability": deltaforge.get("predicted_binder_probability"),
        "affinity_quotable": deltaforge.get("affinity_quotable"),
        "iptm": fold.get("iptm"),
        "ptm": fold.get("ptm"),
        "ipsae": fold.get("ipsae", result.get("ipsae")),
        "peptide_ipsae": fold.get("peptide_ipsae"),
        "mean_plddt": fold.get("mean_plddt"),
        "scorer": deltaforge.get("scorer"),
        "job_id": result.get("jobId") or result.get("job_id"),
        "error": (row or {}).get("error") or "",
    }
    return {key: _cell(values[key]) for key in METRIC_FIELDS}


def output_fieldnames() -> list[str]:
    names = list(SOURCE_COLUMNS) + ["binder_sequence"]
    for prefix in ("brd4", "brd3"):
        names.extend(f"{prefix}_{field}" for field in METRIC_FIELDS)
    return names


def write_output_csv(
    path: Path,
    input_rows: list[dict[str, str]],
    by_key: dict[tuple[str, str], dict[str, Any]],
    brd4_sequence: str,
) -> None:
    fieldnames = output_fieldnames()
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in input_rows:
            sequence = (row.get("sequence") or "").strip()
            try:
                binder = binder_from_sequence(sequence, brd4_sequence) if sequence else ""
            except ValueError:
                binder = ""
            out = {name: row.get(name, "") for name in SOURCE_COLUMNS}
            out["binder_sequence"] = binder
            for prefix in ("brd4", "brd3"):
                metrics = extract_metrics(by_key.get((sequence, prefix)))
                for field, value in metrics.items():
                    out[f"{prefix}_{field}"] = value
            writer.writerow(out)
    tmp.replace(path)


def score_pair(
    *,
    sequence: str,
    model_id: str,
    binder_sequence: str,
    target_key: str,
    target_name: str,
    target_sequence: str,
    scorer: str,
    max_attempts: int,
    backoff_s: float,
) -> dict[str, Any]:
    client = LigandAIClient()
    label = f"{model_id} vs {target_name}"
    print(
        f"  scoring {label} (binder_len={len(binder_sequence)})",
        flush=True,
    )
    try:
        result = score_with_retries(
            client,
            binder_sequence=binder_sequence,
            target_sequence=target_sequence,
            binder_name=model_id,
            target_name=target_name,
            scorer=scorer,
            max_attempts=max_attempts,
            backoff_s=backoff_s,
        )
        print(f"  ok: {label}", flush=True)
        return {
            "sequence": sequence,
            "model_id": model_id,
            "binder_sequence": binder_sequence,
            "target_key": target_key,
            "target_name": target_name,
            "result": compact_result(result),
            "error": None,
        }
    except Exception as exc:  # noqa: BLE001
        detail = _enrich_job_error(client, exc)
        print(f"  ERROR: {label}: {detail}", flush=True)
        return {
            "sequence": sequence,
            "model_id": model_id,
            "binder_sequence": binder_sequence,
            "target_key": target_key,
            "target_name": target_name,
            "result": None,
            "error": detail,
        }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="DeltaForge-score top_0.5pct sequences against Brd4ET and Brd3ET"
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=RFD3_OUT / "top_0.5pct_sequences.csv",
    )
    parser.add_argument(
        "--brd4-fasta",
        type=Path,
        default=SEQUENCES / "Brd4ET.fasta",
    )
    parser.add_argument(
        "--brd3-fasta",
        type=Path,
        default=SEQUENCES / "Brd3ET.fasta",
    )
    parser.add_argument(
        "--checkpoint",
        type=Path,
        default=RFD3_OUT,
        help=(
            "Directory of per-run checkpoint JSON files, or one JSON file. "
            "A directory is loaded and saved as "
            "<run>/top_0.5pct_sequences_deltaforge_checkpoint.json "
            "using the CSV run column."
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=RFD3_OUT / "top_0.5pct_sequences_deltaforge.csv",
    )
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--scorer", default="auto")
    parser.add_argument("--max-attempts", type=int, default=8)
    parser.add_argument("--backoff-s", type=float, default=30.0)
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-score pairs that already have a successful checkpoint result.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    input_path = args.input.expanduser().resolve()
    brd4_fasta = args.brd4_fasta.expanduser().resolve()
    brd3_fasta = args.brd3_fasta.expanduser().resolve()
    checkpoint_path = args.checkpoint.expanduser().resolve()
    output_path = args.output.expanduser().resolve()

    for path, label in (
        (input_path, "Input CSV"),
        (brd4_fasta, "Brd4 FASTA"),
        (brd3_fasta, "Brd3 FASTA"),
    ):
        if not path.exists():
            print(f"{label} not found: {path}", file=sys.stderr)
            return 2

    brd4_name, brd4_sequence = load_target_sequence(brd4_fasta)
    brd3_name, brd3_sequence = load_target_sequence(brd3_fasta)
    input_rows = load_input_rows(input_path)
    runs_by_model = model_run_map(input_rows)
    sequences = unique_sequences(input_rows)
    targets = (
        ("brd4", brd4_name, brd4_sequence),
        ("brd3", brd3_name, brd3_sequence),
    )

    checkpoint_rows = [] if args.force else load_checkpoint(checkpoint_path)
    by_key: dict[tuple[str, str], dict[str, Any]] = {}
    for row in checkpoint_rows:
        sequence = row.get("sequence")
        target_key = row.get("target_key")
        if isinstance(sequence, str) and isinstance(target_key, str):
            by_key[(sequence, target_key)] = row

    pending: list[dict[str, Any]] = []
    for sequence, model_id in sequences:
        binder = binder_from_sequence(sequence, brd4_sequence)
        for target_key, target_name, target_sequence in targets:
            if is_success(by_key.get((sequence, target_key))):
                continue
            pending.append(
                {
                    "sequence": sequence,
                    "model_id": model_id,
                    "binder_sequence": binder,
                    "target_key": target_key,
                    "target_name": target_name,
                    "target_sequence": target_sequence,
                }
            )

    cached = len(sequences) * len(targets) - len(pending)
    print(
        f"Sequences: {len(input_rows)} rows, {len(sequences)} unique; "
        f"targets {brd4_name} ({len(brd4_sequence)} aa), "
        f"{brd3_name} ({len(brd3_sequence)} aa)",
        flush=True,
    )
    print(
        f"Jobs: {len(pending)} to score, {cached} cached; workers={args.workers}",
        flush=True,
    )

    lock = threading.Lock()
    scored = 0
    failed = 0

    def persist() -> None:
        ordered = [
            by_key[(sequence, target_key)]
            for sequence, _model_id in sequences
            for target_key, _name, _seq in targets
            if (sequence, target_key) in by_key
        ]
        if checkpoint_path.is_dir():
            write_checkpoints(checkpoint_path, ordered, runs_by_model)
        else:
            write_json(checkpoint_path, ordered)
        write_output_csv(output_path, input_rows, by_key, brd4_sequence)

    if pending:
        workers = max(1, int(args.workers))
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = [
                pool.submit(
                    score_pair,
                    sequence=job["sequence"],
                    model_id=job["model_id"],
                    binder_sequence=job["binder_sequence"],
                    target_key=job["target_key"],
                    target_name=job["target_name"],
                    target_sequence=job["target_sequence"],
                    scorer=args.scorer,
                    max_attempts=args.max_attempts,
                    backoff_s=args.backoff_s,
                )
                for job in pending
            ]
            for future in as_completed(futures):
                row = future.result()
                key = (row["sequence"], row["target_key"])
                with lock:
                    by_key[key] = row
                    if row.get("error"):
                        failed += 1
                    else:
                        scored += 1
                    persist()
                    done = scored + failed
                    print(
                        f"  progress {done}/{len(pending)} scored={scored} failed={failed}",
                        flush=True,
                    )
    else:
        persist()

    print(
        f"\nDone: scored={scored} failed={failed} cached={cached}",
        flush=True,
    )
    print(f"Checkpoint: {checkpoint_path}", flush=True)
    print(f"CSV: {output_path}", flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
