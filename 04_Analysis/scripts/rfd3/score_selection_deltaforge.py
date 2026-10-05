#!/usr/bin/env python3
"""DeltaForge-score a selected binder table against one or more receptor FASTAs.

Each unique designed sequence is scored once per receptor. Fusion sequences that
end in a receptor sequence are scored as the binder only; the longest matching
suffix is removed. Checkpoint JSON is written per Foundry run under
analysis/outputs/rfd3/<run>/, and a wide CSV is written at the rfd3 output root.
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
REPO = HERE.parents[3]
DELTAFORGE = HERE.parent / "deltaforge"
LIGAND_AI = REPO / "design" / "ligand_ai"
SEQUENCES = REPO / "staging" / "sequences"
RFD3_OUT = REPO / "analysis" / "outputs" / "rfd3"

PROTECTED_NAMES = {
    "top_0.5pct_sequences.csv",
    "top_0.5pct_sequences_deltaforge.csv",
    "top_0.5pct_sequences_deltaforge_checkpoint.json",
    "unassigned_checkpoint.json",
}

REQUIRED_COLUMNS = ("run", "model_id", "sequence")


def refuse_overwrite(path: Path) -> None:
    if path.name in PROTECTED_NAMES:
        raise SystemExit(f"Refusing to overwrite {path}")


def checkpoint_filename(stem: str) -> str:
    return f"{stem}_deltaforge_checkpoint.json"


def unassigned_filename(stem: str) -> str:
    return f"unassigned_{stem}_checkpoint.json"


def summary_filename(stem: str) -> str:
    return f"{stem}_deltaforge.csv"


def binder_sequence(sequence: str, receptor_sequences: list[str]) -> str:
    """Drop the longest receptor suffix. Leave the sequence unchanged otherwise."""
    match = ""
    for receptor in receptor_sequences:
        if receptor and sequence.endswith(receptor) and len(receptor) > len(match):
            match = receptor
    if not match:
        return sequence
    binder = sequence[: -len(match)]
    if not binder:
        raise ValueError("Empty binder after stripping receptor suffix")
    return binder


def load_input_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"No header in {path}")
        missing = [name for name in REQUIRED_COLUMNS if name not in reader.fieldnames]
        if missing:
            raise ValueError(f"{path} is missing columns: {', '.join(missing)}")
        return list(reader.fieldnames), list(reader)


def unique_sequences(rows: list[dict[str, str]]) -> list[tuple[str, str]]:
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


def index_runs(
    rows: list[dict[str, str]],
) -> tuple[dict[str, list[str]], dict[tuple[str, str], str]]:
    """Map each sequence to its runs, and (sequence, run) to that run's model id."""
    runs_by_sequence: dict[str, list[str]] = {}
    model_by_seq_run: dict[tuple[str, str], str] = {}
    for row in rows:
        sequence = (row.get("sequence") or "").strip()
        run = (row.get("run") or "").strip()
        model_id = (row.get("model_id") or "").strip()
        if not sequence or not run:
            continue
        runs = runs_by_sequence.setdefault(sequence, [])
        if run not in runs:
            runs.append(run)
        model_by_seq_run.setdefault((sequence, run), model_id or "binder")
    return runs_by_sequence, model_by_seq_run


def load_json_array(path: Path) -> list[dict[str, Any]]:
    if not path.is_file():
        return []
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, list):
        raise ValueError(f"Expected a JSON array in {path}")
    return [row for row in data if isinstance(row, dict)]


def load_checkpoints(outputs_dir: Path, stem: str) -> list[dict[str, Any]]:
    name = checkpoint_filename(stem)
    rows: list[dict[str, Any]] = []
    for path in sorted(outputs_dir.glob(f"*/{name}")):
        rows.extend(load_json_array(path))
    rows.extend(load_json_array(outputs_dir / unassigned_filename(stem)))
    return rows


def write_run_checkpoints(
    outputs_dir: Path,
    stem: str,
    rows: list[dict[str, Any]],
    runs_by_sequence: dict[str, list[str]],
    model_by_seq_run: dict[tuple[str, str], str],
    write_json,
) -> None:
    grouped: dict[str, list[dict[str, Any]]] = {}
    unassigned: list[dict[str, Any]] = []
    for row in rows:
        sequence = row.get("sequence")
        run_names = runs_by_sequence.get(sequence, []) if isinstance(sequence, str) else []
        if not run_names:
            unassigned.append(dict(row))
            continue
        for run in run_names:
            record = dict(row)
            record["model_id"] = model_by_seq_run.get((sequence, run), row.get("model_id"))
            grouped.setdefault(run, []).append(record)
    filename = checkpoint_filename(stem)
    for run, group in grouped.items():
        path = outputs_dir / run / filename
        refuse_overwrite(path)
        write_json(path, group)
    unassigned_path = outputs_dir / unassigned_filename(stem)
    refuse_overwrite(unassigned_path)
    if unassigned:
        write_json(unassigned_path, unassigned)
    elif unassigned_path.exists():
        unassigned_path.unlink()


def output_fieldnames(input_fields: list[str], receptor_names: list[str], metric_fields) -> list[str]:
    names = list(input_fields)
    if "binder_sequence" not in names:
        names.append("binder_sequence")
    for name in receptor_names:
        for field in metric_fields:
            column = f"{name}_{field}"
            if column not in names:
                names.append(column)
    return names


def write_summary(
    path: Path,
    fieldnames: list[str],
    input_rows: list[dict[str, str]],
    by_key: dict[tuple[str, str], dict[str, Any]],
    receptors: list[tuple[str, str]],
    extract_metrics,
) -> None:
    refuse_overwrite(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    receptor_sequences = [sequence for _name, sequence in receptors]
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in input_rows:
            sequence = (row.get("sequence") or "").strip()
            out = {name: row.get(name, "") for name in fieldnames}
            try:
                out["binder_sequence"] = binder_sequence(sequence, receptor_sequences) if sequence else ""
            except ValueError:
                out["binder_sequence"] = ""
            for name, _sequence in receptors:
                metrics = extract_metrics(by_key.get((sequence, name)))
                for field, value in metrics.items():
                    out[f"{name}_{field}"] = value
            writer.writerow(out)
    tmp.replace(path)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="Selection CSV with run, model_id, and sequence.")
    parser.add_argument(
        "--receptors",
        type=Path,
        nargs="+",
        default=[SEQUENCES / "Brd4ET.fasta", SEQUENCES / "Brd3ET.fasta"],
        help="Receptor FASTA files. The FASTA header is the column prefix.",
    )
    parser.add_argument(
        "--outputs-dir",
        type=Path,
        default=RFD3_OUT,
        help="analysis/outputs/<method> directory. Checkpoints go in <run>/ underneath it.",
    )
    parser.add_argument("--output", type=Path, help="Wide summary CSV. Defaults beside the checkpoints.")
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
    for path in (LIGAND_AI, DELTAFORGE):
        if str(path) not in sys.path:
            sys.path.insert(0, str(path))
    from score_outputs_deltaforge import load_target_sequence  # noqa: E402
    from score_top_sequences_deltaforge import (  # noqa: E402
        METRIC_FIELDS,
        extract_metrics,
        is_success,
        score_pair,
        write_json,
    )

    args = build_parser().parse_args()
    input_path = args.input.expanduser().resolve()
    outputs_dir = args.outputs_dir.expanduser().resolve()
    if not input_path.is_file():
        print(f"Input CSV not found: {input_path}", file=sys.stderr)
        return 2
    stem = input_path.stem
    output_path = (args.output or (outputs_dir / summary_filename(stem))).expanduser().resolve()
    refuse_overwrite(output_path)
    refuse_overwrite(outputs_dir / checkpoint_filename(stem))
    refuse_overwrite(outputs_dir / unassigned_filename(stem))

    receptors: list[tuple[str, str]] = []
    seen_names: set[str] = set()
    for path in args.receptors:
        fasta = path.expanduser().resolve()
        if not fasta.is_file():
            print(f"Receptor FASTA not found: {fasta}", file=sys.stderr)
            return 2
        name, sequence = load_target_sequence(fasta)
        if name in seen_names:
            print(f"Duplicate receptor name {name} from {fasta}", file=sys.stderr)
            return 2
        seen_names.add(name)
        receptors.append((name, sequence))
    receptor_sequences = [sequence for _name, sequence in receptors]

    input_fields, input_rows = load_input_rows(input_path)
    runs_by_sequence, model_by_seq_run = index_runs(input_rows)
    sequences = unique_sequences(input_rows)
    fieldnames = output_fieldnames(input_fields, [name for name, _seq in receptors], METRIC_FIELDS)

    checkpoint_rows = [] if args.force else load_checkpoints(outputs_dir, stem)
    by_key: dict[tuple[str, str], dict[str, Any]] = {}
    for row in checkpoint_rows:
        sequence = row.get("sequence")
        target_key = row.get("target_key")
        if isinstance(sequence, str) and isinstance(target_key, str):
            by_key[(sequence, target_key)] = row

    pending: list[dict[str, Any]] = []
    for sequence, model_id in sequences:
        try:
            binder = binder_sequence(sequence, receptor_sequences)
        except ValueError as exc:
            print(f"{model_id}: {exc}", file=sys.stderr)
            return 2
        for target_key, target_name, target_sequence in (
            (name, name, sequence_) for name, sequence_ in receptors
        ):
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

    cached = len(sequences) * len(receptors) - len(pending)
    receptor_label = ", ".join(f"{name} ({len(sequence)} aa)" for name, sequence in receptors)
    print(
        f"Sequences: {len(input_rows)} rows, {len(sequences)} unique; receptors {receptor_label}",
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
            for target_key, _name, _seq in (
                (name, name, sequence_) for name, sequence_ in receptors
            )
            if (sequence, target_key) in by_key
        ]
        write_run_checkpoints(
            outputs_dir,
            stem,
            ordered,
            runs_by_sequence,
            model_by_seq_run,
            write_json,
        )
        write_summary(output_path, fieldnames, input_rows, by_key, receptors, extract_metrics)

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

    print(f"\nDone: scored={scored} failed={failed} cached={cached}", flush=True)
    print(f"CSV: {output_path}", flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
