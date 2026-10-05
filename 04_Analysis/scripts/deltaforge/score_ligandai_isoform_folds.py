#!/usr/bin/env python3
"""Fold each kept LigandAI design against four ET sequences and record Kd.

Chain A is the isoform sequence. Chain B is the designed binder and is flagged
as a peptide so the platform does not request an MSA for it. This does not use
the binder-scoring fold-and-score endpoint.
"""

from __future__ import annotations

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
OUT = REPO / "analysis" / "outputs" / "ligand_ai" / "isoform_fingerprint"
JOBS = OUT / "top_10pct_predicted_kd_jobs.csv"
SUMMARY = OUT / "top_10pct_cross_isoform_kd.csv"
SEQUENCES = REPO / "staging" / "sequences"
RECEPTOR_FASTAS = (
    SEQUENCES / "Brd2ET.fasta",
    SEQUENCES / "Brd3ET.fasta",
    SEQUENCES / "Brd4ET.fasta",
    SEQUENCES / "BrdTET.fasta",
)
SMOKE = {
    "sequence": "LDENQSLGKVPLSNLRVRDGIYNKGF",
    "model_id": "BRD3ET_rank1",
    "run": "BRD3ET",
    "target_key": "BRD3ET",
    "job_id": "fold_1790181009573_rdv4suh13",
    "kd_nm": 5.8975635,
    "delta_g": -11.226835,
    "iptm": 0.9298084,
    "ptm": 0.940244,
    "ipsae": 0.6728436,
    "peptide_ipsae": 0.6728436,
    "plddt": 85.07092,
    "error": None,
}
METRIC_FIELDS = (
    "kd_nm",
    "delta_g",
    "iptm",
    "ptm",
    "ipsae",
    "peptide_ipsae",
    "plddt",
    "job_id",
    "error",
)

if str(LIGAND_AI) not in sys.path:
    sys.path.insert(0, str(LIGAND_AI))
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from score_outputs_deltaforge import load_target_sequence  # noqa: E402
from src.ligandai_local.config import sync_api_key_to_env  # noqa: E402


def load_jobs() -> tuple[list[str], list[dict[str, str]]]:
    with JOBS.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"No header in {JOBS}")
        return list(reader.fieldnames), list(reader)


def checkpoint_path(run: str) -> Path:
    return OUT / run / "fold_kd_checkpoint.json"


def load_checkpoints() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for path in sorted(OUT.glob("*/fold_kd_checkpoint.json")):
        data = json.loads(path.read_text(encoding="utf-8"))
        if isinstance(data, list):
            rows.extend(row for row in data if isinstance(row, dict))
    return rows


def is_success(row: dict[str, Any] | None) -> bool:
    return bool(row) and row.get("kd_nm") is not None and not row.get("error")


def as_float(value: Any) -> float | None:
    if isinstance(value, bool) or value is None or value == "":
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def metrics_from_payload(payload: dict[str, Any], job_id: str) -> dict[str, Any]:
    result = payload.get("result") if isinstance(payload.get("result"), dict) else {}
    kd = as_float(result.get("predicted_kd") if result.get("predicted_kd") is not None else result.get("predictedKd"))
    delta_g = as_float(result.get("delta_g") if result.get("delta_g") is not None else result.get("deltaG"))
    return {
        "kd_nm": kd,
        "delta_g": delta_g,
        "iptm": as_float(result.get("iptm")),
        "ptm": as_float(result.get("ptm")),
        "ipsae": as_float(result.get("ipsae") if result.get("ipsae") is not None else result.get("peptide_ipsae")),
        "peptide_ipsae": as_float(result.get("peptide_ipsae") if result.get("peptide_ipsae") is not None else result.get("peptideIpsae")),
        "plddt": as_float(result.get("plddt") if result.get("plddt") is not None else result.get("mean_plddt")),
        "job_id": job_id,
        "error": None if kd is not None else "fold completed without predicted_kd",
    }


def release_local_slots(client: Any, job_id: str | None = None) -> None:
    """The SDK counts unfinished folds against a local cap and never clears them."""
    submitted = getattr(client, "submitted_set", None)
    key_hash = getattr(client, "api_key_hash", None)
    if submitted is None or not key_hash:
        return
    if job_id:
        rows = [
            row for row in submitted.list_in_flight(key_hash)
            if str(row.get("job_id") or "") == job_id
        ]
    else:
        rows = list(submitted.list_in_flight(key_hash))
    for row in rows:
        submitted.mark_completed(str(row["submission_hash"]), key_hash)


def score_pair(job: dict[str, Any]) -> dict[str, Any]:
    from ligandai import LigandAI

    label = "%s vs %s" % (job["model_id"], job["target_key"])
    print("  scoring %s" % label, flush=True)
    client = LigandAI()
    try:
        handle = client.fold(
            job["target_sequence"],
            job["binder_sequence"],
            auto_score=True,
            target_gene=job["target_key"],
        )
        handle.wait(timeout=1200, poll_interval=15)
        release_local_slots(client, handle.id)
        payload = client.transport.request("GET", "/api/folding/jobs/%s" % handle.id) or {}
        metrics = metrics_from_payload(payload if isinstance(payload, dict) else {}, handle.id)
    except Exception as exc:  # noqa: BLE001
        print("  ERROR: %s: %s" % (label, exc), flush=True)
        metrics = {field: None for field in METRIC_FIELDS}
        metrics["error"] = "%s: %s" % (type(exc).__name__, exc)
    else:
        if metrics.get("error"):
            print("  ERROR: %s: %s" % (label, metrics["error"]), flush=True)
        else:
            print("  ok: %s kd_nm=%s" % (label, metrics["kd_nm"]), flush=True)
    return {
        "sequence": job["binder_sequence"],
        "model_id": job["model_id"],
        "run": job["run"],
        "target_key": job["target_key"],
        **metrics,
    }


def write_checkpoints(rows: list[dict[str, Any]]) -> None:
    grouped: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        run = row.get("run")
        if isinstance(run, str) and run:
            grouped.setdefault(run, []).append(row)
    for run, group in grouped.items():
        path = checkpoint_path(run)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(group, indent=2, default=str) + "\n", encoding="utf-8")


def fmt(value: Any) -> str:
    if value is None or value == "":
        return ""
    if isinstance(value, float):
        return format(value, ".8g")
    return str(value)


def write_summary(
    fieldnames: list[str],
    input_rows: list[dict[str, str]],
    by_key: dict[tuple[str, str], dict[str, Any]],
    receptor_names: list[str],
) -> None:
    names = list(fieldnames)
    for name in receptor_names:
        for field in METRIC_FIELDS:
            column = "%s_%s" % (name, field)
            if column not in names:
                names.append(column)
    tmp = SUMMARY.with_suffix(".csv.tmp")
    with tmp.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=names, extrasaction="ignore")
        writer.writeheader()
        for row in input_rows:
            sequence = (row.get("sequence") or "").strip()
            out = {name: row.get(name, "") for name in names}
            for name in receptor_names:
                metrics = by_key.get((sequence, name)) or {}
                for field in METRIC_FIELDS:
                    out["%s_%s" % (name, field)] = fmt(metrics.get(field))
            writer.writerow(out)
    tmp.replace(SUMMARY)


def main() -> int:
    if not sync_api_key_to_env():
        print("set LIGANDAI_API_KEY first", file=sys.stderr)
        return 2
    from ligandai import LigandAI

    release_local_slots(LigandAI())
    receptors: list[tuple[str, str]] = []
    for path in RECEPTOR_FASTAS:
        receptors.append(load_target_sequence(path))
    receptor_names = [name for name, _sequence in receptors]
    input_fields, input_rows = load_jobs()

    by_key: dict[tuple[str, str], dict[str, Any]] = {}
    for row in load_checkpoints():
        sequence = row.get("sequence")
        target_key = row.get("target_key")
        if isinstance(sequence, str) and isinstance(target_key, str):
            by_key[(sequence, target_key)] = row
    smoke_key = (SMOKE["sequence"], SMOKE["target_key"])
    if not is_success(by_key.get(smoke_key)):
        by_key[smoke_key] = dict(SMOKE)

    pending: list[dict[str, Any]] = []
    for row in input_rows:
        sequence = (row.get("sequence") or "").strip()
        run = (row.get("run") or row.get("gene") or "").strip()
        model_id = (row.get("model_id") or "").strip()
        for target_key, target_sequence in receptors:
            if is_success(by_key.get((sequence, target_key))):
                continue
            pending.append(
                {
                    "binder_sequence": sequence,
                    "model_id": model_id,
                    "run": run,
                    "target_key": target_key,
                    "target_sequence": target_sequence,
                }
            )

    print(
        "Jobs: %d to score, %d cached; receptors %s"
        % (len(pending), len(input_rows) * len(receptors) - len(pending), ", ".join(receptor_names)),
        flush=True,
    )
    lock = threading.Lock()
    scored = 0
    failed = 0

    def persist() -> None:
        write_checkpoints(list(by_key.values()))
        write_summary(input_fields, input_rows, by_key, receptor_names)

    persist()
    if pending:
        with ThreadPoolExecutor(max_workers=4) as pool:
            futures = [pool.submit(score_pair, job) for job in pending]
            for future in as_completed(futures):
                row = future.result()
                with lock:
                    by_key[(row["sequence"], row["target_key"])] = row
                    if is_success(row):
                        scored += 1
                    else:
                        failed += 1
                    persist()
                    print(
                        "  progress %d/%d scored=%d failed=%d" % (scored + failed, len(pending), scored, failed),
                        flush=True,
                    )
    print("CSV: %s" % SUMMARY, flush=True)
    print("Done: scored=%d failed=%d" % (scored, failed), flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
