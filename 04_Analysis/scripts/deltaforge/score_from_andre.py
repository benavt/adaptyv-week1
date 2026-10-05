#!/usr/bin/env python3
"""DeltaForge-score collaborator peptides against the four ET isoforms.

Sequences are stored with the predicted scores and downloaded structures.
Collaborator run metrics are not written.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
LIGAND_AI = REPO / "design" / "ligand_ai"
SEQUENCES = REPO / "staging" / "sequences"
OUT = REPO / "design" / "ligand_ai" / "from_andre"
CHECKPOINT = OUT / "checkpoint.json"
CACHE = OUT / "_cache"

for path in (LIGAND_AI, HERE):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from score_outputs_deltaforge import load_target_sequence  # noqa: E402
from src.ligandai_local.client import LigandAIClient  # noqa: E402
from src.ligandai_local.config import has_api_key  # noqa: E402

RECEPTORS = (
    ("Brd2ET", SEQUENCES / "Brd2ET.fasta"),
    ("Brd3ET", SEQUENCES / "Brd3ET.fasta"),
    ("Brd4ET", SEQUENCES / "Brd4ET.fasta"),
    ("BrdTET", SEQUENCES / "BrdTET.fasta"),
)

BATCHES: dict[str, tuple[str, ...]] = {
    "seq_batch_1": (
        "VVAGAQLCSVDENEGRQDITEIISPHIHRNI",
        "TNNFHKQACCENVHFDCEDIEDAAFEEYGLSG",
        "VEEDDDLVDAIFIEVLVHQH",
        "SIDGATGVEDEVVSTCVDYDAIIIVEGLE",
        "SILPRELICIPEGFGEHDHDECDD",
        "FATENQHLNRSVILASTAAMSLFPYAYTDDNNTS",
        "GIGDGDCEVAADQVDHIDLDFGGPLGHCDPHTENA",
        "STEELEVYYTDVYAIAGCPLGNAFGIFTPLQPEAEFSTKVSPRIVR",
        "PILLGDNHIQECEWQTDRASDRYIMAGCQAPHSTIN",
        "GNVLESPDGGDIVCCIDIEAETTDGLQDDTDPTPVIAHN",
    ),
    "seq_batch_2": (
        "TVDYCSPAEEEKLEADTEAHL",
        "GSFVKYEDSGFRPNIVAFWLLRQIYFDKESE",
        "EKVKNSYPAMAGNVAEHGVVE",
        "TPAQISTEIDQDTKNANPQH",
        "TPAQIWTEIDQDTKNANWQH",
        "SVVSERDETSVKYPPYNVYIG",
    ),
    "seq_batch_3": (
        "SILPRELICIPEGFGEHDHDECDD",
        "SGATINEAIGCQAEIVLEGETIEDDEIE",
        "TADGRAAYMGERDDTSEAYNRKFLVVLDCKPNFE",
        "SDEVLTGDDVENYNIEAVQS",
        "SDFAGNFEHDCDFKVWEVLDVTAT",
        "PSVIDDNATKVNMWEAALVGFPINYAPEDDI",
        "GCADDADAHGIAELDDKSDDGINV",
        "PSPPTHDIPDLSRILEFKDDTEDADVK",
        "GLEVDAGAITALICFGPFLVEQDEPDDH",
        "GSTAVRVHKGLNRCEGFDVLETGEE",
    ),
    "seq_batch_4": (
        "AFALPCYHTFLSAFSAYKLDELDVYVGYTT",
        "SLQQCIPLYCVLKLQGEYIRAEDIGELLVI",
        "TVDYCSPAEEEKLEADTEAHL",
        "TPAQISTEIDQDTKNANPQH",
        "TPAQIWTEIDQDTKNANWQH",
        "EKVKNSYPAMAGNVAEHGVVE",
        "PLERSELDAGNTHSSVDIEEGLEGFVES",
        "MVVLKSCATVPLPAPNMCHYTPHLSPNSCKPEATLDAET",
        "GPGLEIGEAFINGKSSTLIEEADEPVYNE",
        "TELDITVLIHIENCIPFELDEGVSDV",
    ),
    "seq_batch_5": (
        "ADEPTSQSHDPCNLYEDRIDGGN",
        "TDKLRGGWNTQSTQKFLMECFPEIDEK",
        "TPAQIWTEIDQDTKNANWQH",
        "PSVIDDNATKVNMWEAALVGFPINYAPEDDI",
        "SDAADETALCSAVRASAWFYEYVVLRTAQKPGAC",
        "ATISTEYDQDTKNANPAFEKH",
    ),
}

EXPECTED_LENGTHS = {
    "seq_batch_1": (31, 32, 20, 29, 24, 34, 35, 46, 36, 39),
    "seq_batch_2": (21, 31, 21, 20, 20, 21),
    "seq_batch_3": (24, 28, 34, 20, 24, 31, 24, 27, 28, 25),
    "seq_batch_4": (30, 30, 21, 20, 20, 21, 28, 39, 29, 26),
    "seq_batch_5": (23, 27, 20, 31, 34, 21),
}

SCORE_KEYS = ("kd_nm", "delta_g", "iptm", "ptm", "ipsae", "plddt_mean", "fold_job_id", "error")


def cache_stem(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("utf-8")).hexdigest()[:16]


def cache_paths(sequence: str, isoform: str) -> tuple[Path, Path]:
    directory = CACHE / isoform
    stem = cache_stem(sequence)
    return directory / f"{stem}.pdb", directory / f"{stem}.cif"


def structures_ready(sequence: str, isoform: str) -> bool:
    pdb, cif = cache_paths(sequence, isoform)
    return pdb.is_file() and cif.is_file() and pdb.stat().st_size > 0 and cif.stat().st_size > 0


def as_float(value: Any) -> float | None:
    if isinstance(value, bool) or value is None or value == "":
        return None
    if isinstance(value, list) and value:
        numbers = [as_float(item) for item in value]
        numbers = [item for item in numbers if item is not None]
        if not numbers:
            return None
        return sum(numbers) / len(numbers)
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _sources(payload: dict[str, Any]) -> list[dict[str, Any]]:
    ordered: list[dict[str, Any]] = []
    for key in ("scoring", "deltaforge", "result", "score"):
        inner = payload.get(key)
        if isinstance(inner, dict):
            ordered.append(inner)
    ordered.append(payload)
    return ordered


def _pick(sources: list[dict[str, Any]], keys: tuple[str, ...]) -> Any:
    for source in sources:
        for key in keys:
            if key in source and source[key] is not None:
                return source[key]
    return None


def _fold_ids(payload: Any, found: list[str]) -> None:
    if isinstance(payload, dict):
        for key, value in payload.items():
            if key in {"foldJobId", "fold_job_id", "jobId", "job_id", "id"} and isinstance(value, str):
                if value.startswith("fold_") and value not in found:
                    found.append(value)
            else:
                _fold_ids(value, found)
    elif isinstance(payload, list):
        for item in payload:
            _fold_ids(item, found)


def metrics_from_payload(payload: dict[str, Any]) -> dict[str, Any]:
    sources = _sources(payload)
    kd = as_float(_pick(sources, ("kd_nm", "kdNm", "predicted_kd", "predictedKd", "kd")))
    fold_ids: list[str] = []
    _fold_ids(payload, fold_ids)
    return {
        "kd_nm": kd,
        "delta_g": as_float(_pick(sources, ("delta_g", "deltaG", "dg"))),
        "iptm": as_float(_pick(sources, ("iptm", "iPTM", "fold_iptm", "foldIptm"))),
        "ptm": as_float(_pick(sources, ("ptm", "pTM", "fold_ptm", "foldPtm"))),
        "ipsae": as_float(_pick(sources, ("ipsae", "iPSAE", "peptide_ipsae", "peptideIpsae", "fold_ipsae"))),
        "plddt_mean": as_float(
            _pick(sources, ("plddt_mean", "plddtMean", "mean_plddt", "fold_plddt_mean", "plddt"))
        ),
        "fold_job_id": fold_ids[0] if fold_ids else None,
        "fold_job_ids": fold_ids,
    }


def write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(data, indent=2, default=str) + "\n", encoding="utf-8")
    tmp.replace(path)


def load_checkpoint() -> dict[tuple[str, str], dict[str, Any]]:
    if not CHECKPOINT.exists():
        return {}
    data = json.loads(CHECKPOINT.read_text(encoding="utf-8"))
    if not isinstance(data, list):
        raise ValueError(f"Expected a JSON array in {CHECKPOINT}")
    rows: dict[tuple[str, str], dict[str, Any]] = {}
    for row in data:
        if not isinstance(row, dict):
            continue
        sequence = row.get("sequence")
        isoform = row.get("isoform")
        if isinstance(sequence, str) and isinstance(isoform, str):
            rows[(sequence, isoform)] = row
    return rows


def checkpoint_record(sequence: str, isoform: str, fields: dict[str, Any]) -> dict[str, Any]:
    record = {"sequence": sequence, "isoform": isoform}
    for key in SCORE_KEYS:
        record[key] = fields.get(key)
    return record


def is_complete(row: dict[str, Any] | None) -> bool:
    if not row or row.get("kd_nm") is None or row.get("error"):
        return False
    sequence = row.get("sequence")
    isoform = row.get("isoform")
    if not isinstance(sequence, str) or not isinstance(isoform, str):
        return False
    return structures_ready(sequence, isoform)


def validate_batches() -> list[str]:
    unique: list[str] = []
    seen: set[str] = set()
    for name, sequences in BATCHES.items():
        expected = EXPECTED_LENGTHS[name]
        if len(sequences) != len(expected):
            raise ValueError(f"{name} has {len(sequences)} sequences, expected {len(expected)}")
        for index, sequence in enumerate(sequences):
            if len(sequence) != expected[index]:
                raise ValueError(
                    f"{name} seq_{index + 1:02d} length {len(sequence)} != {expected[index]}"
                )
            if sequence not in seen:
                seen.add(sequence)
                unique.append(sequence)
    if len(unique) != 35:
        raise ValueError(f"Expected 35 unique sequences, found {len(unique)}")
    return unique


def isoform_entry(sequence: str, isoform: str, row: dict[str, Any] | None, seq_id: str) -> dict[str, Any]:
    pdb_rel = f"structures/{isoform}/{seq_id}.pdb"
    cif_rel = f"structures/{isoform}/{seq_id}.cif"
    entry = {
        "kd_nm": None,
        "delta_g": None,
        "iptm": None,
        "ptm": None,
        "ipsae": None,
        "plddt_mean": None,
        "fold_job_id": None,
        "pdb": pdb_rel,
        "cif": cif_rel,
        "error": None,
    }
    if row:
        for key in ("kd_nm", "delta_g", "iptm", "ptm", "ipsae", "plddt_mean", "fold_job_id", "error"):
            entry[key] = row.get(key)
    else:
        entry["error"] = "not scored"
    return entry


def materialize(rows: dict[tuple[str, str], dict[str, Any]]) -> None:
    for batch_name, sequences in BATCHES.items():
        batch_dir = OUT / batch_name
        payload_sequences: list[dict[str, Any]] = []
        for index, sequence in enumerate(sequences, start=1):
            seq_id = f"seq_{index:02d}"
            isoforms: dict[str, Any] = {}
            for isoform, _path in RECEPTORS:
                row = rows.get((sequence, isoform))
                entry = isoform_entry(sequence, isoform, row, seq_id)
                isoforms[isoform] = entry
                if row and structures_ready(sequence, isoform) and not row.get("error"):
                    pdb_src, cif_src = cache_paths(sequence, isoform)
                    for src, relative in (
                        (pdb_src, entry["pdb"]),
                        (cif_src, entry["cif"]),
                    ):
                        dest = batch_dir / relative
                        dest.parent.mkdir(parents=True, exist_ok=True)
                        if not dest.is_file() or dest.stat().st_size != src.stat().st_size:
                            shutil.copyfile(src, dest)
            payload_sequences.append(
                {
                    "id": seq_id,
                    "sequence": sequence,
                    "length": len(sequence),
                    "isoforms": isoforms,
                }
            )
        write_json(batch_dir / "batch.json", {"batch": batch_name, "sequences": payload_sequences})


def save_state(rows: dict[tuple[str, str], dict[str, Any]]) -> None:
    ordered = [
        rows[(sequence, isoform)]
        for sequence in validate_batches()
        for isoform, _path in RECEPTORS
        if (sequence, isoform) in rows
    ]
    write_json(CHECKPOINT, ordered)
    materialize(rows)


def fold_confidence(client: LigandAIClient, fold_job_id: str) -> dict[str, float | None]:
    sdk = client._build_sdk_client()
    payload = sdk.transport.request("GET", f"/api/folding/jobs/{fold_job_id}") or {}
    result = payload.get("result") if isinstance(payload, dict) else None
    if not isinstance(result, dict):
        result = {}
    plddt = result.get("mean_plddt")
    if plddt is None:
        plddt = result.get("plddt")
    return {
        "iptm": as_float(result.get("iptm")),
        "ptm": as_float(result.get("ptm")),
        "ipsae": as_float(result.get("peptide_ipsae") if result.get("peptide_ipsae") is not None else result.get("ipsae")),
        "plddt_mean": as_float(plddt),
    }


def enrich_confidence(client: LigandAIClient, fields: dict[str, Any]) -> None:
    fold_job_id = fields.get("fold_job_id")
    if not isinstance(fold_job_id, str) or not fold_job_id:
        return
    missing = [key for key in ("iptm", "ptm", "plddt_mean") if fields.get(key) is None]
    if not missing and fields.get("ipsae") is not None:
        return
    try:
        confidence = fold_confidence(client, fold_job_id)
    except Exception as exc:  # noqa: BLE001
        print(f"  fold metrics {fold_job_id}: {type(exc).__name__}: {exc}", flush=True)
        return
    for key in missing:
        if confidence.get(key) is not None:
            fields[key] = confidence[key]
    if fields.get("ipsae") is None and confidence.get("ipsae") is not None:
        fields["ipsae"] = confidence["ipsae"]


def download_structures(client: LigandAIClient, sequence: str, isoform: str, fold_ids: list[str]) -> tuple[str | None, str | None]:
    pdb_path, cif_path = cache_paths(sequence, isoform)
    pdb_path.parent.mkdir(parents=True, exist_ok=True)
    errors: list[str] = []
    for fold_id in fold_ids:
        try:
            structures = client.download_folding_job_structures(fold_id)
        except Exception as exc:  # noqa: BLE001
            errors.append(f"{fold_id}: {type(exc).__name__}: {exc}")
            continue
        pdb = structures.get("pdb") or ""
        cif = structures.get("cif") or ""
        if "ATOM" not in pdb or not cif.strip():
            errors.append(f"{fold_id}: structure download was empty")
            continue
        pdb_path.write_text(pdb if pdb.endswith("\n") else pdb + "\n", encoding="utf-8")
        cif_path.write_text(cif if cif.endswith("\n") else cif + "\n", encoding="utf-8")
        return fold_id, None
    return None, ("; ".join(errors) if errors else "no fold job id")


def fold_peptide(
    client: LigandAIClient,
    sequence: str,
    isoform: str,
    target_sequence: str,
    *,
    timeout: float,
) -> dict[str, Any]:
    """Fold the peptide against one isoform and read the auto-score.

    Chain B is marked as a peptide so the platform does not request an MSA
    for it. Chain A keeps MSA on, which hits the cached receptor alignment.
    """
    payload = {
        "model": "boltz2",
        "autoScore": True,
        "templateMode": False,
        "gpuCount": 1,
        "diffusionSamples": 1,
        "targetGeneName": isoform,
        "entities": [
            {
                "type": "protein",
                "chainId": "A",
                "sequence": target_sequence,
                "name": isoform,
                "geneName": isoform,
                "useMsa": True,
                "use_msa": True,
            },
            {
                "type": "protein",
                "chainId": "B",
                "sequence": sequence,
                "name": "peptide",
                "isPeptide": True,
                "useMsa": False,
                "use_msa": False,
            },
        ],
    }
    started = client._post_json(
        "/api/folding/predict",
        payload,
        auth_mode="bearer_and_api_key",
        max_retries=3,
        retry_base_seconds=10.0,
    )
    job_id = started.get("jobId") or started.get("id") or started.get("job_id")
    if not isinstance(job_id, str) or not job_id:
        raise RuntimeError("folding predict did not return a job id")
    deadline = time.time() + timeout
    while True:
        info = client.get_folding_job(job_id)
        status = str(info.get("status") or "").lower()
        if status in {"completed", "complete"}:
            info["id"] = info.get("id") or job_id
            return info
        if status in {"failed", "error", "cancelled"}:
            detail = info.get("error") or info.get("message") or status
            raise RuntimeError(f"DeltaForge job {job_id} ended with status {status} ({detail})")
        if time.time() >= deadline:
            raise TimeoutError(f"DeltaForge job {job_id} did not complete within {timeout:.0f}s (last status: {status})")
        time.sleep(15)


def score_pair(
    sequence: str,
    isoform: str,
    target_sequence: str,
    *,
    max_attempts: int,
    backoff_s: float,
    existing: dict[str, Any] | None,
) -> dict[str, Any]:
    label = f"{sequence[:8]}… vs {isoform}"
    client = LigandAIClient()
    fields: dict[str, Any] = {key: None for key in SCORE_KEYS}
    if existing:
        fields.update({key: existing.get(key) for key in SCORE_KEYS})

    needs_score = fields.get("kd_nm") is None
    if needs_score:
        print(f"  scoring {label}", flush=True)
        payload = None
        last_error = ""
        for attempt in range(1, max_attempts + 1):
            try:
                payload = fold_peptide(
                    client,
                    sequence,
                    isoform,
                    target_sequence,
                    timeout=1200,
                )
                break
            except Exception as exc:  # noqa: BLE001
                last_error = f"{type(exc).__name__}: {exc}"
                text = last_error.lower()
                retryable = any(
                    needle in text
                    for needle in (
                        "msa failure",
                        "timeout",
                        "http 502",
                        "http 503",
                        "http 504",
                        "http 522",
                        "temporarily",
                        "connection",
                        "gpu queue",
                    )
                )
                if not retryable or attempt >= max_attempts:
                    fields["error"] = last_error
                    print(f"  ERROR {label}: {fields['error']}", flush=True)
                    return checkpoint_record(sequence, isoform, fields)
                sleep_for = min(backoff_s * (2 ** min(attempt - 1, 4)), 120.0)
                print(
                    f"    retryable failure (attempt {attempt}/{max_attempts}): {last_error}; sleeping {sleep_for:.0f}s",
                    flush=True,
                )
                time.sleep(sleep_for)
        if payload is None:
            fields["error"] = last_error or "fold did not return a payload"
            return checkpoint_record(sequence, isoform, fields)
        if not isinstance(payload, dict):
            fields["error"] = "score payload was not an object"
            return checkpoint_record(sequence, isoform, fields)
        parsed = metrics_from_payload(payload)
        if parsed["kd_nm"] is None:
            nested = [key for key in ("scoring", "deltaforge", "result", "score") if isinstance(payload.get(key), dict)]
            print(
                f"  no kd in payload keys={sorted(payload)} nested={nested}",
                flush=True,
            )
        fields.update(
            {
                "kd_nm": parsed["kd_nm"],
                "delta_g": parsed["delta_g"],
                "iptm": parsed["iptm"],
                "ptm": parsed["ptm"],
                "ipsae": parsed["ipsae"],
                "plddt_mean": parsed["plddt_mean"],
                "fold_job_id": parsed["fold_job_id"],
                "error": None if parsed["kd_nm"] is not None else "fold completed without kd_nm",
            }
        )
        fold_ids = list(parsed["fold_job_ids"])
    else:
        print(f"  downloading {label}", flush=True)
        fold_ids = []
        if isinstance(fields.get("fold_job_id"), str):
            fold_ids.append(fields["fold_job_id"])

    if fields.get("kd_nm") is None:
        return checkpoint_record(sequence, isoform, fields)

    if not structures_ready(sequence, isoform):
        saved_id, download_error = download_structures(client, sequence, isoform, fold_ids)
        if download_error is not None or not structures_ready(sequence, isoform):
            fields["error"] = download_error or "structure download failed"
            print(f"  ERROR {label}: {fields['error']}", flush=True)
            return checkpoint_record(sequence, isoform, fields)
        if saved_id:
            fields["fold_job_id"] = saved_id

    fields["error"] = None
    enrich_confidence(client, fields)
    print(f"  ok {label} kd_nm={fields['kd_nm']}", flush=True)
    return checkpoint_record(sequence, isoform, fields)


def main() -> int:
    parser = argparse.ArgumentParser(description="DeltaForge-score Andre sequence batches against four ET isoforms")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--max-attempts", type=int, default=8)
    parser.add_argument("--backoff-s", type=float, default=30.0)
    args = parser.parse_args()

    if not has_api_key():
        print("set LIGANDAI_API_KEY or save a key with scripts/set_ligandai_api_key.sh", file=sys.stderr)
        return 2

    unique = validate_batches()
    receptors: list[tuple[str, str]] = []
    for name, path in RECEPTORS:
        _header, sequence = load_target_sequence(path)
        receptors.append((name, sequence))

    rows = load_checkpoint()
    backfill = LigandAIClient()
    backfilled = 0
    for key, row in list(rows.items()):
        if row.get("kd_nm") is None or not row.get("fold_job_id"):
            continue
        if row.get("iptm") is not None and row.get("ptm") is not None and row.get("plddt_mean") is not None:
            continue
        enrich_confidence(backfill, row)
        rows[key] = checkpoint_record(row["sequence"], row["isoform"], row)
        backfilled += 1
    if backfilled:
        print(f"Backfilled fold metrics for {backfilled} pairs", flush=True)
        save_state(rows)

    pending: list[tuple[str, str, str]] = []
    for sequence in unique:
        for isoform, target_sequence in receptors:
            if is_complete(rows.get((sequence, isoform))):
                continue
            pending.append((sequence, isoform, target_sequence))

    print(
        f"Unique sequences: {len(unique)}; jobs {len(pending)} to score, "
        f"{len(unique) * len(receptors) - len(pending)} cached; workers={args.workers}",
        flush=True,
    )
    save_state(rows)

    lock = threading.Lock()
    scored = 0
    failed = 0
    if pending:
        with ThreadPoolExecutor(max_workers=max(1, args.workers)) as pool:
            futures = [
                pool.submit(
                    score_pair,
                    sequence,
                    isoform,
                    target_sequence,
                    max_attempts=args.max_attempts,
                    backoff_s=args.backoff_s,
                    existing=rows.get((sequence, isoform)),
                )
                for sequence, isoform, target_sequence in pending
            ]
            for future in as_completed(futures):
                record = future.result()
                with lock:
                    rows[(record["sequence"], record["isoform"])] = record
                    if is_complete(record):
                        scored += 1
                    else:
                        failed += 1
                    save_state(rows)
                    print(
                        f"  progress {scored + failed}/{len(pending)} scored={scored} failed={failed}",
                        flush=True,
                    )

    save_state(rows)
    print(f"Checkpoint: {CHECKPOINT}", flush=True)
    print(f"Done: scored={scored} failed={failed}", flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
