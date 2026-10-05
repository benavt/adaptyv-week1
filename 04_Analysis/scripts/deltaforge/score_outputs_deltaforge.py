#!/usr/bin/env python3
"""Score outputs/*/sequences.fa binders with DeltaForge against Brd4ET.fasta.

Each FASTA entry is treated as binder+target fusion ending with Brd4ET.
Results are written as deltaforge_scores.json beside each sequences.fa
(resume-safe: successful entries are skipped on re-run).
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any

LIGAND_AI = Path(__file__).resolve().parents[3] / "design" / "ligand_ai"
SEQUENCES = Path(__file__).resolve().parents[3] / "staging" / "sequences"
if str(LIGAND_AI) not in sys.path:
    sys.path.insert(0, str(LIGAND_AI))

from src.ligandai_local.client import LigandAIClient  # noqa: E402


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks)))
            header = line[1:].strip()
            chunks = []
        else:
            chunks.append(line.replace(" ", "").upper())
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def load_target_sequence(path: Path) -> tuple[str, str]:
    records = read_fasta(path)
    if len(records) != 1:
        raise ValueError(f"Expected exactly one sequence in {path}, found {len(records)}")
    name, seq = records[0]
    if not seq:
        raise ValueError(f"Empty target sequence in {path}")
    return name or "Brd4ET", seq


def binder_id(header: str) -> str:
    return header.split(",", 1)[0].strip() or header.strip()


def load_existing(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, list):
        raise ValueError(f"Expected a JSON array in {path}")
    return [row for row in data if isinstance(row, dict)]


def index_by_header(rows: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    out: dict[str, dict[str, Any]] = {}
    for row in rows:
        header = row.get("header")
        if isinstance(header, str) and header:
            out[header] = row
    return out


def write_scores(path: Path, rows: list[dict[str, Any]]) -> None:
    path.write_text(json.dumps(rows, indent=2, default=str) + "\n", encoding="utf-8")


def _is_retryable(exc: BaseException) -> bool:
    text = f"{type(exc).__name__}: {exc}".lower()
    needles = (
        "gpu queue unavailable",
        "redis",
        "ended with status failed",
        "http 522",
        "http 502",
        "http 503",
        "http 504",
        "timeout",
        "temporarily",
        "connection reset",
        "connection aborted",
        "connection error",
        "nodename nor servname",
        "name or service not known",
        "temporary failure in name resolution",
        "network is unreachable",
    )
    return any(n in text for n in needles)


def _enrich_job_error(client: LigandAIClient, exc: BaseException) -> str:
    """Attach folding-job error detail when the SDK only reports status=failed."""
    base = f"{type(exc).__name__}: {exc}"
    text = str(exc)
    marker = "Job "
    if marker not in text or " ended with status " not in text:
        return base
    try:
        job_id = text.split(marker, 1)[1].split(" ", 1)[0].strip()
    except Exception:  # noqa: BLE001
        return base
    if not job_id.startswith("fold_"):
        return base
    try:
        sdk = client._build_sdk_client()
        payload = sdk.transport.request("GET", f"/api/folding/jobs/{job_id}") or {}
        detail = payload.get("error") or payload.get("message")
        if detail:
            return f"{base} ({detail})"
    except Exception:  # noqa: BLE001
        pass
    return base


def score_with_retries(
    client: LigandAIClient,
    *,
    binder_sequence: str,
    target_sequence: str,
    binder_name: str,
    target_name: str,
    scorer: str,
    max_attempts: int,
    backoff_s: float,
) -> Any:
    """Retry scoring. GPU-queue / Redis failures retry indefinitely with capped backoff."""
    attempt = 0
    while True:
        attempt += 1
        try:
            return client.score_deltaforge(
                binder_sequence=binder_sequence,
                target_sequence=target_sequence,
                binder_name=binder_name,
                target_name=target_name,
                scorer=scorer,
            )
        except Exception as exc:  # noqa: BLE001
            detail = _enrich_job_error(client, exc)
            detail_l = detail.lower()
            gpu_down = "gpu queue unavailable" in detail_l or "redis" in detail_l
            if not _is_retryable(exc):
                try:
                    exc.args = (detail,)
                except Exception:  # noqa: BLE001
                    pass
                raise
            if not gpu_down and attempt >= max_attempts:
                try:
                    exc.args = (detail,)
                except Exception:  # noqa: BLE001
                    pass
                raise
            sleep_for = min(backoff_s * (2 ** min(attempt - 1, 4)), 120.0)
            label = "GPU queue" if gpu_down else f"attempt {attempt}/{max_attempts}"
            print(
                f"    retryable failure ({label}): {detail}; sleeping {sleep_for:.0f}s",
                flush=True,
            )
            time.sleep(sleep_for)


def score_directory(
    client: LigandAIClient,
    sequences_fa: Path,
    target_name: str,
    target_sequence: str,
    *,
    fusion_suffix: str,
    scores_name: str,
    scorer: str,
    force: bool,
    max_attempts: int,
    backoff_s: float,
) -> tuple[int, int, int]:
    out_path = sequences_fa.with_name(scores_name)
    existing = load_existing(out_path)
    by_header = index_by_header(existing)
    ordered: list[dict[str, Any]] = []

    scored = 0
    skipped = 0
    failed = 0

    records = read_fasta(sequences_fa)
    print(f"\n=== {sequences_fa.parent.name}: {len(records)} sequences -> {out_path}", flush=True)

    for i, (header, fusion) in enumerate(records, start=1):
        prior = by_header.get(header)
        if (
            not force
            and prior is not None
            and prior.get("result") is not None
            and not prior.get("error")
        ):
            ordered.append(prior)
            skipped += 1
            print(f"  [{i}/{len(records)}] skip (cached): {binder_id(header)}", flush=True)
            continue

        if not fusion.endswith(fusion_suffix):
            row = {
                "header": header,
                "binder_sequence": None,
                "target_name": target_name,
                "error": (
                    "Fusion sequence does not end with expected fusion-suffix "
                    "sequence; cannot extract binder."
                ),
            }
            ordered.append(row)
            failed += 1
            write_scores(out_path, ordered)
            print(f"  [{i}/{len(records)}] FAIL suffix: {binder_id(header)}", flush=True)
            continue

        binder = fusion[: -len(fusion_suffix)]
        if not binder:
            row = {
                "header": header,
                "binder_sequence": "",
                "target_name": target_name,
                "error": "Empty binder after stripping fusion suffix.",
            }
            ordered.append(row)
            failed += 1
            write_scores(out_path, ordered)
            print(f"  [{i}/{len(records)}] FAIL empty binder: {binder_id(header)}", flush=True)
            continue

        name = binder_id(header)
        print(
            f"  [{i}/{len(records)}] scoring {name} (binder_len={len(binder)}) "
            f"vs {target_name} ...",
            flush=True,
        )
        try:
            result = score_with_retries(
                client,
                binder_sequence=binder,
                target_sequence=target_sequence,
                binder_name=name,
                target_name=target_name,
                scorer=scorer,
                max_attempts=max_attempts,
                backoff_s=backoff_s,
            )
            row = {
                "header": header,
                "binder_sequence": binder,
                "target_name": target_name,
                "result": result,
            }
            scored += 1
            print(f"  [{i}/{len(records)}] ok: {name}", flush=True)
        except Exception as exc:  # noqa: BLE001
            detail = _enrich_job_error(client, exc)
            row = {
                "header": header,
                "binder_sequence": binder,
                "target_name": target_name,
                "error": detail,
            }
            failed += 1
            print(f"  [{i}/{len(records)}] ERROR: {name}: {detail}", flush=True)

        ordered.append(row)
        write_scores(out_path, ordered)

    write_scores(out_path, ordered)
    print(
        f"  done {sequences_fa.parent.name}: scored={scored} skipped={skipped} failed={failed}",
        flush=True,
    )
    return scored, skipped, failed


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="DeltaForge-score binders in outputs/*/sequences.fa against Brd4ET.fasta"
    )
    parser.add_argument(
        "--outputs-dir",
        type=Path,
        default=LIGAND_AI / "outputs",
        help="Directory containing */sequences.fa (default: repo outputs/)",
    )
    parser.add_argument(
        "--target-fasta",
        type=Path,
        default=SEQUENCES / "Brd4ET.fasta",
        help="Receptor/target FASTA used for DeltaForge scoring (default: Brd4ET.fasta)",
    )
    parser.add_argument(
        "--fusion-suffix-fasta",
        type=Path,
        default=None,
        help=(
            "FASTA sequence stripped from the C-terminus of each fusion to recover "
            "the binder (default: same as --target-fasta)."
        ),
    )
    parser.add_argument(
        "--scores-name",
        default="deltaforge_scores.json",
        help="Output JSON filename beside each sequences.fa (default: deltaforge_scores.json)",
    )
    parser.add_argument("--scorer", default="auto")
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-score even when a successful result is already cached.",
    )
    parser.add_argument(
        "--only",
        action="append",
        default=[],
        help="Limit to outputs/<name> (repeatable). Default: all with sequences.fa.",
    )
    parser.add_argument(
        "--max-attempts",
        type=int,
        default=8,
        help="Attempts per sequence on retryable GPU/queue failures (default: 8).",
    )
    parser.add_argument(
        "--backoff-s",
        type=float,
        default=30.0,
        help="Initial backoff seconds between retries (doubles each attempt).",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    outputs_dir = args.outputs_dir.expanduser().resolve()
    target_fasta = args.target_fasta.expanduser().resolve()
    fusion_fasta = (
        args.fusion_suffix_fasta.expanduser().resolve()
        if args.fusion_suffix_fasta is not None
        else target_fasta
    )

    if not target_fasta.exists():
        print(f"Target FASTA not found: {target_fasta}", file=sys.stderr)
        return 2
    if not fusion_fasta.exists():
        print(f"Fusion-suffix FASTA not found: {fusion_fasta}", file=sys.stderr)
        return 2
    if not outputs_dir.is_dir():
        print(f"Outputs directory not found: {outputs_dir}", file=sys.stderr)
        return 2

    target_name, target_sequence = load_target_sequence(target_fasta)
    _, fusion_suffix = load_target_sequence(fusion_fasta)
    print(
        f"Target {target_name}: {len(target_sequence)} aa from {target_fasta}",
        flush=True,
    )
    if fusion_fasta != target_fasta:
        print(
            f"Fusion suffix: {len(fusion_suffix)} aa from {fusion_fasta}",
            flush=True,
        )
    print(f"Scores file: {args.scores_name}", flush=True)

    if args.only:
        fasta_paths = [outputs_dir / name / "sequences.fa" for name in args.only]
    else:
        fasta_paths = sorted(outputs_dir.glob("*/sequences.fa"))

    missing = [p for p in fasta_paths if not p.exists()]
    if missing:
        for p in missing:
            print(f"Missing: {p}", file=sys.stderr)
        return 2
    if not fasta_paths:
        print(f"No sequences.fa under {outputs_dir}", file=sys.stderr)
        return 2

    client = LigandAIClient()
    total_scored = total_skipped = total_failed = 0
    for fa in fasta_paths:
        scored, skipped, failed = score_directory(
            client,
            fa,
            target_name,
            target_sequence,
            fusion_suffix=fusion_suffix,
            scores_name=args.scores_name,
            scorer=args.scorer,
            force=args.force,
            max_attempts=args.max_attempts,
            backoff_s=args.backoff_s,
        )
        total_scored += scored
        total_skipped += skipped
        total_failed += failed

    print(
        f"\nAll done: scored={total_scored} skipped={total_skipped} failed={total_failed}",
        flush=True,
    )
    return 1 if total_failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
