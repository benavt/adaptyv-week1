#!/usr/bin/env python3
"""Score the four BRD2ET hotspot designs on the panel folds already run.

Downloads both replicates of each sequence against BRD2ET, BRD3ET, BRD4ET,
and BRDTET, DeltaForge-scores the fold handles, and appends four rows to the
cross-isoform Kd table.
"""

from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
LOG = REPO / "design/ligand_ai/logs/BRD2ET-20260925-114308.log"
WORK = REPO / "design/ligand_ai/workspaces/brd2et_msa_hotspots_20260925"
OUT = REPO / "analysis/outputs/ligand_ai/isoform_fingerprint"
SUMMARY = OUT / "top_10pct_cross_isoform_kd.csv"
CHECKPOINT = OUT / "BRD2ET_msa_hotspots" / "fold_kd_checkpoint.json"
SESSION = "session_parallel_1790350989992_f651cf5f"
RUN = "BRD2ET_msa_hotspots"
ORDER = (
    "LIGTHRALIEVRLKSTTKEGRTLVDKEEKTV",
    "GPEGTIGLCRESLDKKRYNAKEC",
    "IEEVNETQCQFEVKKLDRSKFSGYIGPCTGF",
    "SCRLCEPKRNYNASFRDDIALPFL",
)
PREFIX = {
    "BRD2ET": "Brd2ET",
    "BRD3ET": "BRD3ET",
    "BRD4ET": "Brd4ET",
    "BRDTET": "BrdTET",
}


def panel_legs() -> list[dict]:
    text = LOG.read_text(encoding="utf-8")
    marker = "--- off-target panel ---"
    start = text.index(marker) + len(marker)
    blob = text[start:].lstrip()
    payload, _ = json.JSONDecoder().raw_decode(blob)
    return payload["legs"]


def as_float(value):
    if isinstance(value, bool) or value is None or value == "":
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def score_handle(client, handle: str) -> dict:
    scored = client.deltaforge.batch_score_fold([handle])
    rows = scored.get("results") if isinstance(scored, dict) else None
    row = rows[0] if isinstance(rows, list) and rows else (scored if isinstance(scored, dict) else {})
    metrics = row.get("metrics") if isinstance(row.get("metrics"), dict) else {}
    kd = as_float(metrics.get("kd_nm"))
    if kd is None:
        kd = as_float(row.get("kd_nm"))
    if kd is None:
        kd = as_float(row.get("predicted_kd"))
    delta_g = as_float(metrics.get("delta_g"))
    if delta_g is None:
        delta_g = as_float(row.get("delta_g"))
    reason = None
    if kd is None:
        reason = (
            row.get("affinity_reason")
            or row.get("kd_reason")
            or row.get("error")
            or metrics.get("kd_reason")
            or "no kd_nm"
        )
    return {
        "kd_nm": kd,
        "delta_g": delta_g,
        "raw": scored,
        "error": reason,
    }


def best_design(designs: list[dict], sequence: str) -> dict | None:
    chosen = None
    for design in designs:
        if design.get("sequence") not in (None, sequence):
            continue
        if design.get("coord_degenerate") and chosen is not None:
            continue
        ipsae = as_float(design.get("peptide_ipsae"))
        current = as_float(chosen.get("peptide_ipsae")) if chosen else None
        if chosen is None or (ipsae is not None and (current is None or ipsae >= current)):
            chosen = design
    return chosen


def pdb_text(client, job_id: str, fold_result_id) -> str:
    from ligandai.errors import LigandAIError

    try:
        resp = client.transport.request(
            "GET",
            "/api/folding/jobs/%s/structure/pdb" % job_id,
            expect_json=False,
        )
        text = resp.text if hasattr(resp, "text") else str(resp)
        if "ATOM" in text[:800]:
            return text
    except LigandAIError as exc:
        print("  folding-job pdb failed for %s: %s" % (job_id, exc), flush=True)
    if fold_result_id is None:
        raise RuntimeError("no PDB for %s" % job_id)
    return client.structures.get_pdb(fold_result_id)


def percent_plddt(value):
    number = as_float(value)
    if number is None:
        return None
    if number <= 1.5:
        return number * 100.0
    return number


def main() -> None:
    from ligandai import LigandAI
    from ligandai.errors import LigandAIError

    key = os.environ.get("LIGANDAI_API_KEY")
    if not key:
        sys.exit("set LIGANDAI_API_KEY")
    client = LigandAI(key, timeout=300.0)
    legs = panel_legs()
    records = []
    for leg in legs:
        sequence = leg["sequence"]
        receptor = leg["receptor"]
        rank = ORDER.index(sequence) + 1
        for rep, job_id in enumerate(leg["fold_job_ids"], start=1):
            print("recover %s %s rep %s %s" % (rank, receptor, rep, job_id), flush=True)
            error = None
            try:
                client.folds.recover(job_id, wait=True, timeout=60)
            except LigandAIError as exc:
                error = str(exc)
                print("  recover failed: %s" % error, flush=True)
            chosen = None
            try:
                listing = client.folds.list_designs(job_id)
                chosen = best_design(listing.get("designs") or [], sequence)
            except LigandAIError as exc:
                error = error or str(exc)
                print("  list_designs failed: %s" % exc, flush=True)
            if chosen is None:
                if error is None:
                    try:
                        info = client.jobs.get(job_id)
                        error = info.error_message or info.status or "no design"
                    except LigandAIError as exc:
                        error = str(exc)
                records.append({
                    "sequence": sequence,
                    "receptor": receptor,
                    "rank": rank,
                    "rep": rep,
                    "job_id": job_id,
                    "error": error or "no design",
                })
                _save(records)
                continue
            handle = chosen.get("fold_handle")
            fold_result_id = chosen.get("fold_result_id")
            ipsae = as_float(chosen.get("peptide_ipsae"))
            dest = WORK / "structures" / receptor / ("%d_%s_rep%d.pdb" % (rank, sequence[:8], rep))
            dest.parent.mkdir(parents=True, exist_ok=True)
            pdb = ""
            try:
                pdb = pdb_text(client, job_id, fold_result_id)
            except Exception as exc:  # noqa: BLE001
                error = error or str(exc)
                print("  pdb failed: %s" % exc, flush=True)
            if pdb:
                dest.write_text(pdb if pdb.endswith("\n") else pdb + "\n", encoding="utf-8")
                print("  wrote %s handle=%s ipsae=%s" % (dest.name, handle, ipsae), flush=True)
            scored = score_handle(client, handle) if handle else {
                "kd_nm": None, "delta_g": None, "raw": None, "error": "no fold_handle",
            }
            if scored["error"] is None:
                scored["error"] = chosen.get("kd_reason") if scored["kd_nm"] is None else None
            print("  kd_nm=%s error=%s" % (scored["kd_nm"], scored["error"]), flush=True)
            ptm = None
            raw = scored.get("raw")
            if isinstance(raw, dict) and raw.get("results"):
                ptm = as_float(raw["results"][0].get("ptm"))
            records.append({
                "sequence": sequence,
                "receptor": receptor,
                "rank": rank,
                "rep": rep,
                "job_id": job_id,
                "fold_handle": handle,
                "fold_result_id": fold_result_id,
                "peptide_ipsae": ipsae,
                "iptm": as_float(chosen.get("iptm")),
                "ptm": ptm,
                "plddt": percent_plddt(chosen.get("plddt_mean")),
                "pdb": str(dest.relative_to(REPO)) if pdb else "",
                "kd_nm": scored["kd_nm"],
                "delta_g": scored["delta_g"],
                "error": scored["error"] or "",
                "deltaforge": scored["raw"],
            })
            _save(records)
    append_csv(records)
    print("appended %s" % SUMMARY, flush=True)


def _save(records: list[dict]) -> None:
    CHECKPOINT.parent.mkdir(parents=True, exist_ok=True)
    CHECKPOINT.write_text(json.dumps(records, indent=1, default=str) + "\n", encoding="utf-8")


def better(left: dict, right: dict) -> dict:
    left_ipsae = left.get("peptide_ipsae")
    right_ipsae = right.get("peptide_ipsae")
    if left_ipsae is None:
        return right
    if right_ipsae is None:
        return left
    return left if left_ipsae >= right_ipsae else right


def append_csv(records: list[dict]) -> None:
    with SUMMARY.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)
    rows = [row for row in rows if row.get("run") != RUN]
    by_seq: dict[str, dict[str, dict]] = {}
    for record in records:
        by_seq.setdefault(record["sequence"], {})
        current = by_seq[record["sequence"]].get(record["receptor"])
        by_seq[record["sequence"]][record["receptor"]] = record if current is None else better(current, record)
    for rank, sequence in enumerate(ORDER, start=1):
        row = {name: "" for name in fieldnames}
        row.update({
            "gene": "BRD2ET",
            "session_id": SESSION,
            "rank": str(rank),
            "n_generated": "300",
            "sequence": sequence,
            "length": str(len(sequence)),
            "run": RUN,
            "model_id": "BRD2ET_msa_%d" % rank,
        })
        for receptor, prefix in PREFIX.items():
            chosen = by_seq.get(sequence, {}).get(receptor) or {}
            mapping = {
                "kd_nm": chosen.get("kd_nm"),
                "delta_g": chosen.get("delta_g"),
                "iptm": chosen.get("iptm"),
                "ptm": chosen.get("ptm"),
                "ipsae": chosen.get("peptide_ipsae"),
                "peptide_ipsae": chosen.get("peptide_ipsae"),
                "plddt": chosen.get("plddt"),
                "job_id": chosen.get("job_id"),
                "error": chosen.get("error") or "",
            }
            for field, value in mapping.items():
                column = "%s_%s" % (prefix, field)
                if column not in fieldnames:
                    continue
                if value is None:
                    row[column] = ""
                elif isinstance(value, float):
                    row[column] = format(value, ".8g")
                else:
                    row[column] = str(value)
        rows.append(row)
    tmp = SUMMARY.with_suffix(".csv.tmp")
    with tmp.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    tmp.replace(SUMMARY)


FILL_CHECKPOINT = OUT / "BRD2ET_msa_hotspots" / "fill_blanks_checkpoint.json"
BLANKS = (
    {"sequence": ORDER[0], "rank": 1, "receptor": "BRD2ET", "template": "fold_1790351404275_0yjuojs49"},
    {"sequence": ORDER[2], "rank": 3, "receptor": "BRD2ET", "template": "fold_1790351404275_0yjuojs49"},
    {"sequence": ORDER[1], "rank": 2, "receptor": "BRD4ET", "template": "fold_1790351405329_jz3y3bvzs"},
    {"sequence": ORDER[2], "rank": 3, "receptor": "BRD4ET", "template": "fold_1790351405329_jz3y3bvzs"},
    {"sequence": ORDER[2], "rank": 3, "receptor": "BRDTET", "template": "fold_1790351405882_d121km2v9"},
)


def _receptor_from_sibling(client, template_id: str) -> tuple[str, str]:
    info = client.jobs.get(template_id)
    payload = info.model_dump() if hasattr(info, "model_dump") else dict(info)
    data = payload.get("inputData") or payload.get("input_data") or {}
    entities = data.get("entities") or []
    receptor = next((item for item in entities if item.get("chainId") == "A"), None)
    if not receptor or not receptor.get("sequence"):
        raise RuntimeError("no chain A sequence on %s" % template_id)
    gene = receptor.get("geneName") or data.get("targetGeneName") or ""
    return receptor["sequence"], gene


def _fill_key(item: dict) -> str:
    return "%s|%s" % (item["sequence"], item["receptor"])


def _load_fill() -> dict:
    if not FILL_CHECKPOINT.exists():
        return {}
    rows = json.loads(FILL_CHECKPOINT.read_text(encoding="utf-8"))
    return {_fill_key(row): row for row in rows}


def _save_fill(rows: dict) -> None:
    FILL_CHECKPOINT.parent.mkdir(parents=True, exist_ok=True)
    FILL_CHECKPOINT.write_text(
        json.dumps(list(rows.values()), indent=1, default=str) + "\n",
        encoding="utf-8",
    )


def _wait_status(client, job_id: str, timeout: float = 3600.0) -> str:
    import time
    from ligandai.errors import LigandAIError

    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        try:
            info = client.jobs.get(job_id)
        except LigandAIError as exc:
            print("  status %s: %s" % (job_id, exc), flush=True)
            time.sleep(20)
            continue
        status = str(info.status or "").lower()
        print("  %s %s" % (job_id, status), flush=True)
        if status in {"completed", "complete", "failed", "cancelled", "error"}:
            return status
        time.sleep(20)
    raise TimeoutError(job_id)


def _score_stored(client, item: dict, row: dict) -> None:
    """Score a finished fold from its stored PDB when no design row exists yet."""
    from ligandai.errors import LigandAIError

    job_id = row["job_id"]
    body = client.transport.request("GET", "/api/folding/jobs/%s" % job_id) or {}
    result = body.get("result") if isinstance(body.get("result"), dict) else {}
    pdb = result.get("pdbContent") or result.get("pdb_content") or ""
    if not isinstance(pdb, str):
        pdb = ""
    if "ATOM" not in pdb:
        print("  inline pdb len=%d start=%r" % (len(pdb), pdb[:80]), flush=True)
        resp = client.transport.request(
            "GET",
            "/api/folding/jobs/%s/structure/pdb" % job_id,
            expect_json=False,
        )
        pdb = resp.text if hasattr(resp, "text") else str(resp)
    if "ATOM" not in pdb:
        raise RuntimeError("no ATOM records for %s (start %r)" % (job_id, pdb[:80]))
    dest = (
        WORK / "structures" / item["receptor"]
        / ("%d_%s_rep1.pdb" % (item["rank"], item["sequence"][:8]))
    )
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(pdb if pdb.endswith("\n") else pdb + "\n", encoding="utf-8")
    try:
        scored = client.deltaforge.score_pdb(
            pdb_content=pdb,
            receptor_chains=["A"],
            peptide_chain="B",
        )
    except LigandAIError as exc:
        row["error"] = str(exc)
        row["pdb"] = str(dest.relative_to(REPO))
        return
    kd = as_float(getattr(scored, "kd_nm", None))
    plddt_raw = result.get("mean_plddt")
    if plddt_raw is None:
        plddt_raw = result.get("plddt")
    if isinstance(plddt_raw, list) and plddt_raw:
        plddt_raw = sum(float(value) for value in plddt_raw) / len(plddt_raw)
    row.update({
        "peptide_ipsae": as_float(result.get("peptideIpsae") or result.get("peptide_ipsae") or result.get("ipsae")),
        "iptm": as_float(result.get("iptm")),
        "ptm": as_float(result.get("ptm")),
        "plddt": percent_plddt(plddt_raw),
        "pdb": str(dest.relative_to(REPO)),
        "kd_nm": kd,
        "delta_g": as_float(getattr(scored, "dg", None)),
        "error": "" if kd is not None else "no kd_nm",
    })


def _designs_or_recover(client, job_id: str) -> dict:
    import time
    from ligandai.errors import LigandAIError

    last = {}
    for attempt in range(6):
        try:
            last = client.folds.list_designs(job_id) or {}
        except LigandAIError as exc:
            print("  list_designs %s: %s" % (job_id, exc), flush=True)
            last = {}
        if last.get("designs"):
            return last
        print("  recovering %s (attempt %d)" % (job_id, attempt + 1), flush=True)
        try:
            client.folds.recover(job_id, wait=True, timeout=90)
        except LigandAIError as exc:
            print("  recover %s: %s" % (job_id, exc), flush=True)
        time.sleep(5)
    return last


def _patch_csv(updates: list[dict]) -> None:
    with SUMMARY.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)
    by_seq = {
        row["sequence"]: row
        for row in rows
        if row.get("run") == RUN
    }
    for update in updates:
        row = by_seq.get(update["sequence"])
        if row is None:
            raise RuntimeError("missing CSV row for %s" % update["sequence"])
        prefix = PREFIX[update["receptor"]]
        values = {
            "kd_nm": update.get("kd_nm"),
            "delta_g": update.get("delta_g"),
            "iptm": update.get("iptm"),
            "ptm": update.get("ptm"),
            "ipsae": update.get("peptide_ipsae"),
            "peptide_ipsae": update.get("peptide_ipsae"),
            "plddt": update.get("plddt"),
            "job_id": update.get("job_id"),
            "error": "" if update.get("kd_nm") is not None else (update.get("error") or ""),
        }
        for field, value in values.items():
            column = "%s_%s" % (prefix, field)
            if value is None:
                row[column] = ""
            elif isinstance(value, float):
                row[column] = format(value, ".8g")
            else:
                row[column] = str(value)
    tmp = SUMMARY.with_suffix(".csv.tmp")
    with tmp.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    tmp.replace(SUMMARY)


def fill_blanks() -> None:
    import time
    from ligandai import LigandAI
    from ligandai.errors import LigandAIError

    key = os.environ.get("LIGANDAI_API_KEY")
    if not key:
        sys.exit("set LIGANDAI_API_KEY")
    client = LigandAI(key, timeout=300.0)
    receptors: dict[str, tuple[str, str]] = {}
    for template in {item["template"] for item in BLANKS}:
        sequence, gene = _receptor_from_sibling(client, template)
        receptors[template] = (sequence, gene)
        print("template %s gene=%s length=%d" % (template, gene, len(sequence)), flush=True)

    state = _load_fill()
    for item in BLANKS:
        slot = _fill_key(item)
        row = state.get(slot) or dict(item)
        if row.get("kd_nm") is not None and row.get("pdb"):
            state[slot] = row
            continue
        if not row.get("job_id"):
            receptor_seq, gene = receptors[item["template"]]
            print("folding %s x %s" % (item["sequence"][:8], item["receptor"]), flush=True)
            job = client.fold(
                receptor_seq,
                item["sequence"],
                diffusion_samples=4,
                sampling_steps=50,
                msa_enabled=True,
                auto_score=False,
                force_resubmit=True,
                target_gene=gene or None,
            )
            row["job_id"] = job.id
            row["submitted_at"] = time.time()
            print("  submitted %s" % job.id, flush=True)
        state[slot] = row
        _save_fill(state)

    for item in BLANKS:
        slot = _fill_key(item)
        row = state[slot]
        if row.get("kd_nm") is not None and row.get("pdb"):
            continue
        job_id = row["job_id"]
        print("waiting %s %s" % (item["receptor"], job_id), flush=True)
        status = _wait_status(client, job_id)
        row["status"] = status
        if status not in {"completed", "complete"}:
            info = client.jobs.get(job_id)
            row["error"] = info.error_message or status
            row["kd_nm"] = None
            state[slot] = row
            _save_fill(state)
            print("  failed: %s" % row["error"], flush=True)
            continue
        try:
            listing = client.folds.list_designs(job_id) or {}
        except LigandAIError as exc:
            print("  list_designs %s: %s" % (job_id, exc), flush=True)
            listing = {}
        chosen = best_design(listing.get("designs") or [], item["sequence"])
        if chosen is None:
            _score_stored(client, item, row)
            state[slot] = row
            _save_fill(state)
            print(
                "  stored %s kd_nm=%s ipsae=%s"
                % (item["receptor"], row.get("kd_nm"), row.get("peptide_ipsae")),
                flush=True,
            )
            continue
        handle = chosen.get("fold_handle")
        fold_result_id = chosen.get("fold_result_id")
        dest = (
            WORK / "structures" / item["receptor"]
            / ("%d_%s_rep1.pdb" % (item["rank"], item["sequence"][:8]))
        )
        dest.parent.mkdir(parents=True, exist_ok=True)
        pdb = pdb_text(client, job_id, fold_result_id)
        dest.write_text(pdb if pdb.endswith("\n") else pdb + "\n", encoding="utf-8")
        scored = score_handle(client, handle) if handle else {
            "kd_nm": None, "delta_g": None, "raw": None, "error": "no fold_handle",
        }
        ptm = None
        raw = scored.get("raw")
        if isinstance(raw, dict) and raw.get("results"):
            ptm = as_float(raw["results"][0].get("ptm"))
        row.update({
            "fold_handle": handle,
            "fold_result_id": fold_result_id,
            "peptide_ipsae": as_float(chosen.get("peptide_ipsae")),
            "iptm": as_float(chosen.get("iptm")),
            "ptm": ptm,
            "plddt": percent_plddt(chosen.get("plddt_mean")),
            "pdb": str(dest.relative_to(REPO)),
            "kd_nm": scored["kd_nm"],
            "delta_g": scored["delta_g"],
            "error": "" if scored["kd_nm"] is not None else (scored["error"] or "no kd_nm"),
            "deltaforge": scored["raw"],
        })
        state[slot] = row
        _save_fill(state)
        print(
            "  %s kd_nm=%s ipsae=%s"
            % (dest.name, row["kd_nm"], row["peptide_ipsae"]),
            flush=True,
        )

    _patch_csv(list(state.values()))
    print("patched %s" % SUMMARY, flush=True)


def inspect_fills() -> None:
    from ligandai import LigandAI
    from ligandai.errors import LigandAIError

    key = os.environ.get("LIGANDAI_API_KEY")
    if not key:
        sys.exit("set LIGANDAI_API_KEY")
    client = LigandAI(key, timeout=60.0)
    state = _load_fill()
    for row in state.values():
        if row.get("kd_nm") is not None:
            continue
        job_id = row["job_id"]
        info = client.jobs.get(job_id)
        print(
            row["sequence"][:8], row["receptor"], job_id,
            "status", info.status, "error", info.error_message,
            flush=True,
        )
        try:
            body = client.transport.request("GET", "/api/folding/jobs/%s" % job_id) or {}
        except LigandAIError as exc:
            print("  folding", type(exc).__name__, str(exc)[:180], flush=True)
            continue
        result = body.get("result") if isinstance(body.get("result"), dict) else {}
        print(
            "  has", body.get("hasStructure") or body.get("has_structure"),
            "modal", body.get("modalStage"),
            "result_keys", list(result)[:8],
            flush=True,
        )


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "fill":
        fill_blanks()
    elif len(sys.argv) > 1 and sys.argv[1] == "inspect":
        inspect_fills()
    else:
        main()
