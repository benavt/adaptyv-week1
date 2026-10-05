#!/usr/bin/env python3
"""Snapshot campaign recovery and scoring coverage without launching scorers.

Reads curated manifests and original status evidence, never assigns biological
success from a scheduler state or a numeric display field. Remote facts require
an explicitly collected observation; this script does not contact Mimas.
"""
import csv
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "03_Filtering/workflows/campaign_completion_20261005"
REPORTS = ROOT / "03_Filtering/workflows/campaign_top10_egfr_20261005/reports.json"
FLAGS = {"failed_before_prediction", "prepared_not_submitted"}


def evidence(path):
    p = ROOT / path
    return {"path": path, "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}


def native_status(run):
    path = ROOT / run["status_evidence"]
    if not run.get("status_locator"):
        return run.get("recorded_status", "unknown"), {}
    doc = json.loads(path.read_text())
    value = doc
    for token in run["status_locator"].split("/")[1:]:
        value = value[token.replace("~1", "/").replace("~0", "~")]
    if isinstance(value, list):
        statuses = {r.get("status", "unknown") for r in value}
        status = next(iter(statuses)) if len(statuses) == 1 else "mixed_outcomes" if statuses else "unknown"
    else:
        status = value
    return status, doc


def numeric(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def link(base, target, label=None):
    import os
    return f"[{label or target}](<{os.path.relpath(ROOT / target, base)}>)"


def main():
    BASE.mkdir(parents=True, exist_ok=True)
    reports = json.loads(REPORTS.read_text())
    observed = json.loads((BASE / "mimas_observation.json").read_text())
    recovery = []
    all_runs = []
    for manifest in sorted((ROOT / "03_Filtering/campaigns").glob("*/manifest.json")):
        campaign = json.loads(manifest.read_text())
        for run in campaign.get("runs", []):
            status, doc = native_status(run)
            row = dict(campaign_id=campaign["campaign_id"], run_id=run["run_id"], recorded_status=status,
                       status_locator=run.get("status_locator", ""), source=evidence(run["status_evidence"]),
                       superseded_by=run.get("superseded_by"), completed_models=doc.get("completed_models"),
                       failed_requests=doc.get("failed_requests"), reason=doc.get("reason"))
            all_runs.append(row)
            if status not in FLAGS:
                continue
            if not run.get("status_locator"):
                row["reason"] = (ROOT / run["status_evidence"]).read_text().strip()
            if run.get("supersession_evidence"):
                row["supersession_source"] = evidence(run["supersession_evidence"])
            if status == "failed_before_prediction":
                row["reason"] = "; ".join(sorted({v.get("error", "unrecorded") for v in doc.get("outcomes", [])}))
            row["analysis_complete"] = False
            row["recovery_run"] = observed["run_id"]
            row["recovery_state"] = observed["state"]
            row["next_action"] = "Collect and validate the already-generated retry outputs if pursuing AFSample3; reuse the existing paired protein-only Boltz2 ranking and score the five selected sequence pairs with pKa and DeltaForge."
            recovery.append(row)
    coverage = []
    for cid, report in reports.items():
        selected = report["top_models"]
        models = [s for pair in selected for s in pair["species"].values()]
        df_count = sum(numeric(m["metrics"].get("DeltaForge_dG")) and numeric(m["metrics"].get("DeltaForge_Kd_nM")) for m in models)
        pka_count = sum(all(numeric(m["metrics"].get(k)) for k in ["HIS37_pKa", "HIS100_pKa"]) for m in models)
        row = {"campaign_id": cid, "paired_candidates": report["paired_count"], "selected_sequences": len(selected),
               "selected_complexes": len(models), "reported_deltaforge_complexes": df_count,
               "reported_histidine_pka_complexes": pka_count,
               "binding_ranking": "available" if selected else "pending_linker_and_EGFR_folding",
               "analysis_completion": "not_established", "source": evidence(str(REPORTS.relative_to(ROOT))),
               "selected": [{"design_id": x["id"], "label": x["label"], "mean_ipSAE": x["mean_ipSAE"],
                              "species": {s: {k: v[k] for k in ["model", "model_sha256", "binder_sequence", "egfr_sequence", "binder_chain", "egfr_chain"]} for s, v in x["species"].items()}} for x in selected]}
        if cid == "dd2_rf3_top5pct_20261003":
            row["next_action"] = "Finish 100 RFD3 diffusions for each of the ten distinct parents, then the queued MPNN retry. Verify one-chain DD2-linker-DD1, exact domains, 5-20 residues and H/K acceptance. Fold accepted fusions separately with Human and Mouse EGFR; rank by Boltz2 mean ipSAE, then run pKa and DeltaForge."
        elif cid == "egfr_d3_references_md":
            row["next_action"] = "Resolve the reported Human glycan bond geometry issue and choose chemistry-aware preparation retaining the original glycans and binder disulfides. The current protein-only pKa batch runner rejects these inputs. Score both species with pKa and DeltaForge on documented supported chemistry."
        elif cid == "boltz2_human_rfd3_egfr_20261002":
            row["next_action"] = "Use the updated 10,000-pair ranking to select ten sequences (20 complexes). Audit retained JustHISpKa checkpoints by exact model hash. Prepare a fresh focused pKa/DeltaForge batch. The suspended 17,386-complex batch has deliberately removed PROPKA outputs; resuming it regenerates them."
        elif cid == "dd1_062_propka_20261004":
            row["next_action"] = "Reuse the validated Human bound PROPKA/JustHISpKa result with its geometry warnings and exact hashes; add Mouse bound pKa. Retain existing paired DeltaForge provenance and publish the joined two-species result."
        elif df_count == len(models) and pka_count == len(models):
            row["next_action"] = "Validate existing DeltaForge raw response, chain mapping and input hashes. Historical displayed JustHISpKa values are receptor-only controls; obtain or reuse validated bound-complex pKa results before certifying analysis completion. Keep receptor-only controls labelled."
        else:
            row["next_action"] = "Use the selected native Human/Mouse complexes for a focused pKa and DeltaForge batch, reusing only matching validated results. Join scores by stable design ID, species and exact input hash; resolve failed or conditional-only pKa sites."
        coverage.append(row)
        dest = ROOT / "03_Filtering/campaigns" / cid / "ANALYSIS_COMPLETION.md"
        base = dest.parent
        content = f"# Analysis completion: {cid}\n\nSnapshot: {datetime.now(timezone.utc).isoformat()}. **Analysis completion: not established.**\n\n"
        content += "The user's completion definition requires top binding sequences identified by Boltz2 ipSAE, with pKa and DeltaForge predictions for those sequences. "
        content += "For this paired-species pass, use mean Human/Mouse ipSAE and ten distinct binder sequences, or all eligible sequences if fewer than ten. This is a completion criterion, not evidence of experimentally confirmed binding.\n\n"
        content += f"| Gate | Source-backed snapshot |\n| --- | --- |\n| Paired Boltz2 candidates | {row['paired_candidates']} |\n| Selected sequences / complexes | {len(selected)} / {len(models)} |\n| DeltaForge ΔG and Kd displayed | {df_count}/{len(models)} complexes |\n| HIS37/HIS100 pKa displayed | {pka_count}/{len(models)} complexes; scope and preparation must be checked |\n| Bound pKa and provenance acceptance | Not yet certified for the entire selected cohort |\n\n"
        content += f"Next action: {row['next_action']}\n\n"
        content += "Numeric coverage is a report-field audit, not fresh scorer execution or certification. Historical receptor-only pKa values do not establish a binder effect. Completed stage records retain their original meaning; the final analysis criterion is tracked separately.\n\n"
        content += link(base, "docs/CAMPAIGN_ANALYSIS_COMPLETION.md", "Completion policy and recovery tracker") + " · "
        content += link(base, str(REPORTS.relative_to(ROOT)), "Frozen ranking and source structures") + "\n"
        dest.write_text(content)
    sources = [evidence("registry/campaign_runs.csv"), evidence(str(REPORTS.relative_to(ROOT))),
               evidence(str((BASE / "mimas_observation.json").relative_to(ROOT))),
               evidence("03_Filtering/JustHIpKA/run_top20_dd1_hispka.py"),
               evidence("03_Filtering/JustHIpKA/run_y32d_hispka.py"),
               evidence("03_Filtering/JustHIpKA/runs/2026-10-04_DD1_062_Human_JustHISpKa_bound_002/RESULTS.md"),
               evidence("01_Staging/2026-10-04-DD1_DD2_Stitch_EGFR_AFS3_submit_004/RUN_RECORD.md"),
               evidence("01_Staging/2026-10-04-DD1_DD2_Stitch_EGFR_AFS3_submit_005/RUN_RECORD.md"),
               evidence("01_Staging/runs/2026-10-04_human_rfd3_rank1_pka_batch_002/manifest.json"),
               evidence("03_Filtering/campaigns/boltz2_human_rfd3_egfr_20261002/analyses/pka_batch/2026-10-04_human_rfd3_rank1_pka_batch_002/propka_removal_20261005T090922.json")]
    out = {"schema_version": 1, "snapshot_at": datetime.now(timezone.utc).isoformat(),
           "policy": "docs/CAMPAIGN_ANALYSIS_COMPLETION.md", "sources": sources,
           "flagged_runs": recovery, "all_recorded_runs": all_runs, "campaign_coverage": coverage,
           "remote_observation": observed,
           "additional_campaign": {"campaign_id": "ligand_ai_egfr_20261005", "analysis_completion": "not_established", "next_action": "Collect and validate submitted paired Boltz2 outputs; establish full candidate ranking, then run pKa and DeltaForge on the top ten distinct sequence pairs. Latest prediction completeness is not assessed in this snapshot."}}
    (BASE / "audit.json").write_text(json.dumps(out, indent=2) + "\n")
    with (BASE / "recovery_runs.csv").open("w", newline="") as f:
        fields = ["campaign_id", "run_id", "recorded_status", "reason", "superseded_by", "recovery_run", "recovery_state", "analysis_complete", "next_action"]
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(recovery)
    print(json.dumps({"flagged_runs": len(recovery), "campaigns": len(coverage), "analysis_complete_certified": 0, "audit": str((BASE / "audit.json").relative_to(ROOT))}))


if __name__ == "__main__":
    main()
