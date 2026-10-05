#!/usr/bin/env python3
"""Build campaign navigation without moving scientific files or changing IDs.

Sources are curated 03_Filtering/campaigns/*/manifest.json files. All paths in
these manifests are project-relative. --check validates references and verifies
that generated pages/indexes match the sources; it writes nothing.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import re
import sys
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]
CAMPAIGN_COLUMNS = ["campaign_id", "title", "home", "manifest_source", "source_sha256", "current_outputs", "notes"]
RUN_COLUMNS = ["campaign_id", "run_id", "phase", "run_path", "recorded_status", "status_evidence", "status_locator", "superseded_by", "output_path", "prepared_requests", "submitted_requests", "failed_requests", "completed_models"]
LANDING_PAGES = ["README.md", "AGENTS.md", "01_Staging/README.md", "01_Staging/references/README.md", "01_Staging/runs/README.md", "01_Staging/EGFR/README.md", "01_Staging/Putative_DD1/README.md", "01_Staging/frozen-tiger-ice/README.md", "03_Filtering/README.md", "03_Filtering/Refolding/README.md", "03_Filtering/workflows/README.md", "03_Filtering/DD1_ASP_20261003/README.md", "03_Filtering/PyMol_top20/README.md", "03_Filtering/JustHIpKA/README.md", "tools/README.md", "results/README.md", "docs/PROJECT_ORGANIZATION.md", "docs/STRUCTURE_INFERENCE.md", "docs/MIMAS_SKILLS.md", "docs/WORKFLOW_SPEC.md", "docs/remote-sync.md", "registry/README.md"]


def project_path(root, value):
    p = PurePosixPath(value)
    if not value or p.is_absolute() or ".." in p.parts or "\\" in value:
        raise ValueError(f"Expected a project-relative path: {value!r}")
    result = root / p
    if not result.resolve().is_relative_to(root.resolve()):
        raise ValueError(f"Reference escapes project: {value}")
    return result


def existing(root, value):
    p = project_path(root, value)
    if not p.exists():
        raise ValueError(f"Missing navigation reference: {value}")
    return p


def pointer(document, locator):
    if not locator.startswith("/"):
        raise ValueError(f"Expected JSON pointer: {locator}")
    for part in locator[1:].split("/"):
        key = part.replace("~1", "/").replace("~0", "~")
        document = document[int(key)] if isinstance(document, list) else document[key]
    return document


def markdown_link(label, target, page):
    import os
    relative = os.path.relpath(target, page.parent).replace("\\", "/")
    return f"[{label}](<{relative}>)"


def cell(value):
    return str(value).replace("|", "\\|").replace("\n", " ")


def csv_text(rows, columns):
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue()


def load_campaigns(root):
    campaigns, seen = [], set()
    for source in sorted((root / "03_Filtering/campaigns").glob("*/manifest.json")):
        campaign = json.loads(source.read_text())
        cid = campaign["campaign_id"]
        if campaign.get("schema_version") != 1 or cid != source.parent.name or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", cid):
            raise ValueError(f"Invalid campaign identity/schema: {source}")
        if cid in seen:
            raise ValueError(f"Duplicate campaign ID: {cid}")
        seen.add(cid)
        for key in ("title", "summary", "locations", "current_outputs", "notes"):
            if key not in campaign:
                raise ValueError(f"Missing {key} in {source}")
        for entry in campaign["locations"] + campaign["current_outputs"]:
            existing(root, entry["path"])
            existing(root, entry["evidence_path"])
        for entry in campaign.get("documented_unavailable", []):
            existing(root, entry["evidence_path"])
            if project_path(root, entry["expected_path"]).exists():
                raise ValueError(f"Previously unavailable output now exists; review authority: {entry['expected_path']}")
        for run in campaign.get("runs", []):
            existing(root, run["path"])
            existing(root, run["status_evidence"])
            if run.get("output_path"):
                existing(root, run["output_path"])
        campaigns.append((source, campaign))
    if not campaigns:
        raise ValueError("No campaign manifests found under 03_Filtering/campaigns/")
    return campaigns


def run_rows(root, source, campaign):
    runs = campaign.get("runs", [])
    by_id = {r["run_id"]: r for r in runs}
    if len(by_id) != len(runs):
        raise ValueError(f"Duplicate run ID in {source}")
    rows = []
    for run in runs:
        evidence = existing(root, run["status_evidence"])
        document = json.loads(evidence.read_text()) if evidence.suffix == ".json" else None
        locator = run.get("status_locator", "")
        status = pointer(document, locator) if locator else run.get("recorded_status", "unknown")
        if isinstance(status, list):
            if any(not isinstance(outcome, dict) or not isinstance(outcome.get("status"), str) for outcome in status):
                raise ValueError(f"Invalid per-job status evidence: {source}/{run['run_id']}")
            statuses = {outcome["status"] for outcome in status}
            status = next(iter(statuses)) if len(statuses) == 1 else "mixed_outcomes" if statuses else "unknown"
        if not isinstance(status, str):
            raise ValueError(f"Status must be a string: {source}/{run['run_id']}")
        successor = run.get("superseded_by", "")
        if successor:
            if successor not in by_id:
                raise ValueError(f"Unknown successor: {successor}")
            proof = existing(root, run["supersession_evidence"])
            if successor not in proof.read_text():
                raise ValueError(f"Successor not named in evidence: {proof}")
        trail, current = set(), run["run_id"]
        while current:
            if current in trail:
                raise ValueError(f"Supersession cycle: {source}/{current}")
            trail.add(current)
            current = by_id[current].get("superseded_by", "")
        row = {"campaign_id": campaign["campaign_id"], "run_id": run["run_id"], "phase": run["phase"], "run_path": run["path"], "recorded_status": status, "status_evidence": run["status_evidence"], "status_locator": locator, "superseded_by": successor, "output_path": run.get("output_path", "")}
        for key in RUN_COLUMNS[-4:]:
            row[key] = document.get(key, "") if isinstance(document, dict) else ""
        rows.append(row)
    return rows


def campaign_page(root, source, campaign, runs):
    page = source.parent / "README.md"
    link = lambda label, path: markdown_link(label, root / path, page)
    lines = [f"# {campaign['title']}", "", campaign["summary"], "", f"Campaign ID: `{campaign['campaign_id']}`. Paths and report authority come from the linked evidence; recorded execution status is a local observation, not a live cluster check.", "", f"Navigation source: {link('manifest.json', source.relative_to(root).as_posix())}. Regenerate with `python3 scripts/build_project_navigation.py` from the project root; edit the manifest rather than this generated page.", "", "## Locations", "", "| Purpose | Location | Evidence |", "| --- | --- | --- |"]
    metadata = source.parent / "CAMPAIGN_METADATA.md"
    if metadata.is_file():
        lines[6:6] = [f"Run metadata and unfinished work: {markdown_link('CAMPAIGN_METADATA.md', metadata, page)}. This curated snapshot includes targets, parameters, jobs, filtering, evidence hashes and remaining work.", ""]
    for entry in campaign["locations"]:
        lines.append(f"| {cell(entry['role'])} | {link(entry['label'], entry['path'])} | {link('source', entry['evidence_path'])} |")
    lines += ["", "## Current available outputs", ""]
    if campaign["current_outputs"]:
        lines += ["These are available navigation targets. Their scientific scope and validation remain those recorded by their source documents.", "", "| Output | File | Authority evidence |", "| --- | --- | --- |"]
        for entry in campaign["current_outputs"]:
            lines.append(f"| {cell(entry['label'])} | {link(Path(entry['path']).name, entry['path'])} | {link('source', entry['evidence_path'])} |")
    local_links = []
    for group in ("run_display", "pptx"):
        directory = page.parent / group
        if directory.is_dir():
            for child in sorted(directory.iterdir()):
                if child.is_symlink() and child.exists():
                    local_links.append((group, child))
    if local_links:
        lines += ["", "## Campaign-local links", "", "These symlinks provide campaign-local entry points while preserving the original artifact paths.", "", "| Type | Campaign link | Original target |", "| --- | --- | --- |"]
        for group, child in local_links:
            target = child.resolve()
            lines.append(f"| {group} | {markdown_link(child.name, child, page)} | `{target.relative_to(root)}` |")
    if not campaign["current_outputs"] and not local_links:
        lines += ["No finalized output is designated by this navigation manifest."]
    if runs:
        lines += ["", "## Recorded runs and attempts", "", "| Run / attempt | Phase | Recorded status | Superseded by | Evidence |", "| --- | --- | --- | --- | --- |"]
        for run in runs:
            lines.append(f"| {link(run['run_id'], run['run_path'])} | {cell(run['phase'])} | {cell(run['recorded_status'])} | {cell(run['superseded_by'] or 'N/A')} | {link('source', run['status_evidence'])} `{run['status_locator']}` |")
        lines += ["", "Preparation, transfer, accepted submission, collection, and verification are separate outcomes. Blank counts in the registry attempt index mean unrecorded; zero means explicitly recorded as zero."]
    if campaign.get("documented_unavailable"):
        lines += ["", "## Documented outputs unavailable locally", "", "These files are named in existing documentation but were absent when this index was built. They are not designated as available or complete.", ""]
        for entry in campaign["documented_unavailable"]:
            lines.append(f"- `{entry['expected_path']}` — {link('documentation', entry['evidence_path'])}.")
    lines += ["", "## Scope and historical exceptions", ""]
    lines += [f"- {note}" for note in campaign["notes"]]
    lines += ["", f"Browse {link('all campaigns', '03_Filtering/campaigns/README.md')}, {link('staging', '01_Staging/README.md')}, or {link('model/design registry', 'registry/README.md')}.", ""]
    return "\n".join(lines)


def render(root):
    campaigns = load_campaigns(root)
    outputs, summary, attempts, evidence_hashes = {}, [], [], {}
    index = root / "03_Filtering/campaigns/README.md"
    lines = ["# Campaign navigation", "", "Start here to follow a biological question across staging, predictions, scoring, selections, and reports. Existing scientific files stay at their original paths. This generated index extends the model/design registry; it does not create new model IDs or duplicate coordinates.", "", "| Campaign | Scope | Run metadata |", "| --- | --- | --- |"]
    for source, campaign in campaigns:
        runs = run_rows(root, source, campaign)
        attempts.extend(runs)
        page = source.parent / "README.md"
        outputs[page.relative_to(root).as_posix()] = campaign_page(root, source, campaign, runs)
        metadata = source.parent / "CAMPAIGN_METADATA.md"
        metadata_link = markdown_link('Status, parameters and remaining work', metadata, index) if metadata.is_file() else 'Not prepared'
        lines.append(f"| {markdown_link(campaign['title'], page, index)} | {cell(campaign['summary'])} | {metadata_link} |")
        summary.append({"campaign_id": campaign["campaign_id"], "title": campaign["title"], "home": page.parent.relative_to(root).as_posix(), "manifest_source": source.relative_to(root).as_posix(), "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(), "current_outputs": ";".join(e["path"] for e in campaign["current_outputs"]), "notes": " ".join(campaign["notes"])})
        evidence = {source.relative_to(root).as_posix()}
        if metadata.is_file():
            evidence.add(metadata.relative_to(root).as_posix())
        evidence.update(e["evidence_path"] for e in campaign["locations"] + campaign["current_outputs"] + campaign.get("documented_unavailable", []))
        evidence.update(r["status_evidence"] for r in campaign.get("runs", []))
        evidence.update(r["supersession_evidence"] for r in campaign.get("runs", []) if r.get("supersession_evidence"))
        for rel in evidence:
            evidence_hashes[rel] = hashlib.sha256(existing(root, rel).read_bytes()).hexdigest()
    method = root / "docs/CAMPAIGN_METADATA_EXTRACTION.md"
    if method.is_file():
        lines += ["", markdown_link('Metadata extraction, snapshot inventory and update procedure', method, index) + ". Companions are curated snapshots; the navigation builder does not refresh their scientific status."]
        evidence_hashes[method.relative_to(root).as_posix()] = hashlib.sha256(method.read_bytes()).hexdigest()
    lines += ["", "Machine-readable indexes: [campaigns.csv](../../registry/campaigns.csv) and [campaign_runs.csv](../../registry/campaign_runs.csv).", "", "Sources are each campaign's `manifest.json`; regenerate pages and registry indexes with `python3 scripts/build_project_navigation.py`. Check all references and generated content with `python3 scripts/build_project_navigation.py --check`. Status is local evidence at rebuild time, not live monitoring.", ""]
    outputs["03_Filtering/campaigns/README.md"] = "\n".join(lines)
    outputs["registry/campaigns.csv"] = csv_text(summary, CAMPAIGN_COLUMNS)
    outputs["registry/campaign_runs.csv"] = csv_text(attempts, RUN_COLUMNS)
    outputs["registry/navigation_manifest.json"] = json.dumps({"schema_version": 1, "generator": "scripts/build_project_navigation.py", "campaigns": len(campaigns), "recorded_runs": len(attempts), "evidence_sha256": dict(sorted(evidence_hashes.items())), "policy": "Additive navigation only. Originals, model IDs, scores and runtime paths are unchanged. Status is recorded evidence, not live state."}, indent=2) + "\n"
    return outputs


def check_markdown(root, pages, planned=()):
    errors, checked = [], 0
    planned_paths = {project_path(root, p).resolve() for p in planned}
    for rel, content in pages.items():
        page = project_path(root, rel)
        # Ignore fenced examples; check actual inline links, including angle paths.
        content = re.sub(r"```.*?```", "", content, flags=re.S)
        for match in re.finditer(r"!?\[[^\]]*\]\((<[^>]*>|[^\s)]+)(?:\s+\"[^\"]*\")?\)", content):
            target = match.group(1).strip("<>")
            url = urlsplit(target)
            if url.scheme or url.netloc or not url.path:
                continue
            resolved = (page.parent / unquote(url.path)).resolve()
            checked += 1
            if not resolved.is_relative_to(root.resolve()) or (not resolved.exists() and resolved not in planned_paths):
                errors.append(f"Broken local link in {rel}: {target}")
    return checked, errors


def build(root, check=False):
    root = root.resolve()
    outputs = render(root)
    pages = {p: v for p, v in outputs.items() if p.endswith(".md")}
    for rel in LANDING_PAGES:
        p = project_path(root, rel)
        if p.is_file():
            pages[rel] = p.read_text()
    for p in (root / "03_Filtering/campaigns").glob("*/CAMPAIGN_METADATA.md"):
        pages[p.relative_to(root).as_posix()] = p.read_text()
    method = root / "docs/CAMPAIGN_METADATA_EXTRACTION.md"
    if method.is_file():
        pages[method.relative_to(root).as_posix()] = method.read_text()
    checked, errors = check_markdown(root, pages, outputs if not check else ())
    if check:
        for rel, content in outputs.items():
            p = project_path(root, rel)
            if not p.is_file() or p.read_text() != content:
                errors.append(f"Navigation output missing or stale: {rel}")
    if errors:
        raise ValueError("\n".join(errors))
    if not check:
        # Validate all inputs/links before writing any generated output.
        for rel, content in outputs.items():
            p = project_path(root, rel)
            if p.is_symlink():
                raise ValueError(f"Refusing to overwrite a navigation symlink: {rel}")
        for rel, content in outputs.items():
            p = project_path(root, rel)
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text(content)
    return {"status": "passed", "generated_files": len(outputs), "local_markdown_links_checked": checked, "mode": "check" if check else "build"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--check", action="store_true", help="Validate paths and generated content without writing files")
    args = parser.parse_args()
    try:
        print(json.dumps(build(args.root, args.check), indent=2))
    except (OSError, ValueError, KeyError, IndexError, TypeError) as exc:
        print(f"Navigation validation failed: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
