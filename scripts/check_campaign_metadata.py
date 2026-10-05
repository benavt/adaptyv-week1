#!/usr/bin/env python3
"""Check curated campaign Markdown, evidence freshness and local links, read-only.

Hashes establish snapshot freshness, not scientific validity. Directory counts
and remote state must be re-extracted as described in the shared method document.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re
import sys

from build_project_navigation import ROOT, check_markdown, project_path

SECTIONS = (
    "Scope and targets", "Inputs and parameters", "Run and job ledger",
    "Filtering and selection", "Available outputs and validation",
    "Unfinished work and unknowns", "Evidence catalog",
)
EVIDENCE = re.compile(
    r"^\| (E\d{2}) \| \[([^\]]+)\]\(<([^>]+)>\) \| (.*?) \| `([0-9a-f]{64})` \|$",
    re.M,
)


def check(root: Path) -> dict:
    root = root.resolve()
    manifests = sorted((root / "03_Filtering/campaigns").glob("*/manifest.json"))
    errors, pages, checked_sources = [], {}, 0
    if not manifests:
        errors.append("No campaign manifests found")
    for manifest in manifests:
        campaign_id = json.loads(manifest.read_text())["campaign_id"]
        page = manifest.parent / "CAMPAIGN_METADATA.md"
        if not page.is_file():
            errors.append(f"Missing metadata: {page.relative_to(root)}")
            continue
        content = page.read_text()
        rel = page.relative_to(root).as_posix()
        pages[rel] = content
        if f"| Campaign ID | `{campaign_id}` |" not in content:
            errors.append(f"Campaign ID mismatch: {rel}")
        headings = re.findall(r"^## (.+)$", content, re.M)
        if headings != list(SECTIONS):
            errors.append(f"Missing or nonuniform sections: {rel}")
        evidence = EVIDENCE.findall(content)
        identifiers = [row[0] for row in evidence]
        if not evidence or len(set(identifiers)) != len(identifiers):
            errors.append(f"Missing/duplicate evidence identifiers: {rel}")
        unknown = set(re.findall(r"\bE\d{2}\b", content)) - set(identifiers)
        if unknown:
            errors.append(f"Undefined evidence in {rel}: {sorted(unknown)}")
        for identifier, source, target, locator, expected_hash in evidence:
            try:
                path = project_path(root, source)
                if (page.parent / target).resolve() != path.resolve():
                    errors.append(f"Evidence link/label mismatch: {rel}/{identifier}")
                actual = hashlib.sha256(path.read_bytes()).hexdigest()
                checked_sources += 1
                if actual != expected_hash:
                    errors.append(f"Changed evidence; re-extract metadata: {rel}/{identifier}: {source}")
                if not locator.strip():
                    errors.append(f"Missing evidence locator: {rel}/{identifier}")
            except (OSError, ValueError) as exc:
                errors.append(f"Invalid evidence: {rel}/{identifier}: {exc}")
    method = root / "docs/CAMPAIGN_METADATA_EXTRACTION.md"
    if method.is_file():
        pages[method.relative_to(root).as_posix()] = method.read_text()
    else:
        errors.append("Missing metadata extraction documentation")
    links, link_errors = check_markdown(root, pages)
    errors.extend(link_errors)
    if errors:
        raise ValueError("\n".join(errors))
    return {"status": "passed", "campaigns": len(manifests),
            "evidence_file_hashes_checked": checked_sources,
            "local_markdown_links_checked": links,
            "validation_level": "documentation structure, source freshness and links; no scientific or remote runtime validation"}


if __name__ == "__main__":
    try:
        print(json.dumps(check(ROOT), indent=2))
    except (OSError, ValueError, KeyError) as exc:
        print(str(exc), file=sys.stderr)
        sys.exit(1)
