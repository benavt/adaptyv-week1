#!/usr/bin/env python3
"""Check catalog coverage, integrity, links, identities, and statistics sources."""
from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
import os
from pathlib import Path

from build_model_registry import ROOT, Registry, is_structure


def rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def validate(root, output, verify_hashes=True):
    root, output = root.resolve(), output.resolve()
    manifest = json.loads((output / "registry_manifest.json").read_text())
    errors = []
    models = rows(output / "models.csv")
    designs = {r["design_id"]: r for r in rows(output / "designs.csv")}
    sequences = {r["sequence_id"]: r for r in rows(output / "sequences.csv")}
    model_ids = {r["model_id"] for r in models}
    if manifest.get("status") != "complete":
        errors.append("Registry generation is incomplete.")
    if len(model_ids) != len(models):
        errors.append("Duplicate model IDs.")
    discovered = {p.relative_to(root).as_posix() for p in Registry(root, output).discover() if is_structure(p)}
    registered = {r["source_path"] for r in models}
    if discovered != registered:
        errors.append(f"Coverage mismatch: {len(discovered - registered)} missing, {len(registered - discovered)} stale structure files.")
    statistical_observations = 0
    ensemble_files = 0
    warning_counts = collections.Counter()
    species_conflicts = []
    for row in models:
        home = root / row["home"]
        meta = json.loads((home / "metadata.json").read_text())
        warning_counts.update(meta["warnings"])
        stats = json.loads((home / "statistics.json").read_text())
        source = root / row["source_path"]
        link = home / source.name
        if not link.is_symlink() or link.resolve() != source.resolve():
            errors.append(f"Wrong structure link: {row['model_id']}")
        if meta["model_id"] != row["model_id"] or stats["model_id"] != row["model_id"]:
            errors.append(f"Metadata identity mismatch: {row['model_id']}")
        if meta["design_id"]:
            did = meta["design_id"]
            if did not in designs:
                errors.append(f"Unregistered design: {did}")
            else:
                binder = [c for c in meta["chains"] if c["role"] == "binder"]
                # A target-only preparation can belong to a design without containing its binder.
                if binder and any(c["sequence_id"] != designs[did]["sequence_id"] for c in binder):
                    errors.append(f"Binder/design sequence mismatch: {row['model_id']}")
        for chain in meta["chains"]:
            sid = chain["sequence_id"]
            if sid not in sequences or chain["sequence"] != sequences[sid]["sequence"]:
                errors.append(f"Chain sequence missing/conflicting: {row['model_id']}")
        for parent in meta["provenance"]["parent_model_ids"]:
            if parent["model_id"] not in model_ids:
                errors.append(f"Missing parent: {row['model_id']}")
        for observation in stats["observations"]:
            if not (root / observation["source_path"]).is_file():
                errors.append(f"Missing statistic source: {row['model_id']}")
            if "row_number" not in observation and "json_pointer" not in observation and "values" not in observation:
                errors.append(f"Statistic has no source locator or value: {row['model_id']}")
            declared = observation.get("values", {}).get("species") if isinstance(observation.get("values"), dict) else None
            if declared and meta["target_species"] and declared.capitalize() != meta["target_species"]:
                species_conflicts.append({"model_id": row["model_id"], "statistics_source": observation["source_path"], "row_number": observation.get("row_number"), "recorded_species": declared, "sequence_verified_species": meta["target_species"]})
        statistical_observations += len(stats["observations"])
        ensemble_files += len(meta["coordinate_model_numbers"]) > 1
    for did, row in designs.items():
        home = root / row["home"]
        meta = json.loads((home / "metadata.json").read_text())
        if meta["sequence_id"] not in sequences:
            errors.append(f"Missing design sequence: {did}")
        for mid in meta["model_ids"]:
            link = home / "models" / mid
            if mid not in model_ids or not link.is_symlink() or not link.exists():
                errors.append(f"Broken design/model link: {did}/{mid}")
    broken_links = []
    for directory, subdirs, files in os.walk(output, followlinks=False):
        for name in subdirs + files:
            p = Path(directory) / name
            if p.is_symlink() and (not p.exists() or Path(os.readlink(p)).is_absolute()):
                broken_links.append(str(p.relative_to(output)))
    if broken_links:
        errors.append(f"Broken or non-relative links: {broken_links[:10]}")
    checked_sources = 0
    if verify_hashes:
        for row in rows(output / "source_files.csv"):
            p = root / row["source_path"]
            if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest() != row["sha256"]:
                errors.append(f"Source content changed since indexing: {row['source_path']}")
            checked_sources += 1
    report = {"status": "passed" if not errors else "failed", "models": len(models), "designs": len(designs), "unique_sequences": len(sequences), "discovered_structure_files": len(discovered), "statistics_observations": statistical_observations, "ensemble_files": ensemble_files, "source_hashes_checked": checked_sources, "broken_links": len(broken_links), "warning_counts": dict(warning_counts), "statistics_species_conflicts": species_conflicts, "errors": errors, "note": "Passing verifies registry integrity; recorded scientific caveats and empty artifacts remain flagged, not scientifically validated or silently corrected."}
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--skip-source-hashes", action="store_true")
    args = parser.parse_args()
    output = args.output or args.root / "registry"
    report = validate(args.root, output, not args.skip_source_hashes)
    (output / "validation.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    raise SystemExit(0 if report["status"] == "passed" else 1)


if __name__ == "__main__":
    main()
