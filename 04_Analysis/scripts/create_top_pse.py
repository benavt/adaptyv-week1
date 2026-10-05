#!/usr/bin/env python3
"""PyMOL sessions for the unique top-0.5% models of each completed run.

One session per run loads the folded RF3 models, colors chain A purple and
chain B cyan, and highlights hotspot, avoid, and flexible residues from the
design specification. Input residue ids (A30) are mapped through
extra.sampled_contig onto the output chain and residue number.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

FOUNDRY = Path(__file__).resolve().parents[2] / "design" / "rfd3" / "foundry"
RUNS = FOUNDRY / "runs"
DEFAULT_CSV = FOUNDRY / "reports" / "top_0.5pct_sequences.csv"
PYMOL = Path("/opt/homebrew/bin/pymol")

BACKBONE_RE = re.compile(r"^(.*)_b\d+_d\d+$")
DIFFUSED_RE = re.compile(r"^(\d+)P$")
RANGE_RE = re.compile(r"^([A-Za-z]+)(\d+)-(\d+)$")
RESIDUE_RE = re.compile(r"^([A-Za-z]+)(\d+)$")
AVOID_KEYS = ("select_avoid", "select_negative_hotspots")
CHAIN_LETTERS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"


def backbone_id(model_id: str) -> str:
    match = BACKBONE_RE.match(model_id)
    return match.group(1) if match else model_id


def parse_sampled_contig(sampled: str) -> tuple[dict[str, tuple[str, int]], str]:
    """Map input residue ids onto (output chain, 1-based residue number).

    Components split by /0 become chains A, B, ... in order. A token like
    108P is that many diffused residues with no input id. A144 or A1-67 are
    fixed target residues numbered consecutively on that output chain.
    The target chain is the first component that contains a fixed residue.
    """
    components = [part.strip(" ,") for part in sampled.split("/0")]
    components = [part for part in components if part]
    mapping: dict[str, tuple[str, int]] = {}
    target_chain: str | None = None
    for index, component in enumerate(components):
        if index >= len(CHAIN_LETTERS):
            raise ValueError(f"Too many contig components in {sampled}")
        chain = CHAIN_LETTERS[index]
        resi = 1
        for token in (piece.strip() for piece in component.split(",")):
            if not token:
                continue
            diffused = DIFFUSED_RE.fullmatch(token)
            if diffused:
                resi += int(diffused.group(1))
                continue
            span = RANGE_RE.fullmatch(token)
            if span:
                letter, start_s, end_s = span.groups()
                start, end = int(start_s), int(end_s)
                step = 1 if end >= start else -1
                for number in range(start, end + step, step):
                    mapping[f"{letter}{number}"] = (chain, resi)
                    resi += 1
                target_chain = target_chain or chain
                continue
            single = RESIDUE_RE.fullmatch(token)
            if single:
                letter, number = single.group(1), int(single.group(2))
                mapping[f"{letter}{number}"] = (chain, resi)
                resi += 1
                target_chain = target_chain or chain
                continue
            raise ValueError(f"Unrecognized sampled_contig token {token!r} in {sampled}")
    if target_chain is None:
        target_chain = "B" if len(components) > 1 else "A"
    return mapping, target_chain


def is_bkbn(value: Any) -> bool:
    return isinstance(value, str) and value.strip().upper() == "BKBN"


def is_fully_flexible(value: Any) -> bool:
    if value is None:
        return True
    if isinstance(value, list):
        return len(value) == 0
    text = str(value).strip()
    return text in {"", "[]"}


def residue_map(value: Any) -> dict[str, Any]:
    if isinstance(value, dict):
        return value
    return {}


def find_spec_json(run_dir: Path, backbone: str) -> Path:
    direct = [
        run_dir / "rfd3" / "outputs" / f"{backbone}.json",
        RUNS / "NS1BctdRNA" / "outputs" / f"{backbone}.json",
    ]
    for path in direct:
        if path.is_file():
            return path
    hits = sorted(RUNS.glob(f"*/rfd3/outputs/{backbone}.json"))
    hits += sorted(RUNS.glob(f"*/outputs/{backbone}.json"))
    if not hits:
        raise FileNotFoundError(f"No design JSON for {backbone}")
    return hits[0]


def load_models(csv_path: Path, only: str | None) -> dict[str, list[str]]:
    by_run: dict[str, list[str]] = defaultdict(list)
    seen: dict[str, set[str]] = defaultdict(set)
    with csv_path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            run = row["run"]
            if only and run != only:
                continue
            model_id = row["model_id"]
            if model_id in seen[run]:
                continue
            seen[run].add(model_id)
            by_run[run].append(model_id)
    return dict(by_run)


def mapped_keys(
    keys: dict[str, Any],
    mapping: dict[str, tuple[str, int]],
    *,
    predicate: Any | None = None,
) -> list[tuple[str, int]]:
    chosen: list[tuple[str, int]] = []
    missing: list[str] = []
    for key, value in keys.items():
        if predicate is not None and not predicate(value):
            continue
        hit = mapping.get(key)
        if hit is None:
            missing.append(key)
            continue
        chosen.append(hit)
    if missing:
        print(f"  unmapped residues: {', '.join(missing)}", file=sys.stderr)
    return chosen


def selection(obj: str, residues: list[tuple[str, int]]) -> str:
    by_chain: dict[str, list[int]] = defaultdict(list)
    for chain, resi in residues:
        if resi not in by_chain[chain]:
            by_chain[chain].append(resi)
    parts = []
    for chain, numbers in by_chain.items():
        resi = "+".join(str(number) for number in sorted(numbers))
        parts.append(f"({obj} and chain {chain} and resi {resi})")
    return " or ".join(parts)


class ModelView:
    def __init__(
        self,
        model_id: str,
        cif: Path,
        target_chain: str,
        hotspots: list[tuple[str, int]],
        avoid: list[tuple[str, int]],
        flex_full: list[tuple[str, int]],
        flex_bkbn: list[tuple[str, int]],
    ) -> None:
        self.model_id = model_id
        self.cif = cif
        self.target_chain = target_chain
        self.hotspots = hotspots
        self.avoid = avoid
        self.flex_full = flex_full
        self.flex_bkbn = flex_bkbn


def build_model(run_dir: Path, model_id: str) -> ModelView:
    cif = run_dir / "rf3" / "outputs" / model_id / f"{model_id}_model.cif"
    if not cif.is_file():
        raise FileNotFoundError(cif)
    spec_path = find_spec_json(run_dir, backbone_id(model_id))
    spec = json.loads(spec_path.read_text())["specification"]
    sampled = spec["extra"]["sampled_contig"]
    mapping, target_chain = parse_sampled_contig(sampled)
    hotspots = mapped_keys(residue_map(spec.get("select_hotspots")), mapping)
    avoid_src: dict[str, Any] = {}
    for key in AVOID_KEYS:
        avoid_src.update(residue_map(spec.get(key)))
    avoid = mapped_keys(avoid_src, mapping)
    fixed = spec.get("select_fixed_atoms")
    flex_full = mapped_keys(residue_map(fixed), mapping, predicate=is_fully_flexible)
    flex_bkbn = mapped_keys(residue_map(fixed), mapping, predicate=is_bkbn)
    return ModelView(model_id, cif, target_chain, hotspots, avoid, flex_full, flex_bkbn)


def selection_parts(models: list[ModelView], attr: str) -> list[str]:
    parts = []
    for model in models:
        residues = getattr(model, attr)
        if residues:
            parts.append(selection(model.model_id, residues))
    return parts


def write_pml(run: str, models: list[ModelView], pml_path: Path, pse_path: Path) -> None:
    """One Python block so PyMOL's startup command queue cannot drop the save."""
    body = [
        "from pymol import cmd",
        "cmd.reinitialize()",
        "cmd.bg_color('white')",
        "cmd.set('ray_opaque_background', 1)",
        "cmd.set('cartoon_fancy_helices', 1)",
        "cmd.hide('everything')",
    ]
    for model in models:
        body.append(f"cmd.load({str(model.cif.resolve())!r}, {model.model_id!r})")
    body.extend(
        [
            "cmd.hide('everything')",
            "cmd.show('cartoon')",
            "cmd.color('purple', 'chain A')",
            "cmd.color('cyan', 'chain B')",
        ]
    )

    def add_select(name: str, attr: str, extra: list[str]) -> None:
        parts = selection_parts(models, attr)
        expr = " or ".join(parts) if parts else "none"
        body.append(f"cmd.select({name!r}, {expr!r})")
        if parts:
            body.extend(extra)

    add_select(
        "flex_full",
        "flex_full",
        ["cmd.show('sticks', 'flex_full')", "cmd.color('blue', 'flex_full')"],
    )
    add_select(
        "avoid",
        "avoid",
        ["cmd.show('sticks', 'avoid')", "cmd.color('red', 'avoid')"],
    )
    add_select(
        "hotspots",
        "hotspots",
        ["cmd.show('sticks', 'hotspots')", "cmd.color('green', 'hotspots')"],
    )
    bkbn = selection_parts(models, "flex_bkbn")
    bkbn_expr = " or ".join(bkbn) if bkbn else "none"
    body.append(f"cmd.select('flex_bkbn', {bkbn_expr!r})")
    if bkbn:
        body.append("cmd.color('blue', 'flex_bkbn and name n+ca+c+o')")

    anchor = models[0]
    for model in models[1:]:
        mobile = f"{model.model_id} and chain {model.target_chain}"
        target = f"{anchor.model_id} and chain {anchor.target_chain}"
        body.append(f"cmd.align({mobile!r}, {target!r})")
    zoom = f"{anchor.model_id} and chain {anchor.target_chain}"
    body.append(f"cmd.zoom({zoom!r})")
    body.append(f"cmd.save({str(pse_path.resolve())!r})")

    lines = [
        f"# Top 0.5% models for {run}",
        "python",
        *body,
        "python end",
        "quit",
    ]
    pml_path.parent.mkdir(parents=True, exist_ok=True)
    pml_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_pymol(pml_path: Path) -> None:
    if not PYMOL.is_file():
        raise FileNotFoundError(f"PyMOL not found at {PYMOL}")
    completed = subprocess.run(
        [str(PYMOL), "-c", "-q", str(pml_path)],
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        sys.stderr.write(completed.stdout[-4000:])
        sys.stderr.write(completed.stderr[-4000:])
        raise RuntimeError(f"PyMOL failed on {pml_path} (exit {completed.returncode})")
    pse_hint = pml_path.with_suffix(".pse")
    if not pse_hint.is_file():
        sys.stderr.write(completed.stdout[-4000:])
        sys.stderr.write(completed.stderr[-4000:])


def build_run(run: str, model_ids: list[str]) -> Path:
    run_dir = RUNS / run
    models = [build_model(run_dir, model_id) for model_id in model_ids]
    display = run_dir / "run_display"
    pml_path = display / f"{run}_top_0.5pct.pml"
    pse_path = display / f"{run}_top_0.5pct.pse"
    first = models[0]
    print(
        f"{run}: {len(models)} models, target chain {first.target_chain}, "
        f"hotspots {len(first.hotspots)}, avoid {len(first.avoid)}, "
        f"flex [] {len(first.flex_full)}, flex BKBN {len(first.flex_bkbn)}"
    )
    write_pml(run, models, pml_path, pse_path)
    run_pymol(pml_path)
    if not pse_path.is_file() or pse_path.stat().st_size == 0:
        raise RuntimeError(f"PyMOL did not write {pse_path}")
    print(f"  wrote {pse_path} ({pse_path.stat().st_size} bytes)")
    return pse_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    parser.add_argument("--run", help="Build only this run directory name")
    args = parser.parse_args()
    runs = load_models(args.csv, args.run)
    if not runs:
        raise SystemExit(f"No models found in {args.csv}")
    for run, model_ids in runs.items():
        build_run(run, model_ids)


if __name__ == "__main__":
    main()
