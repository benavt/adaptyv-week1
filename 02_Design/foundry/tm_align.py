#!/usr/bin/env python3
"""TM-score / RMSD between two identical-sequence structures via USalign."""

from __future__ import annotations

import argparse
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ALIGNED_RE = re.compile(
    r"Aligned length=\s*(\d+)\s*,\s*RMSD=\s*([0-9.]+).*?=\s*([0-9.]+)",
    re.I,
)
TM_RE = re.compile(r"TM-score=\s*([0-9.]+)", re.I)
STRUCT_EXTS = {".cif", ".pdb", ".cif.gz", ".pdb.gz"}
MACHO_MAGICS = {
    b"\xfe\xed\xfa\xce",
    b"\xfe\xed\xfa\xcf",
    b"\xce\xfa\xed\xfe",
    b"\xcf\xfa\xed\xfe",
    b"\xca\xfe\xba\xbe",
    b"\xbe\xba\xfe\xca",
}


def _is_runnable(path: Path) -> bool:
    if not path.is_file() or not os.access(path, os.X_OK):
        return False
    try:
        magic = path.read_bytes()[:4]
    except OSError:
        return False
    if magic in MACHO_MAGICS:
        return False
    return True


def resolve_usalign(explicit: str | None = None) -> str:
    candidates: list[Path] = []
    if explicit:
        candidates.append(Path(explicit).expanduser())
    env = os.environ.get("USALIGN")
    if env:
        candidates.append(Path(env).expanduser())
    which = shutil.which("USalign")
    if which:
        candidates.append(Path(which))
    here = Path(__file__).resolve().parent / "USalign"
    candidates.append(here)

    seen: set[Path] = set()
    for cand in candidates:
        try:
            resolved = cand.resolve()
        except OSError:
            continue
        if resolved in seen:
            continue
        seen.add(resolved)
        if _is_runnable(resolved):
            return str(resolved)
    raise FileNotFoundError(
        "USalign not found or not runnable on this OS. "
        "Put a Linux build on PATH, set USALIGN, or pass --usalign. "
        "The repo-local ./USalign is a macOS binary and will not run here."
    )


def _needs_decompress(path: Path) -> bool:
    name = path.name.lower()
    return name.endswith(".gz")


def _open_for_usalign(path: Path, tmpdir: str) -> str:
    path = path.resolve()
    if not path.is_file():
        raise FileNotFoundError(f"structure not found: {path}")
    if not _needs_decompress(path):
        return str(path)
    dest = Path(tmpdir) / path.with_suffix("").name
    if dest.suffix.lower() not in {".cif", ".pdb"}:
        dest = dest.with_suffix(dest.suffix + ".cif")
    with gzip.open(path, "rb") as src, dest.open("wb") as out:
        shutil.copyfileobj(src, out)
    return str(dest)


def parse_usalign_stdout(text: str) -> dict:
    aligned = ALIGNED_RE.search(text)
    tms = TM_RE.findall(text)
    if not aligned or len(tms) < 1:
        raise ValueError(f"could not parse USalign output:\n{text}")
    tm1 = float(tms[0])
    tm2 = float(tms[1]) if len(tms) > 1 else tm1
    return {
        "aligned_length": int(aligned.group(1)),
        "rmsd": float(aligned.group(2)),
        "seq_id": float(aligned.group(3)),
        "tm_score_1": tm1,
        "tm_score_2": tm2,
    }


def compute_structure_similarity(
    pdb_file1: str | os.PathLike,
    pdb_file2: str | os.PathLike,
    multimer: bool = True,
    usalign: str = "USalign",
) -> dict:
    """Run USalign and return TM-score / RMSD for two structures."""
    binary = resolve_usalign(usalign)
    file1 = Path(pdb_file1)
    file2 = Path(pdb_file2)
    with tempfile.TemporaryDirectory(prefix="tm_align_") as tmpdir:
        path1 = _open_for_usalign(file1, tmpdir)
        path2 = _open_for_usalign(file2, tmpdir)
        cmd = [binary, path1, path2, "-mol", "prot"]
        if multimer:
            cmd.extend(["-mm", "1", "-ter", "1"])
        cmd.extend(["-het", "1"])
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        detail = (result.stderr or result.stdout or "").strip()
        raise RuntimeError(
            f"USalign failed ({result.returncode}) for {file1} vs {file2}"
            + (f": {detail}" if detail else "")
        )
    parsed = parse_usalign_stdout(result.stdout)
    parsed["model1"] = str(file1)
    parsed["model2"] = str(file2)
    parsed["command"] = cmd
    return parsed


def collect_structures(root: Path) -> list[Path]:
    found: list[Path] = []
    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        name = path.name.lower()
        if any(name.endswith(ext) for ext in STRUCT_EXTS):
            found.append(path.resolve())
    return found


TSV_FIELDS = (
    "model1",
    "model2",
    "tm_score_1",
    "tm_score_2",
    "rmsd",
    "aligned_length",
    "seq_id",
)

RF3_CONFIDENCE_FIELDS = (
    "overall_plddt",
    "ptm",
    "iptm",
    "ranking_score",
    "overall_pae",
    "overall_pde",
    "has_clash",
    "chain_ptm",
)

_RF3_SKIP_KEYS = {
    "atom_chain_ids",
    "atom_plddts",
    "pae",
    "token_chain_ids",
    "token_res_ids",
}
_RF3_SKIP_PREFIXES = ("chain_pair_",)


def rf3_summary_path(model: str | os.PathLike) -> Path:
    """Sibling RF3 summary JSON for a ranked or sample *_model.cif."""
    path = Path(model)
    name = path.name
    for suffix in (".cif.gz", ".pdb.gz", ".cif", ".pdb"):
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
            break
    if name.endswith("_model"):
        name = name[: -len("_model")]
    return path.parent / f"{name}_summary_confidences.json"


def _format_confidence_value(value) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        return f"{value:.4f}"
    return str(value)


def load_rf3_confidences(model: str | os.PathLike) -> dict[str, str]:
    """Flat scalars from a sibling *_summary_confidences.json; missing file -> blanks."""
    out = {key: "" for key in RF3_CONFIDENCE_FIELDS}
    sidecar = rf3_summary_path(model)
    if not sidecar.is_file():
        return out
    try:
        data = json.loads(sidecar.read_text())
    except (OSError, json.JSONDecodeError):
        return out
    if not isinstance(data, dict):
        return out
    for key, value in data.items():
        if key in _RF3_SKIP_KEYS or key.startswith(_RF3_SKIP_PREFIXES):
            continue
        if key == "chain_ptm" and isinstance(value, list):
            if value and isinstance(value[0], list):
                continue
            out[key] = ",".join(_format_confidence_value(item) for item in value)
            continue
        if isinstance(value, (int, float, bool, str)) or value is None:
            out[key] = _format_confidence_value(value)
    return out


def output_fields(rows: list[dict]) -> tuple[str, ...]:
    extras: list[str] = []
    reserved = set(TSV_FIELDS) | set(RF3_CONFIDENCE_FIELDS) | {"command"}
    for row in rows:
        for key in row:
            if key in reserved or key in extras:
                continue
            extras.append(key)
    return TSV_FIELDS + RF3_CONFIDENCE_FIELDS + tuple(extras)


def _tsv_cell(row: dict, key: str) -> str:
    if key not in row or row[key] == "" or row[key] is None:
        return ""
    value = row[key]
    if key in ("model1", "model2", "aligned_length", "has_clash", "chain_ptm"):
        return str(value)
    if isinstance(value, float):
        return f"{value:.4f}"
    return str(value)


def print_tsv(rows: list[dict]) -> None:
    fields = output_fields(rows)
    print("\t".join(fields))
    for row in rows:
        print("\t".join(_tsv_cell(row, key) for key in fields))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="TM-score / RMSD via USalign for identical-sequence structures."
    )
    parser.add_argument("models", nargs="*", help="two structure files (cif/pdb)")
    parser.add_argument(
        "--dir",
        dest="directory",
        help="directory of structures to score against --ref (recursive)",
    )
    parser.add_argument("--ref", help="reference structure (required with --dir)")
    parser.add_argument(
        "--monomer",
        action="store_true",
        help="omit -mm 1 -ter 1 (default is multimer / complex)",
    )
    parser.add_argument(
        "--usalign",
        default=os.environ.get("USALIGN"),
        help="USalign binary (default: $USALIGN or USalign on PATH)",
    )
    args = parser.parse_args(argv)
    multimer = not args.monomer

    try:
        binary = resolve_usalign(args.usalign)
        if args.directory:
            if not args.ref:
                parser.error("--dir requires --ref")
            if args.models:
                parser.error("do not pass model files together with --dir")
            root = Path(args.directory)
            if not root.is_dir():
                raise NotADirectoryError(f"not a directory: {root}")
            ref = Path(args.ref).resolve()
            if not ref.is_file():
                raise FileNotFoundError(f"reference not found: {ref}")
            rows = []
            for model in collect_structures(root):
                if model == ref:
                    continue
                row = compute_structure_similarity(
                    model, ref, multimer=multimer, usalign=binary
                )
                row.update(load_rf3_confidences(model))
                rows.append(row)
            if not rows:
                raise FileNotFoundError(f"no cif/pdb files under {root}")
            print_tsv(rows)
            return 0

        if len(args.models) != 2:
            parser.error("provide two structure files, or --dir and --ref")
        row = compute_structure_similarity(
            args.models[0], args.models[1], multimer=multimer, usalign=binary
        )
        row.update(load_rf3_confidences(args.models[0]))
        print_tsv([row])
        return 0
    except (FileNotFoundError, NotADirectoryError, RuntimeError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
