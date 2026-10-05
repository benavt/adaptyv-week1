"""
Minimal PyMOL display scripts for single-chain centroid PDB models.

Unlike the case-study writers, these do not need master_alignment.csv or
CSP/classification coloring — just a clean cartoon for interactive view capture.
"""

from __future__ import annotations

from pathlib import Path


def _pdb_load_path(pdb_path: str | Path) -> str:
    """Return a path string suitable for PyMOL ``load`` commands."""
    return str(Path(pdb_path).resolve())


def write_centroid_display_script(
    structure_pdb: str | Path,
    out_path: str | Path,
    *,
    object_name: str = "structure",
) -> Path:
    """
    Write a minimal cartoon display ``.pml`` for one centroid PDB.

    Style: gray30 cartoon, hydrogens hidden, white background.
    """
    pdb_ref = _pdb_load_path(structure_pdb)
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    lines = [
        "reinitialize",
        f"load {pdb_ref}, {object_name}",
        f"hide everything, {object_name}",
        f"show cartoon, {object_name}",
        f"color gray30, {object_name}",
        f"hide (hydro and {object_name})",
        "bg_rgb white",
        f"# Centroid display: gray30 cartoon for {Path(structure_pdb).name}",
    ]

    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out
