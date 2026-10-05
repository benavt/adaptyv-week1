"""
Generate the three PyMOL scripts used by case-study asset rendering.

Source: scripts/visualize.py
  - write_pymol_color_csp_mask_script (700-800)
  - write_pymol_occlusion_script (1521-1662)
  - write_pymol_csp_classification_script (2063-2273)

These writers accept explicit PDB paths and precomputed residue sets (from CSV),
avoiding CSPResult, sequence alignment, and PDB B-factor rewriting.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, Set

from .colors import CLASSIFICATION_COLORS, hex_to_rgb01
from .csv_adapter import MasterAlignmentData


def _pdb_load_path(pdb_path: str | Path) -> str:
    """Return a path string suitable for PyMOL ``load`` commands."""
    return str(Path(pdb_path).resolve())


def write_color_csp_mask_script(
    data: MasterAlignmentData,
    structure_pdb: str | Path,
    out_path: str | Path,
    *,
    object_name: str = "structure",
) -> Path:
    """
    CSP mask panel: receptor gray30, significant CSP residues red, ligand cyan.

    Equivalent to pipeline ``color_csp_mask.pml``.
    """
    pdb_ref = _pdb_load_path(structure_pdb)
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    lines = [
        "reinitialize",
        f"load {pdb_ref}, {object_name}",
        f"hide everything, {object_name}",
        f"show cartoon, {object_name}",
        f"color cyan, {object_name} and chain {data.ligand_chain}",
        f"color gray30, {object_name} and chain {data.receptor_chain}",
    ]
    for res_num in sorted(data.significant_residues):
        lines.append(
            f"color red, {object_name} and chain {data.receptor_chain} and resi {res_num}"
        )
    lines.append(f"set cartoon_transparency, 0.2, {object_name}")
    lines.append("# CSP mask: cyan=ligand, gray30=non-significant receptor, red=significant CSP")
    lines.append(f"# Significant residues: {len(data.significant_residues)}")

    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out


def write_occlusion_script(
    data: MasterAlignmentData,
    occlusion_pdb: str | Path,
    out_path: str | Path,
    *,
    object_name: str = "structure",
) -> Path:
    """
    Binding-site panel: receptor gray30, binding-site union red, ligand cyan.

    Equivalent to pipeline ``color_occlusion.pml``.
    """
    pdb_ref = _pdb_load_path(occlusion_pdb)
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    lines = [
        "reinitialize",
        f"load {pdb_ref}, {object_name}",
        f"hide everything, {object_name}",
        f"show cartoon, {object_name}",
        f"color gray30, {object_name} and chain {data.receptor_chain}",
    ]
    for res_num in sorted(data.binding_site_residues):
        lines.append(
            f"color red, {object_name} and chain {data.receptor_chain} and resi {res_num}"
        )
    lines.append(f"color cyan, {object_name} and chain {data.ligand_chain}")
    lines.append(f"set cartoon_transparency, 0.2, {object_name}")
    lines.append("# Binding site: occluded OR CA filter OR interaction OR sub-2A contact")
    lines.append(f"# Binding site residues: {len(data.binding_site_residues)}")

    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out


def write_classification_script(
    data: MasterAlignmentData,
    structure_pdb: str | Path,
    out_path: str | Path,
    *,
    object_name: str = "structure",
) -> Path:
    """
    TP/FP/TN/FN classification panel with custom PyMOL colors.

    Equivalent to pipeline ``csp_classification_original.pml``.
    """
    pdb_ref = _pdb_load_path(structure_pdb)
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    lines = [
        "reinitialize",
        f"load {pdb_ref}, {object_name}",
        f"hide everything, {object_name}",
        f"show cartoon, {object_name}",
        f"color cyan, {object_name} and chain {data.ligand_chain}",
        f"set cartoon_tube_radius, 0.45, {object_name} and chain {data.ligand_chain}",
        f"color gray, {object_name} and chain {data.receptor_chain}",
    ]

    for label, color_name in (
        ("TP", "tp_color"),
        ("TN", "tn_color"),
        ("FP", "fp_color"),
        ("FN", "fn_color"),
    ):
        rgb = hex_to_rgb01(CLASSIFICATION_COLORS[label])
        lines.append(
            f"set_color {color_name}, [{rgb[0]:.4f}, {rgb[1]:.4f}, {rgb[2]:.4f}]"
        )

    grouped = _group_by_classification(data.classifications)
    for cls in ("TP", "FP", "TN", "FN"):
        residues = grouped.get(cls, [])
        if not residues:
            continue
        lines.append(f"# {cls}: {len(residues)} residues")
        color_name = f"{cls.lower()}_color"
        for res_num in sorted(residues):
            lines.append(
                f"color {color_name}, {object_name} and chain {data.receptor_chain} and resi {res_num}"
            )

    lines.extend(
        [
            "set cartoon_transparency, 0.2, structure",
            "set cartoon_fancy_helices, 1",
            "set cartoon_ring_mode, 1",
            "# Classification colors match pipeline bar-plot legend",
            f"# TP ({CLASSIFICATION_COLORS['TP']}): Sig. CSP in binding site",
            f"# FP ({CLASSIFICATION_COLORS['FP']}): Sig. CSP -- allosteric",
            f"# TN ({CLASSIFICATION_COLORS['TN']}): low CSP -- allosteric",
            f"# FN ({CLASSIFICATION_COLORS['FN']}): low CSP in binding site",
        ]
    )

    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out


def write_all_case_study_scripts(
    data: MasterAlignmentData,
    *,
    structure_pdb: str | Path,
    occlusion_pdb: str | Path,
    output_dir: str | Path,
) -> Dict[str, Path]:
    """
    Write all three .pml scripts into ``output_dir``.

    Returns dict mapping logical name -> path.
    """
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    return {
        "csp_mask": write_color_csp_mask_script(
            data, structure_pdb, out_dir / "color_csp_mask.pml"
        ),
        "occlusion": write_occlusion_script(data, occlusion_pdb, out_dir / "color_occlusion.pml"),
        "classification": write_classification_script(
            data, structure_pdb, out_dir / "csp_classification_original.pml"
        ),
    }


def _group_by_classification(classifications: Dict[int, str]) -> Dict[str, list[int]]:
    grouped: Dict[str, list[int]] = {"TP": [], "FP": [], "TN": [], "FN": []}
    for res_num, cls in classifications.items():
        if cls in grouped:
            grouped[cls].append(res_num)
    return grouped
