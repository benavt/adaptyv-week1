"""
Binding-site and classification logic for PyMOL coloring.

Source:
  - scripts/merge_csv.py compute_classification (lines 173-248)
  - scripts/visualize.py write_pymol_occlusion_script binding union (lines 1604-1635)
"""

from __future__ import annotations

from typing import Any, Dict, Iterable, List, Set


def to_bool(value: Any) -> bool:
    """Parse common bool-like CSV / JSON values."""
    if value is None:
        return False
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    return text in ("true", "1", "yes", "y", "t")


def is_significant(row: Dict[str, str], *, field: str = "significant") -> bool:
    """Return True when a significance column is truthy (1 / True / yes)."""
    raw = (row.get(field) or "").strip()
    if not raw:
        return False
    if raw.isdigit():
        return int(raw) == 1
    return to_bool(raw)


def is_binding_site_residue(row: Dict[str, str]) -> bool:
    """
    Ground truth for binding site: occluded OR CA-distance filter OR interaction
    OR min inter-chain atom distance < 2 A.

    Matches merge_csv.compute_classification positive-strategy and
    visualize.write_pymol_occlusion_script residue union.
    """
    if to_bool(row.get("is_occluded_occlusion", "")):
        return True
    if to_bool(row.get("passes_filter_distance", "")):
        return True
    if to_bool(row.get("has_hbond_interaction", "")):
        return True
    if to_bool(row.get("has_charge_complement_interaction", "")):
        return True
    if to_bool(row.get("has_pi_contact_interaction", "")):
        return True
    if to_bool(row.get("passes_sub_2A_filter_any_atom", "")):
        return True

    min_any = (row.get("min_any_atom_distance") or "").strip()
    if min_any:
        try:
            return float(min_any) < 2.0
        except ValueError:
            pass
    return False


def compute_classification(row: Dict[str, str], *, significant_field: str = "significant") -> str:
    """
    Return TP/FP/TN/FN for a master_alignment row, or '' when significance is missing.
    """
    sig_raw = (row.get(significant_field) or "").strip()
    if not sig_raw:
        return ""
    if sig_raw.isdigit():
        significant = int(sig_raw) == 1
    elif to_bool(sig_raw):
        significant = True
    else:
        return ""

    positive = is_binding_site_residue(row)
    if significant and positive:
        return "TP"
    if significant and not positive:
        return "FP"
    if not significant and not positive:
        return "TN"
    return "FN"


def binding_site_residue_numbers(rows: Iterable[Dict[str, str]]) -> Set[int]:
    """Collect holo_resi values that belong to the binding-site union."""
    out: Set[int] = set()
    for row in rows:
        resi = parse_holo_resi(row)
        if resi is not None and is_binding_site_residue(row):
            out.add(resi)
    return out


def significant_residue_numbers(
    rows: Iterable[Dict[str, str]],
    *,
    significant_field: str = "significant",
) -> Set[int]:
    """Collect holo_resi values flagged significant for CSP mask coloring."""
    out: Set[int] = set()
    for row in rows:
        resi = parse_holo_resi(row)
        if resi is not None and is_significant(row, field=significant_field):
            out.add(resi)
    return out


def classification_by_residue(rows: Iterable[Dict[str, str]]) -> Dict[int, str]:
    """
    Map holo_resi -> TP/FP/TN/FN using the precomputed ``classification`` column
    when present; otherwise derive from significance + binding-site flags.
    """
    out: Dict[int, str] = {}
    for row in rows:
        resi = parse_holo_resi(row)
        if resi is None:
            continue
        cls = (row.get("classification") or "").strip()
        if cls not in ("TP", "FP", "TN", "FN"):
            cls = compute_classification(row)
        if cls in ("TP", "FP", "TN", "FN"):
            out[resi] = cls
    return out


def parse_holo_resi(row: Dict[str, str]) -> int | None:
    """Parse ``holo_resi`` (PDB residue number used as PyMOL ``resi``)."""
    raw = (row.get("holo_resi") or "").strip()
    if not raw:
        return None
    try:
        return int(float(raw))
    except (ValueError, TypeError):
        return None


def infer_receptor_chain(rows: List[Dict[str, str]], default: str = "A") -> str:
    """Read receptor chain from CSV ``chain`` column when available."""
    for row in rows:
        chain = (row.get("chain") or "").strip()
        if chain:
            return chain
    return default
