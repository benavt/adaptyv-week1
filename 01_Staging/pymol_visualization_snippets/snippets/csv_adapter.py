"""
Load master_alignment.csv into structures used by PyMOL script writers.

The CSV is the primary portable input: holo_resi is treated as PyMOL resi
(already sequence-aligned in the CSP_UBQ pipeline).
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Set

from .binding_site import (
    binding_site_residue_numbers,
    classification_by_residue,
    infer_receptor_chain,
    significant_residue_numbers,
)


@dataclass
class MasterAlignmentData:
    """Parsed master_alignment.csv ready for PyMOL script generation."""

    rows: List[Dict[str, str]]
    receptor_chain: str
    ligand_chain: str
    significant_residues: Set[int] = field(default_factory=set)
    binding_site_residues: Set[int] = field(default_factory=set)
    classifications: Dict[int, str] = field(default_factory=dict)

    @property
    def holo_pdb(self) -> str:
        for row in self.rows:
            value = (row.get("holo_pdb") or "").strip()
            if value:
                return value
        return ""


def load_master_alignment(
    csv_path: str | Path,
    *,
    receptor_chain: str | None = None,
    ligand_chain: str | None = None,
    significant_field: str = "significant",
) -> MasterAlignmentData:
    """
    Read master_alignment.csv and derive residue sets for the three case-study panels.

    Args:
        csv_path: Path to master_alignment.csv
        receptor_chain: Override chain ID for protein (default: CSV ``chain`` or A)
        ligand_chain: Override chain ID for peptide/ligand (default: B)
        significant_field: Column for CSP mask significance (default: significant)
    """
    path = Path(csv_path)
    if not path.exists():
        raise FileNotFoundError(f"master_alignment CSV not found: {path}")

    with path.open("r", newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    if not rows:
        raise ValueError(f"No data rows in {path}")

    rec_chain = receptor_chain or infer_receptor_chain(rows, default="A")
    lig_chain = ligand_chain or "B"

    return MasterAlignmentData(
        rows=rows,
        receptor_chain=rec_chain,
        ligand_chain=lig_chain,
        significant_residues=significant_residue_numbers(rows, significant_field=significant_field),
        binding_site_residues=binding_site_residue_numbers(rows),
        classifications=classification_by_residue(rows),
    )
