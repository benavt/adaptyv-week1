#!/usr/bin/env python3
"""PyMOL sessions for the five most isoform-specific binders.

Each session loads the Boltz-2 complex of those binders against the isoform
they were ranked for. Chain A (ET domain) is green and chain B (binder) is
purple. The five receptors are aligned.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
LIGAND_AI = REPO / "design" / "ligand_ai"
PYMOL = Path("/opt/homebrew/bin/pymol")
TOP_N = 5

if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from plot_ligandai_cross_isoform_kd import (  # noqa: E402
    INPUT_PATH,
    ISOFORMS,
    OUT,
    load_rows,
    top_specificity,
)

SESSION_DIR = OUT / "top5_specificity_pymol"
PDB_DIR = SESSION_DIR / "pdbs"
SOURCES_PATH = SESSION_DIR / "sources.txt"
SEQUENCES = REPO / "staging" / "sequences"
AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M",
}
ET_LENGTH_LIMIT = 120


def job_column(prefix: str) -> str:
    return "%s_job_id" % prefix


def rows_by_sequence(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    chosen: dict[str, dict[str, str]] = {}
    for row in rows:
        sequence = row["sequence"].strip()
        if sequence in chosen:
            raise SystemExit("Duplicate sequence in %s: %s" % (INPUT_PATH, sequence))
        chosen[sequence] = row
    return chosen


def source_phrase(row: dict[str, str]) -> str:
    run = row.get("run", "").strip()
    model_id = row.get("model_id", "").strip()
    gene = row.get("gene", "").strip()
    if "msa_hotspot" in run.lower():
        return (
            "ligand_ai, BRD2ET MSA-hotspot workspace (`%s`, model `%s`), "
            "not a standard isoform campaign" % (run, model_id)
        )
    return "ligand_ai, designed in %s as %s" % (gene, model_id)


def source_lines(
    ranked: dict[str, list[dict[str, object]]],
    by_sequence: dict[str, dict[str, str]],
) -> list[str]:
    lines: list[str] = []
    for _gene, _prefix, label in ISOFORMS:
        lines.append(label)
        lines.append("")
        for item in ranked[label]:
            if int(item["rank"]) > TOP_N:
                continue
            sequence = str(item["sequence"])
            phrase = source_phrase(by_sequence[sequence])
            lines.append("- %s %s — %s" % (item["rank"], sequence, phrase))
        lines.append("")
    return lines


def ligand_client():
    sys.path.insert(0, str(LIGAND_AI))
    if not (LIGAND_AI / "src" / "ligandai_local").is_dir():
        raise SystemExit("LigandAI package not found at %s" % LIGAND_AI)
    from src.ligandai_local.client import LigandAIClient

    return LigandAIClient()


def download_pdbs(
    client,
    ranked: dict[str, list[dict[str, object]]],
    by_sequence: dict[str, dict[str, str]],
) -> dict[str, list[Path]]:
    written: dict[str, list[Path]] = {}
    for _gene, prefix, label in ISOFORMS:
        column = job_column(prefix)
        paths: list[Path] = []
        for item in ranked[label]:
            rank = int(item["rank"])
            if rank > TOP_N:
                continue
            sequence = str(item["sequence"])
            row = by_sequence[sequence]
            job_id = (row.get(column) or "").strip()
            if not job_id:
                raise SystemExit("%s rank %s is missing %s" % (label, rank, column))
            path = PDB_DIR / label / ("rank%d.pdb" % rank)
            path.parent.mkdir(parents=True, exist_ok=True)
            if path.is_file() and path.stat().st_size > 100 and "ATOM" in path.read_text(encoding="utf-8", errors="replace"):
                print("  cached %s" % path.relative_to(OUT), flush=True)
            else:
                text = client.get_folding_job_structure(job_id, "pdb")
                if not isinstance(text, str) or "ATOM" not in text:
                    raise RuntimeError("%s returned no pdb structure" % job_id)
                path.write_text(text if text.endswith("\n") else text + "\n", encoding="utf-8")
                print("  wrote %s (%s bytes)" % (path.relative_to(OUT), path.stat().st_size), flush=True)
            paths.append(path)
        written[label] = paths
    return written


def fasta_sequence(label: str) -> str:
    path = SEQUENCES / ("%s.fasta" % label)
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    return "".join(line for line in lines if not line.startswith(">"))


def chain_a_sequence(pdb_path: Path) -> dict[int, str]:
    residues: dict[int, str] = {}
    for line in pdb_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("ATOM") or line[12:16].strip() != "CA" or line[21] != "A":
            continue
        residues[int(line[22:26])] = AA3.get(line[17:20].strip(), "X")
    return residues


def et_residue_span(pdb_path: Path, label: str) -> tuple[int, int] | None:
    """ET residues inside a long chain A, or None when chain A is already the domain."""
    residues = chain_a_sequence(pdb_path)
    if len(residues) <= ET_LENGTH_LIMIT:
        return None
    sequence = "".join(residues[index] for index in sorted(residues))
    offset = sequence.find(fasta_sequence(label))
    if offset < 0:
        raise SystemExit("%s chain A is %s residues and does not contain the ET sequence" % (pdb_path, len(residues)))
    numbers = sorted(residues)
    start = numbers[offset]
    end = numbers[offset + len(fasta_sequence(label)) - 1]
    print("  %s chain A is %s residues; ET is %s-%s" % (pdb_path.name, len(residues), start, end), flush=True)
    return start, end


def write_session(label: str, pdbs: list[Path]) -> Path:
    pse_path = SESSION_DIR / ("%s_top5.pse" % label)
    pml_path = SESSION_DIR / ("%s_top5.pml" % label)
    lines = [
        "from pymol import cmd",
        "cmd.reinitialize()",
    ]
    for index, pdb in enumerate(pdbs, start=1):
        obj = "rank%d" % index
        lines.append("cmd.load(%r, %r)" % (str(pdb.resolve()), obj))
        lines.append("cmd.hide('everything', %r)" % obj)
        lines.append("cmd.show('cartoon', %r)" % obj)
        span = et_residue_span(pdb, label)
        if span is None:
            lines.append("cmd.color('green', %r)" % ("%s and chain A" % obj))
        else:
            start, end = span
            lines.append("cmd.color('gray70', %r)" % ("%s and chain A" % obj))
            lines.append(
                "cmd.color('green', %r)"
                % ("%s and chain A and resi %s-%s" % (obj, start, end))
            )
        lines.append("cmd.color('purple', %r)" % ("%s and chain B" % obj))
        if index > 1:
            lines.append(
                "cmd.align(%r, %r)"
                % ("%s and chain A" % obj, "rank1 and chain A")
            )
    lines.append("cmd.orient('rank1 and chain A')")
    lines.append("cmd.zoom('rank1', 8)")
    lines.append("cmd.save(%r)" % str(pse_path.resolve()))
    pml_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    if not PYMOL.is_file():
        raise SystemExit("PyMOL not found: %s" % PYMOL)
    completed = subprocess.run(
        [str(PYMOL), "-c", "-q", str(pml_path)],
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0 or not pse_path.is_file():
        detail = (completed.stderr or completed.stdout or "").strip()
        raise SystemExit("PyMOL failed for %s:\n%s" % (label, detail))
    print("  saved %s" % pse_path.relative_to(OUT), flush=True)
    return pse_path


def main() -> None:
    rows = load_rows(INPUT_PATH)
    by_sequence = rows_by_sequence(rows)
    ranked = top_specificity(rows)
    SESSION_DIR.mkdir(parents=True, exist_ok=True)
    print("Downloading target-isoform folds", flush=True)
    client = ligand_client()
    pdbs = download_pdbs(client, ranked, by_sequence)
    print("Writing PyMOL sessions", flush=True)
    for _gene, _prefix, label in ISOFORMS:
        write_session(label, pdbs[label])
    text = "\n".join(source_lines(ranked, by_sequence)).rstrip() + "\n"
    SOURCES_PATH.write_text(text, encoding="utf-8")
    print(text, end="")
    print("Saved sources to %s" % SOURCES_PATH)


if __name__ == "__main__":
    main()
