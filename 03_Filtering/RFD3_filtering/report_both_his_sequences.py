#!/usr/bin/env python3
"""Report DD1 binder sequences whose ASPs are near both EGFR HIS37 and HIS100."""
import csv
import gzip
from pathlib import Path
from Bio.PDB import MMCIFParser

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
INFILE = HERE / "out/dd1_asp_his_proximity.csv"
OUTFILE = HERE / "out/dd1_asp_near_both_his_sequences.csv"


def sequence(path):
    structure = MMCIFParser(QUIET=True).get_structure(path.stem, gzip.open(path, "rt"))
    chains = list(structure[0].get_chains())
    receptor = max(chains, key=lambda c: len(list(c.get_residues())))
    binder = next(c for c in chains if c.id != receptor.id)
    aa = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
    return "".join(aa.get(r.resname, "X") for r in binder if r.id[0] == " ")


def main():
    rows = []
    with INFILE.open() as handle:
        for row in csv.DictReader(handle):
            if row["has_asp_near_both"] != "True":
                continue
            path = Path(row["structure_path"])
            row["binder_sequence"] = sequence(path)
            rows.append({key: row[key] for key in ["run", "design_id", "binder_sequence", "asp_near_his37", "asp_near_his100", "structure_path"]})
    fields = ["run", "design_id", "binder_sequence", "asp_near_his37", "asp_near_his100", "structure_path"]
    with OUTFILE.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} sequences to {OUTFILE}")


if __name__ == "__main__":
    main()
