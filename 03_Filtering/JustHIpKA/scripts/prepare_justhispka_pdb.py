#!/usr/bin/env python3
"""Create a conservative protein-only PDB suitable for JustHISpKa."""
from __future__ import annotations
import argparse, shutil
from pathlib import Path

STANDARD = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}

def prepare(source: Path, output: Path) -> tuple[int, int]:
    output.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, output.parent / (source.name + ".original"))
    kept = skipped = 0
    seen_alt = set()
    with source.open() as inp, output.open("w") as out:
        for line in inp:
            if not line.startswith("ATOM  "):
                continue
            name = line[17:20].strip().upper()
            if name not in STANDARD:
                skipped += 1
                continue
            alt = line[16:17]
            if alt not in (" ", "A"):
                continue
            chain = line[21:22] or " "
            if chain == " ":
                line = line[:21] + "A" + line[22:]
            key = (line[21:22], line[22:27], line[12:16].strip())
            if alt == "A":
                seen_alt.add(key)
            out.write(line.rstrip("\n") + "\n")
            kept += 1
        out.write("TER\nEND\n")
    return kept, skipped

def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("source", type=Path)
    p.add_argument("output", type=Path)
    args = p.parse_args()
    kept, skipped = prepare(args.source, args.output)
    print(f"wrote {args.output}: {kept} ATOM records; skipped {skipped} non-standard records")

if __name__ == "__main__":
    main()
