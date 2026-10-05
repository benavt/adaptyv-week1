#!/usr/bin/env python3
"""Fail unless the four EGFR disulfides and protonated His37/His100 are present.

pdb2gmx keeps the GROMACS names CYS and HISH in the topology. The Amber
building blocks are still CYX (SG type S, no HG) and HIP (HD1 and HE2).
"""
import re
import sys
from pathlib import Path

top_path = Path(sys.argv[1] if len(sys.argv) > 1 else "topol.top")
log_path = Path(sys.argv[2] if len(sys.argv) > 2 else "pdb2gmx.out")
top = top_path.read_text().splitlines()
log = log_path.read_text() if log_path.exists() else ""

in_atoms = False
in_bonds = False
sgs = {}
atoms_by_res = {}
bonds = []
for line in top:
    if line.startswith("[ atoms ]"):
        in_atoms, in_bonds = True, False
        continue
    if line.startswith("[ bonds ]"):
        in_atoms, in_bonds = False, True
        continue
    if line.startswith("["):
        in_atoms = in_bonds = False
        continue
    if not line.strip() or line.lstrip().startswith(";"):
        continue
    parts = line.split()
    if in_atoms and len(parts) >= 5 and parts[0].isdigit():
        nr = int(parts[0])
        atype = parts[1]
        resnr = int(parts[2])
        resname = parts[3]
        atom = parts[4]
        atoms_by_res.setdefault(resnr, {"name": resname, "atoms": {}})
        atoms_by_res[resnr]["atoms"][atom] = atype
        if atom == "SG":
            sgs[nr] = resnr
    elif in_bonds and len(parts) >= 2 and parts[0].isdigit():
        bonds.append((int(parts[0]), int(parts[1])))

pairs = []
for a, b in bonds:
    if a in sgs and b in sgs:
        pairs.append(tuple(sorted((sgs[a], sgs[b]))))
expected = {(4, 29), (137, 166), (173, 182), (177, 190)}
errors = []
if set(pairs) != expected:
    errors.append("SG-SG pairs %s != %s" % (sorted(pairs), sorted(expected)))
for resnr in sorted({r for pair in expected for r in pair}):
    info = atoms_by_res.get(resnr, {"name": None, "atoms": {}})
    sg_type = info["atoms"].get("SG")
    if sg_type != "S" or "HG" in info["atoms"]:
        errors.append(
            "residue %s %s SG type %s is not a CYX disulfide"
            % (resnr, info["name"], sg_type)
        )
for resnr in (37, 100):
    info = atoms_by_res.get(resnr, {"name": None, "atoms": {}})
    names = info["atoms"]
    if info["name"] not in ("HISH", "HIP") or "HD1" not in names or "HE2" not in names:
        errors.append(
            "residue %s %s is not protonated HIP (need HISH/HIP with HD1 and HE2)"
            % (resnr, info["name"])
        )
for match in re.finditer(r"Will use HIS[ED] for residue (\d+)", log):
    if int(match.group(1)) in (37, 100):
        errors.append(match.group(0))
if errors:
    print("Protonation/disulfide check failed:", file=sys.stderr)
    for err in errors:
        print(" ", err, file=sys.stderr)
    sys.exit(1)
print("Disulfides: CYX parameters on", sorted(expected))
print("His37 and His100 are protonated (HIP parameters)")
