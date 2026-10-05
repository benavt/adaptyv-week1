#!/bin/bash
# Solvate and ionize one EGFR segment. Outputs go to Human/prep or Mouse/prep.
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib.sh"
[[ $# -eq 1 ]] || usage
resolve_system "$1"
enable_log "$WORK/prep_structure.log" "$1"

echo "System: $SYS"
echo "Work directory: $WORK"
mkdir -p "$WORK"

awk 'substr($0,1,4)=="ATOM" && substr($0,22,1)=="A"' \
  "$ROOT/$SRC_PDB" > "$WORK/$PROT_PDB"
printf 'END\n' >> "$WORK/$PROT_PDB"

if grep -E 'NAG|BMA|MAN' "$WORK/$PROT_PDB"; then
  echo "Glycan records remain in $WORK/$PROT_PDB" >&2
  exit 1
fi
atoms=$(grep -c '^ATOM' "$WORK/$PROT_PDB")
if [[ "$atoms" -ne "$EXPECT_ATOMS" ]]; then
  echo "Expected $EXPECT_ATOMS ATOM records in $PROT_PDB, found $atoms" >&2
  exit 1
fi
echo "Protein PDB: $WORK/$PROT_PDB ($atoms atoms, chain A only)"

# HISH is the GROMACS name for doubly protonated histidine (Amber HIP).
# Columns 18-21 hold the 4-character residue name; column 22 stays the chain ID.
awk '
  substr($0,1,4)=="ATOM" {
    resn = substr($0,18,3)
    resi = substr($0,23,4) + 0
    if (resn == "HIS" && (resi == 37 || resi == 100)) {
      $0 = substr($0,1,17) "HISH" substr($0,22)
    }
  }
  { print }
' "$WORK/$PROT_PDB" > "$WORK/$PROT_PDB.hish"
mv "$WORK/$PROT_PDB.hish" "$WORK/$PROT_PDB"
echo "Renamed His37 and His100 to HISH in $PROT_PDB"


cp "$SCRIPT_DIR/specbond.dat" "$WORK/specbond.dat"
cp "$SCRIPT_DIR/check.com" "$SCRIPT_DIR/check2.com" "$WORK/"

load_gmx
cd "$WORK"

PDB=$PROT_PDB
FF="amber99sb-ildn"
WATER="tip3p"

set +e
$GMX pdb2gmx -f "$PDB" -o processed.gro -p topol.top -ff "$FF" -water "$WATER" -ignh > pdb2gmx.out 2>&1
pdb_status=$?
set -e
cat pdb2gmx.out
if [[ "$pdb_status" -ne 0 ]]; then
  exit "$pdb_status"
fi
python3 "$SCRIPT_DIR/check_topology.py" topol.top pdb2gmx.out

$GMX editconf -f processed.gro -o boxed.gro -c -d 1.0 -bt cubic
$GMX solvate -cp boxed.gro -o solvated.gro -p topol.top
$GMX grompp -f "$COMMON/ions.mdp" -c solvated.gro -p topol.top -o ions.tpr
echo "SOL" | $GMX genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.05 -seed 12345

echo "Box (last line of solv_ions.gro):"
tail -1 solv_ions.gro
echo "Next: bash prep/prep_equil.sh $SYS"
