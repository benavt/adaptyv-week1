#!/bin/bash
# Energy minimization, 500 ps NVT, 500 ps NPT. Run from EGFR_D3 after prep_structure.sh.
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib.sh"
[[ $# -eq 1 ]] || usage
resolve_system "$1"
enable_log "$WORK/prep_equil.log" "$1"

if [[ ! -f "$WORK/solv_ions.gro" || ! -f "$WORK/topol.top" ]]; then
  echo "Missing $WORK/solv_ions.gro or topol.top. Run: bash prep/prep_structure.sh $SYS" >&2
  exit 1
fi

echo "System: $SYS"
echo "Work directory: $WORK"
load_gmx
cd "$WORK"

# Physical GPU 3. After this assignment GROMACS sees that device as GPU 0.
export CUDA_VISIBLE_DEVICES=3
echo "Using physical GPU ${CUDA_VISIBLE_DEVICES}"

$GMX grompp -f "$COMMON/minim.mdp" -c solv_ions.gro -p topol.top -o em.tpr
$GMX mdrun -deffnm em -nb gpu -pme gpu -v

$GMX grompp -f "$COMMON/nvt.mdp" -c em.gro -r em.gro -p topol.top -o nvt.tpr
$GMX mdrun -deffnm nvt -nb gpu -pme gpu -v

$GMX grompp -f "$COMMON/npt.mdp" -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
$GMX mdrun -deffnm npt -nb gpu -pme gpu -v

cp npt.gro equil.gro

SUMMARY_FILE="system_summary.txt"
{
  echo "--------------------------------------------------"
  echo "System summary generated on: $(date)"
  echo "--------------------------------------------------"
  echo "Directory: $(pwd)"
} >> "$SUMMARY_FILE"
ATOM_COUNT=$(grep -v '^;' solv_ions.gro | grep -v '^$' | tail -n +3 | head -n -1 | wc -l)
WATER_COUNT=$(grep -c "SOL" solv_ions.gro)
ION_COUNT=$(grep -cE "NA|CL" solv_ions.gro)
printf "Atoms: %d, Waters: %d, Ions: %d\n\n" \
  "$ATOM_COUNT" "$WATER_COUNT" "$ION_COUNT" >> "$SUMMARY_FILE"
echo "System composition written to $WORK/$SUMMARY_FILE"
echo "Next: bash prep/setup_md_10ns.sh $SYS"
