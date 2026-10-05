#!/bin/bash
# Build md_10ns.tpr and the SLURM script in Human/md_10ns or Mouse/md_10ns.
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib.sh"
[[ $# -eq 1 ]] || usage
resolve_system "$1"
enable_log "$OUT/setup_md_10ns.log" "$1"

for need in npt.gro npt.cpt topol.top; do
  if [[ ! -f "$WORK/$need" ]]; then
    echo "Missing $WORK/$need. Run prep_structure.sh and prep_equil.sh for $SYS first." >&2
    exit 1
  fi
done

echo "System: $SYS"
echo "Output directory: $OUT"
mkdir -p "$OUT"
cp "$WORK/npt.gro" "$WORK/npt.cpt" "$WORK/topol.top" "$OUT/"
if [[ -f "$WORK/posre.itp" ]]; then
  cp "$WORK/posre.itp" "$OUT/"
fi
cp "$COMMON/md_10ns.mdp" "$OUT/"

load_gmx
cd "$OUT"
$GMX grompp -f md_10ns.mdp -c npt.gro -t npt.cpt -p topol.top -o md_10ns.tpr

cat > run_md_10ns.slurm << EOF
#!/bin/bash
#SBATCH --job-name=${JOB}
#SBATCH --output=md_10ns.out
#SBATCH --error=md_10ns.err
#SBATCH --partition=debug
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=32G

cd "\${SLURM_SUBMIT_DIR}"

source /gpfs/u/barn/VCBM/shared/local/gromacs-2025/bin/GMXRC
GMX=gmx_mpi

# Physical GPU 3. GROMACS sees it as device 0.
export CUDA_VISIBLE_DEVICES=3
echo "Using physical GPU ${CUDA_VISIBLE_DEVICES}"

\$GMX mdrun -deffnm md_10ns -ntomp \$SLURM_CPUS_PER_TASK -nb gpu -pme gpu -v
EOF
chmod 755 run_md_10ns.slurm

echo "tpr: $OUT/md_10ns.tpr"
echo "Submit from EGFR_D3 with:"
echo "  cd $SYS/md_10ns && sbatch run_md_10ns.slurm"
echo "The job is pinned to physical GPU 3."
