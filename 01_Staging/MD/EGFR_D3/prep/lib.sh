#!/bin/bash
# Shared paths for the EGFR D3 10 ns prep scripts.
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
COMMON="$SCRIPT_DIR/common"
GMXRC=/gpfs/u/barn/VCBM/shared/local/gromacs-2025/bin/GMXRC

usage() {
  echo "Usage: bash prep/$(basename "$0") Human|Mouse" >&2
  exit 1
}

resolve_system() {
  local arg
  arg=$(printf '%s' "${1:-}" | tr '[:upper:]' '[:lower:]')
  case "$arg" in
    human)
      SYS=Human
      SRC_PDB=Human_EGFR_seg_medoid.pdb
      PROT_PDB=Human_EGFR_seg_protein.pdb
      JOB=Human-EGFR-10ns
      EXPECT_ATOMS=1480
      ;;
    mouse)
      SYS=Mouse
      SRC_PDB=Mouse_EGFR_seg_medoid.pdb
      PROT_PDB=Mouse_EGFR_seg_protein.pdb
      JOB=Mouse-EGFR-10ns
      EXPECT_ATOMS=1503
      ;;
    *)
      usage
      ;;
  esac
  WORK="$ROOT/$SYS/prep"
  OUT="$ROOT/$SYS/md_10ns"
}

load_gmx() {
  if ! command -v gmx &> /dev/null; then
    # shellcheck disable=SC1090
    # GMXRC references unset shell and GMXLDLIB. nounset would abort the source.
    set +u
    source "$GMXRC"
    set -u
    GMX=gmx_mpi
  else
    GMX=$(command -v gmx)
  fi
  export GMXLIB="/gpfs/u/barn/VCBM/shared/local/gromacs-2025/share/gromacs/top"
  echo "Using GROMACS at: $GMX"
}

# Re-run this script once the work/output directory exists, saving a log there.
enable_log() {
  local logfile=$1
  shift
  if [[ -z "${EGFR_PREP_LOGGING:-}" ]]; then
    mkdir -p "$(dirname "$logfile")"
    export EGFR_PREP_LOGGING=1
    set +e
    bash "$0" "$@" 2>&1 | tee "$logfile"
    local status=${PIPESTATUS[0]}
    set -e
    exit "$status"
  fi
}
