#!/bin/zsh
set -euo pipefail

cd "$(dirname "$0")/../../../design/ligand_ai"
python -m src.ligandai_local.cli deltaforge "$@"
