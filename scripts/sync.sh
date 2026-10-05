#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
export PROJECT_ROOT
# Help remains available before local setup. No config is loaded for help.
if [[ $# == 0 || ${1:-} == -h || ${1:-} == --help ]]; then
  exec python3 "$PROJECT_ROOT/scripts/remote_sync.py" --help
fi
ENV_FILE=${SYNC_ENV_FILE:-"$PROJECT_ROOT/.env.local"}
if [[ ! -f "$ENV_FILE" ]]; then
  echo 'Missing local config. Copy .env.example to .env.local and fill it in.' >&2
  exit 2
fi
set -a
# shellcheck disable=SC1090
source "$ENV_FILE"
set +a
exec python3 "$PROJECT_ROOT/scripts/remote_sync.py" "$@"
