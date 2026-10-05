# Week1 EGFR workspace

Eleven historical campaigns: `python3 scripts/workspace.py history`.
Query one campaign with `history CAMPAIGN` and `scripts/artifacts.py list --campaign CAMPAIGN`.
All campaign-owned historical packages are under campaigns/<id>/<stage>/.
Shared tools and evidence remain in reference/. See [history](docs/HISTORY.md).

New executions use campaigns/<id>/runs/YYYY-MM-DD_label_NNN with immutable
inputs and explicit scientific settings. Historical imports are not frozen runs.
Use `show` and one selected stage guide before new work. Native launch adapters
still require an explicit dry run and validation before scientific computation.
