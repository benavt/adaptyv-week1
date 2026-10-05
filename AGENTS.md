# Week1 EGFR agent entry point

Start with `python3 scripts/workspace.py show` and overview.md. Select one
campaign/stage and load only its pinned guide with `show --stage STAGE`.
Read docs/CAMPAIGN_WORKSPACE.md. New runs use the same skeleton as Week2.
Use `history [CAMPAIGN]` for historical work; keep its files and IDs intact.
Historical scientific requirements are preserved in
docs/LEGACY_SCIENTIFIC_WORKFLOW.md; read it when invoking legacy workflows.

Prepare inputs locally, submit requested GPU jobs through sbatch on Mimas,
and score/analyze selected verified results on the Mac. Slurm owns scheduling.
Native adapters for the new run paths must be validated before scientific jobs.
Record exact constructs, numbering, chemistry, chain roles, species, engine,
settings, seed and checkpoints. Do not apply TNFA assumptions to EGFR.
Explicit user settings override defaults; generic Foundry requests use 1,000
RFD3 designs, with required molecular inputs and runtime still validated.

Use date-first YYYY-MM-DD run IDs with New York dates and offset/Z timestamps.
Changed inputs/settings create immutable successor snapshots and parent links.
Keep shared code/guide pins deliberate. Candidate records, source-linked scores
and frozen selections must trace exact sequences to their generation runs.
Do not invent missing historical lineage or overwrite prior evidence.

Fleet resolves private paths, machines, environments, selected downloads and
archives. Keep bulk artifacts outside Git; never delete unique/unverified data.
For one task use a work/YYYY-MM-DD_task-slug branch and its separate worktree.
main is the stable control checkout; keep commits small and purpose-specific.
Read Fleet docs/repository-conventions.md for adoption and storage rules.
Ask a targeted question when a missing scientific choice or disruptive move
needs user input. Preserve current historical paths until their consumers are audited.
