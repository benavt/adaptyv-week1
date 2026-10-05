# Week1 — EGFR

Independent scientific project `EGFR` for Adaptyv Challenge 01. Prepare inputs
and analyze results on the Mac; run requested GPU work through Slurm on Mimas.
This is the primary source repository. Its historical predecessor is retained
as a read-only GitHub archive; original scientific files remain in place locally.

Start with `python3 scripts/workspace.py show` and [overview.md](overview.md).
Use `history` to browse the 11 existing campaigns. Read
[workspace instructions](docs/WORKSPACE.md) before new work.

```text
project.json                    Target and chosen pipeline responsibilities
workflow.lock.json              Deliberately pinned shared code and guides
campaigns/<campaign>/
  campaign.json                 Scientific question and recipe
  runs/<YYYY-MM-DD_label_NNN>/
    run.json                    Immutable input/settings identity and parent link
    01_staging/                 Frozen native inputs
    02_design/                  Design results
    03_filtering/               Predictions, scores and filtering evidence
    04_analysis/                Comparisons, plots and reports
    state.json                  Recorded progress
    receipts/                   Verified transport and execution observations
registry/                       Candidate identity and source-backed evidence
drafts/                         Editable inputs
submissions/                    Frozen reviewed selection packages
local/                          Private bindings and migration receipts
```

The same new-run skeleton is used by Week2. The four historical stage folders
remain unchanged for existing EGFR campaigns; they are compatibility locations,
not a second system for new runs. Existing IDs and evidence are preserved.
Native execution adapters for the new run paths still require validation before
scientific compute. No inference was submitted during organization.

Git tracks reviewed source, docs, pins and campaign configuration. Large models,
structures, predictions, renderings, logs, registry observations and private
settings stay outside Git. Fleet resolves their locations, verifies selected
downloads, and copies stable remote folders to the registered archive without
deleting working sources. A Git clone alone does not restore scientific data
or private bindings. Source with embedded site details remains local pending
parameterization; see the private migration receipt rather than publishing it.

Use stable `main`, a `work/YYYY-MM-DD_task-slug` branch for each bounded task,
and a separate task worktree. Changed inputs/settings create a successor run;
update workflow pins deliberately. Preserve explicit user scientific choices.

Read [challenge](docs/CHALLENGE.md) and [historical navigation](docs/HISTORY.md).
