# Week1 — EGFR

Independent scientific project `EGFR` for Adaptyv Challenge 01. Prepare inputs
and analyze results on the Mac; run requested GPU work through Slurm on Mimas.
This is the primary source repository. Its historical predecessor is retained
as a read-only GitHub archive; historical scientific files are organized by campaign locally.

Start with `python3 scripts/workspace.py show` and [overview.md](../overview.md).
Use `history` to browse the 11 existing campaigns. Read
[workspace instructions](WORKSPACE.md) before new work.

```text
project.json                    Target and chosen pipeline responsibilities
workflow.lock.json              Deliberately pinned shared code and guides
campaigns/<campaign>/
  campaign.json                 Scientific question and recipe
  01_staging/ ... 04_analysis/   Preserved historical packages
  catalog.json / artifacts.jsonl Original paths, file counts and SHA-256 checksums
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

Week1 and Week2 use the same new-run skeleton. Week1 also has eleven imported
historical campaigns, each with the four stage folders. Shared targets, archived
helpers, engine checkouts and original registry evidence remain in `reference/`.
Use `python3 scripts/artifacts.py list` for compact campaign counts, and
`locate ORIGINAL_PATH --verify` to locate an old path and report its current checksum.
The private map also resolves retired navigation shortcuts without duplicate data.

Native execution adapters for the new run paths still require validation before
scientific compute. No inference was submitted during organization.

Git tracks reviewed source, docs, pins, campaign configuration and the explicitly published README presentation package. Large models,
structures, predictions, renderings, logs, registry observations and private
settings stay outside Git. Fleet resolves their locations, verifies selected
downloads, and copies stable remote folders to the registered archive without
deleting working sources. A Git clone alone does not restore scientific data
or private bindings. Source with embedded site details remains local pending
parameterization; see the private migration receipt rather than publishing it.

Use stable `main`, a `work/YYYY-MM-DD_task-slug` branch for each bounded task,
and a separate task worktree. Changed inputs/settings create a successor run;
update workflow pins deliberately. Preserve explicit user scientific choices.

Read [challenge](CHALLENGE.md) and [historical navigation](HISTORY.md).
