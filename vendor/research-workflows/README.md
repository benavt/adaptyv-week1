# Research workflows

Shared tooling and workflow contracts for research projects. The canonical
Mac checkout is `~/Code/research-workflows`; this repository is the global
coordination home for portable defaults, provenance and execution boundaries.
Each scientific project owns its target, campaigns, settings and domain rules.
Fleet is the machine inventory and SSH/Tailscale access layer: it resolves
approved hosts, bindings and access, but does not own project workflow policy.

## Project ownership and adoption

| Layer | Owns | Does not own |
| --- | --- | --- |
| `research-workflows` | Cross-project workflow contracts, portable defaults, provenance rules and thin adapter interfaces | Study-specific targets, hypotheses, thresholds or campaign history |
| Scientific project | Its campaigns, scientific settings, project adapter/overrides and pinned workflow revision | Copies of global policy or private machine credentials |
| Fleet | Machine inventory and SSH/Tailscale bindings and access | Research decisions, project lifecycle policy or scientific execution semantics |

Research repositories should pin this repository as a Git submodule at
`external/research-workflows`, using a deliberate commit update reviewed with
the project adapter. This is the adoption policy, not a claim that every
existing repository has already added the submodule. Historical data and runs
stay where they are unless their consumers and dependencies are audited.

The global contract supplies defaults; the project adapter supplies explicit
project-specific choices and may override defaults with a documented reason.
User instructions and validated scientific settings take precedence. Keep
adapter code thin: translate the shared contract to the project workflow and
resolve host/access information through Fleet. Do not duplicate global policy
in every `AGENTS.md`; use the short pointer recipe below.

For each adopting project, keep `AGENTS.md` concise and add this pointer:

```text
For shared workflow policy, start at external/research-workflows/README.md;
read only the selected contract/stage and this project's adapter/overrides.
Use Fleet only to resolve machine and SSH/Tailscale access. Do not move frozen
scientific paths or change project pins without auditing their consumers.
```

Adopt the submodule and migrate project adapters deliberately. Preserve
historical research paths and data until consumers and dependencies are audited.

A project uses `scripts/workspace.py`. Start with `show`, select one campaign and
stage, then read only that guide. `new` freezes explicit inputs/settings; `verify`
checks their integrity. Upload and selected collection preview by default;
`--apply` executes SSH/rsync transport, checks hashes and saves private receipts.
The CLI does not submit compute jobs or choose scientific settings.

One campaign describes a scientific question. One run freezes an attempt and its
parent runs. One candidate identifies an exact generated sequence and its source,
parents and measurements. A frozen selection combines candidates across campaigns
in ranked order. Read [sequence provenance](docs/provenance.md) when registering,
scoring or selecting sequences; a run ID alone does not explain a final molecule.

## Navigation

- `scripts/workflow.py`: portable Python 3.9+ project entry point.
- `stages/`: selected stage responsibilities and acceptance gates.
- `standards/`: figures, molecular views and presentations.
- `docs/provenance.md`: candidate lineage, observations and selection drafts.
- `docs/execution-contract.md`: portable execution, monitoring and evidence boundaries.
- `examples/`: configuration and provenance templates; placeholders are not evidence.
- `tests/`: integrity, pins, lineage, evidence and transport checks.

Changed inputs/settings create a new immutable date-first run, with New York
calendar dates and timezone-explicit timestamps. Parent references record origins;
selected inputs must be materialized explicitly. Progress and transport receipts
remain separate from frozen inputs. Projects keep actual locations, bulk artifacts
and private references outside their Git allowlist.

Upstream engines and operational guide repositories retain independent identities.
This release provides organization, verified transport and sequence records.
Native preparation, runtime/checkpoint identity and scientific acceptance require
separate validation before execution. Recorded generation metadata is not proof
that a job ran correctly. A selection package remains a draft until reviewed and
submitted through the competition process.

Run `python3 -m unittest discover -s tests` before committing changes. Update project
pins deliberately after validated shared changes.

For generic new runs, use [pipeline defaults](docs/defaults.md), then freeze the
effective user/recipe/adapter settings. Slurm schedules jobs submitted by the
project adapter. For remote execution, follow the
[execution contract](docs/execution-contract.md); Fleet supplies machine and
SSH/Tailscale access resolution only.
