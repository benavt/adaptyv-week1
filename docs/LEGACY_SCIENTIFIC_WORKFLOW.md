# Project instructions

This project develops and evaluates EGFR binders. Preserve original scientific
files and unrelated work. Use stable design/run IDs and source-backed provenance.

## Mimas design and folding tasks

For navigation and new artifact placement, read
[docs/PROJECT_ORGANIZATION.md](docs/PROJECT_ORGANIZATION.md),
[01_Staging/README.md](01_Staging/README.md) and
[03_Filtering/README.md](03_Filtering/README.md). Existing scientific paths and
pinned submodule locations remain stable. Use the campaign pages and registry
indexes to find cross-stage artifacts; do not infer authority from filenames.

Use [docs/MIMAS_SKILLS.md](docs/MIMAS_SKILLS.md) to route design and folding
requests. These operational guides are pinned Git submodules, not automatically
installed Codex skills or local prediction environments.

- For binder generation/design on Mimas, read
  [Mimas Design Skills](02_Design/Mimas_Design_Skills/README.md) and its
  [remote agent guide](02_Design/Mimas_Design_Skills/remote_agent_guide.md), then
  the selected platform's input validation and CLI references. Routes include
  BoltzGen, BindCraft, BindCraft2 and Foundry RFD3 + MPNN.
- For folding/refolding on Mimas, read the structure inference workflow below,
  [Mimas Folding Skills](01_Staging/Mimas_Folding_Skills/README.md) and its
  [remote agent guide](01_Staging/Mimas_Folding_Skills/remote_agent_guide.md),
  then the chosen endpoint's native input and runtime references. Do not apply
  the single-protein smoke campaign to a multichain EGFR/binder complex.
- Initialize missing submodules at the project's pinned commits using the
  command in the routing guide. Check actual remote runtime paths before
  execution; preparation, transfer and completed prediction are distinct.
- Keep new scientific artifacts and run records in this project's run
  directories, outside the skill source trees. Record the submodule commit
  along with the actual runtime/version/checkpoint and submission receipt.

## Structure inference requests

For any request to prepare, submit, rerun, import, or analyze a structure
prediction, read [docs/STRUCTURE_INFERENCE.md](docs/STRUCTURE_INFERENCE.md) first.
Use its intake, routing, validation, execution, and result-collection workflow.
Also read [docs/DATA_CONTRACT.md](docs/DATA_CONTRACT.md) and
[docs/WORKFLOW_SPEC.md](docs/WORKFLOW_SPEC.md).

The preparation references bundled in
[Mimas Folding Skills](01_Staging/Mimas_Folding_Skills/README.md), especially
[input workflows](01_Staging/Mimas_Folding_Skills/input_workflows.md) and
[input validation](01_Staging/Mimas_Folding_Skills/input_validation.md), supply
AIS Stager preparation recipes and acceptance checks. Its `workspace/`, `src/`,
CLI commands, and historical examples belong to the separate AIS Stager
repository, not this checkout. Resolve that dependency before using its APIs;
do not assume the copied docs install code or a prediction runtime. Blueprint
features are proposals. Follow this project's paths for project artifacts.

Record organism/isoform, actual sequences, constructs and numbering, mutations,
tags/linkers, chain order/copies, engine/destination, seeds, template/MSA policy,
chemistry, and source evidence. Resolve only missing choices that change the
molecule or requested execution. Preserve an exact supplied construct. A fusion
is one chain; an assembly has separate chains. DD2 can be an extension of DD1
in one chain: do not infer topology from a design family or filename.

Inspect JSON contents to identify native AF3 versus AlphaFold Server dialects.
Expand copy counts/ID lists when comparing JSON and FASTA. Check chemistry,
dependencies, sequence agreement, and actual chain roles at each tool boundary.
Never silently drop chemistry, change a sequence, substitute a source, or treat
a syntax check as runtime validation. Never use display-remapped structures
for scoring without a recorded mapping.

Use fresh run directories under `01_Staging/`, existing engine-specific runners,
and `scripts/sync.sh` for remote transfer. Input preparation, transfer, submission,
and completed prediction are separate outcomes. If execution is requested,
continue through collection and verification when the runtime is available;
otherwise retain the prepared inputs and report the concrete missing dependency.
Do not execute the optional `sync.sh run` command without inspecting its configured
action and matching it to the requested job. Keep credentials private.

Complete a run record using [docs/RUN_RECORD_TEMPLATE.md](docs/RUN_RECORD_TEMPLATE.md).
Retain raw structures, confidence/PAE sidecars, logs, manifests, and per-job
outcomes. Refresh and validate the model registry after importing new scientific
artifacts. Report validation level and unresolved issues precisely.

<!-- fleet-week1-additive-workspace -->
## Shared campaign workspace

For a concise entry point use `python3 scripts/workspace.py history` for existing
campaigns, or `show` for future shared-workflow runs. Read
[docs/CAMPAIGN_WORKSPACE.md](docs/CAMPAIGN_WORKSPACE.md). Historical scientific
paths, IDs, registry, native runner contracts and instructions above remain
preserved. Do not treat historical manifests as complete frozen-run records.
New campaigns choose scientific methods explicitly; no TNFA defaults apply.
Shared workflow versions and separate private guide releases are pinned in
workflow.lock.json. Live submodule checkouts have not been reset or updated.
<!-- /fleet-week1-additive-workspace -->
