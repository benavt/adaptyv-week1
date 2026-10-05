# Project organization and path compatibility

## Entry points

- [Staging](../01_Staging/README.md): reusable references, input preparation and handoffs.
- [Filtering](../03_Filtering/README.md): scores, selections, rendering and reports.
- [Campaigns](../03_Filtering/campaigns/README.md): navigation across stages and tools.
- [Registry](../registry/README.md): exact design/model identity and sourced statistics.
- [Software](../tools/README.md): installed dependencies and pinned source checkouts.
- [Raw results](../results/README.md): native remote downloads and historical exceptions.

Existing scientific paths, run/design IDs, pinned submodule locations and runtime
contracts are preserved. The campaign pages are an additive navigation layer;
they link existing artifacts without copying coordinates or creating symlink
trees that could be discovered as duplicate scientific inputs.

## Placement of new work

| Content | Preferred location | Compatibility rule |
| --- | --- | --- |
| New reusable reference collections | `01_Staging/references/<reference_id>/` | Existing EGFR sources remain in `01_Staging/EGFR/` |
| New input preparation | `01_Staging/runs/<run_id>/` | Use only with runners supporting that input parent; existing runners retain `01_Staging/<run_id>/` |
| Design generation | Existing `02_Design/<platform>/runs/<run_id>/` contract | Do not change the Foundry contract |
| Remote predictions and logs | Existing sync destinations in `results/` and `logs/` | Preserve native route and run-relative paths; historical staging/refolding packages stay in place |
| New campaign-specific analyses | `03_Filtering/campaigns/<campaign_id>/analyses/<method>/<analysis_id>/` | Existing tool/campaign output paths remain authoritative |
| New selection records | Campaign `selections/<selection_id>/` | Record thresholds, version, retained/rejected counts and input hashes |
| New display assets | Campaign `renders/<render_id>/` | Preserve camera and transformation mappings; scoring uses native coordinates |
| New reports | Campaign `reports/<report_id>/` | Record exact source tables/renders; new revisions get distinct names |
| Reusable code | `04_Analysis/` | Existing scripts remain callable at their original paths |
| New software / environment / cache | `tools/software/`, `tools/environments/`, `tools/cache/` | Existing environments and pinned submodules are retained; environments must be recreated rather than moved |

Create only needed folders. A scorer's `raw/` contains original scorer responses
and calculation logs; the source prediction package remains in its actual raw
result home. Source snapshots in `provenance/` preserve evidence and do not
redefine the authoritative original file.

For a new folding run, choose the supported layout before execution and use a
private per-run sync config. Set `LOCAL_INPUT_ROOT=01_Staging/runs` for new
nested runs or `LOCAL_INPUT_ROOT=01_Staging` for legacy root-level runs. Pass the
unchanged single run ID to `scripts/sync.sh` and select the config with
`SYNC_ENV_FILE`. Do not edit shared private configuration to make a navigation
change. The wrapper's generic default remains `inputs/` for compatibility.

## Campaign navigation source

Each campaign has a curated `manifest.json` with schema version 1 and these fields:

- `campaign_id`: stable navigation identity, matching its folder. This groups related scientific runs; it does not replace their exact run/design/model IDs.
- `title`, `summary`, `notes`: human-readable scope, caveats and historical exceptions.
- `locations`: entries with `role`, `label`, project-relative `path` and `evidence_path`.
- `current_outputs`: available output entries with `label`, `path` and authority `evidence_path`. This is an explicit selection, not a claim that every artifact is final.
- `documented_unavailable` (optional): `expected_path` and `evidence_path` for outputs named in documentation but absent locally. These are not available outputs.
- `runs` (optional): exact `run_id`, `phase`, `path`, `status_evidence`, optional JSON `status_locator`, optional `recorded_status` for textual evidence, optional `output_path` and explicit `superseded_by` / `supersession_evidence`.

All artifact/evidence paths resolve from the project root, including status
sources. A JSON status pointer reads the source value during rebuild. A pointer
to per-job `outcomes` records their shared status only if all agree; different
statuses produce `mixed_outcomes`, and an empty list leaves status `unknown`. Textual
status labels are curated claims backed by the linked document. Unknown status
stays `unknown`; absence of evidence does not imply failure or completion.
Supersession is recorded only when explicitly named in evidence. A larger
attempt number does not automatically supersede an earlier successful result.

Preparation, transfer, accepted submission, collection and verification are
distinct phases. Counts are copied only when explicitly recorded by the status
source; blank is unrecorded and zero is an explicit zero. Status is a local
observation, not live remote monitoring. This task did not contact a cluster.

Regenerate navigation after reviewing source evidence:

```bash
python3 scripts/build_project_navigation.py
python3 scripts/build_project_navigation.py --check
```

The builder writes campaign READMEs, the campaign landing page,
`registry/campaigns.csv`, `registry/campaign_runs.csv` and
`registry/navigation_manifest.json`. The latter retains evidence hashes. The
read-only check validates paths, local documentation links, explicit successor
references and generated freshness. It does not run predictions or recalculate
scientific scores. The existing model builder preserves these navigation indexes.

For imported scientific artifacts, separately refresh and validate the model
registry using the [workflow](WORKFLOW_SPEC.md). Navigation-only edits do not
require re-indexing every coordinate file.

New analysis tables still need explicit schema/source adapters in the scientific
registry. A campaign navigation link makes a table discoverable; it does not
automatically import its values as model statistics. Existing metric/scorer
adapters retain their original paths and schemas.

The [compatibility audit](../registry/navigation_compatibility.json) records the
organization change's original-path preservation, model-ID comparison and link
checks. It does not revalidate scientific predictions or remote runtime state.

## Future physical migration

No existing scientific file was moved as part of this organization. Before any
future move, inventory consumers across scripts, CSV/JSON provenance, PML/PSE,
deck source metadata, docs, runtime configurations, sync roots and registry
annotations. Embedded PSE/deck references require explicit inspection; a text
search cannot prove they are safe.

Model IDs use original project-relative source paths, so moving a source needs
an explicit registry identity migration, including annotations and parent links.
Retain source hashes, a before/after path map and old/new model IDs. Update
consumers and discovery rules, refresh/validate the registry, and verify affected
runner/report entry points before removing old paths. A convenience symlink is
not a complete migration: the uploader rejects symlinks and the scientific
registry deliberately excludes symlink inputs.
