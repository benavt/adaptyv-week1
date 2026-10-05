# Run record: `<design_id>`

Copy this file into the relevant run directory and replace every bracketed value.

## Identity

- Design ID: `[DD1_062_Y32D]`
- Design family: `[DD1/DD2]`
- Parent design: `[ID or N/A]`
- Mutation(s): `[Y32D]`
- Operator: `[name]`
- Date/time (America/New_York): `[YYYY-MM-DD HH:MM]`
- Status: `[planned/running/complete/blocked]`
- Scientific run ID / attempt ID: `[exact existing ID]`
- Campaign navigation ID and manifest: `[ID/path or N/A]`
- Actual staging / raw result / analysis roots: `[paths]`
- Status evidence and observation time: `[receipt/report path, timestamp; not live state]`
- Supersedes / superseded by: `[exact ID and evidence or N/A]`
- Authoritative outputs / historical outputs: `[exact paths and scope]`

## Inputs

- Binder sequence: `[full sequence]`
- Human EGFR source: `[path, accession/version]`
- Mouse EGFR source: `[path, accession/version]`
- Structure inputs: `[paths]`
- Chain mapping: `binder=A; EGFR=B` unless explicitly remapped for display
- Residue numbering notes: `[notes]`

## Structure inference (when applicable)

- Request scope: `[prepare/predict/rerun/compare/import]`
- Engine/destination and runtime/version/checkpoint: `[values]`
- AIS Stager checkout/package version, if used: `[path/version or N/A]`
- Intake, source snapshots and input/output hashes: `[paths]`
- Ordered chains: `[role, species/isoform, source, full sequence, length, copies]`
- Constructs/mutations/tags/linkers and numbering maps: `[details/paths]`
- Seeds, sampling/recycles and per-chain template/MSA policy: `[values]`
- Ligands/PTMs/glycans/bonds/custom CCDs, retained or omitted: `[details]`
- Source authority/fallbacks and batch/grid policy: `[details or N/A]`
- Validation report and level: `[syntax/sequence/topology/chemistry/runtime parsing]`
- Runtime input/dependency paths and immutable remote snapshot: `[paths]`
- Preparation/upload/submission commands and environment: `[exact commands]`
- Job ID, status, logs and frozen configuration: `[values/paths]`
- Expected/prepared/submitted/completed/failed/missing counts: `[counts]`
- Raw structures, confidence/PAE, result index and manifests: `[paths]`
- Verified output chain mapping and sequence/species matching: `[evidence]`
- Registry refresh/validation outcome: `[status/paths or N/A]`
- Unresolved issues or information loss: `[details or none]`

## Design and filtering

- Foundry run ID: `[ID or N/A]`
- YAML/PDB inputs: `[paths]`
- Checkpoints: `[paths/versions]`
- Filter scripts and versions: `[paths/commits]`
- Filter thresholds: `[thresholds and units]`
- Retained/rejected counts: `[counts]`

## Scoring

| Method | Input | Command/API | Version/checkpoint | Output | Status/notes |
|---|---|---|---|---|---|
| Boltz2 | `[path]` | `[command]` | `[version]` | `[path]` | `[notes]` |
| DeltaForge | `[path]` | `[command/API]` | `[version]` | `[path]` | `[notes]` |
| JustHISpKa | `[path]` | `[command/server]` | `[version]` | `[path]` | `[notes]` |
| USalign | `[paths]` | `[command]` | `[binary/version]` | `[path]` | `[reference orientation]` |

## Ranking

- Selection spec/version: `[e.g. TOP20_PLAN.md v1]`
- Cutoffs: `[values]`
- Composite score definition: `[formula or N/A]`
- Rank: `[number or N/A]`

## Reporting assets

- Human PML/PNG: `[path]`
- Mouse PML/PNG: `[path]`
- Slide/deck: `[path]`
- Raw data sources: `[paths]`

## QA and exceptions

- [ ] Inputs and sequences verified
- [ ] Chain mapping verified at every tool boundary
- [ ] Human/Mouse pairing verified
- [ ] Missing values labelled
- [ ] Proxy/placeholder structures labelled
- [ ] Rendered images inspected
- [ ] Final deck opens and has expected slide count

Exceptions, manual edits, and unresolved issues:

`[free text]`
