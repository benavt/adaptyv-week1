# Canonical EGFR–DD1/DD2 workflow

## Purpose

This is the canonical path for taking a binder design or mutant from a staged sequence to a reproducible Human/Mouse EGFR evaluation and report. Every stage should leave inputs, commands, outputs, and a run record behind.

## Stage map

1. **Stage and validate inputs** in `01_Staging/`.
2. **Generate designs** with Foundry RFD3/MPNN/RF3 in `02_Design/`.
3. **Filter structural candidates** with RFD3 sanity filters and interface/proximity filters.
4. **Refold and index paired complexes** with Boltz2; preserve Human/Mouse pairing by sequence, not directory name.
5. **Score candidates** with DeltaForge, JustHISpKa, and structural comparison tools as applicable.
6. **Rank candidates** using the versioned selection heuristic recorded in the run record.
7. **Render structures** with the project PyMOL conventions.
8. **Build reports/slides** from normalized CSV/JSON outputs and rendered assets.
9. **Perform QA and archive** the final deliverables and provenance.

## Mimas execution guides

For design generation and folding/refolding on Mimas, follow
[MIMAS_SKILLS.md](MIMAS_SKILLS.md). It routes agents to the pinned Design Skills
and Folding Skills submodules, their remote access/native input guides and
project artifact locations. Apply the readiness gates below and record both
reference commits and actual runtime identities; reading a guide or preparing
inputs does not establish completed execution.

## Stage gates

Structure prediction and input handoff requests use
[STRUCTURE_INFERENCE.md](STRUCTURE_INFERENCE.md), including AIS Stager source
selection, constructs/copies, output dialects, chemistry/dependencies, and
engine-specific runtime validation. Its intake and execution record extends
Gate 1 and the refolding stage. Select the requested engine rather than treating
Boltz2 as mandatory for every inference job; the established paired EGFR
evaluation still uses the Boltz2 route below.

### Gate 1 — Input readiness

- Binder sequence and mutation notation are recorded.
- Human and Mouse EGFR sequences are explicitly identified.
- Source structures are copied without overwriting originals.
- Chain IDs, residue numbers, and intended binder/target roles are recorded.
- PDB/CIF atom naming is compatible with the next tool.

### Gate 2 — Design readiness

- Foundry YAML and PDB agree on the target and hotspot residues.
- The intended designed chain is explicit.
- Run ID is unique and the job manifest is retained.
- Checkpoint, environment, Slurm job ID, and source inputs are recorded.

### Gate 3 — Filter readiness

- Filter thresholds and distance definitions are recorded.
- A passing filter is not described as evidence of binding.
- All retained and rejected counts are saved.
- Human/Mouse results are paired by stable design ID and sequence.

### Gate 4 — Scoring readiness

- Each external or local score has a source file and method/version.
- Missing values are represented as `N/A` or blank according to the data contract; they are never silently substituted.
- Chain arguments are checked before DeltaForge or metric extraction.
- Proxy, placeholder, or failed protonation preparations are labelled.

### Gate 5 — Reporting readiness

- Both species are rendered with the appropriate fixed views.
- Required score rows and units are present.
- The sequence, color key, and method note are included.
- Slide dimensions preserve image aspect ratios.
- The final deck opens and contains the expected number of slides.

## Canonical project paths

See [PROJECT_ORGANIZATION.md](PROJECT_ORGANIZATION.md) for campaign navigation,
new-work layouts, source authority and path compatibility. Existing artifacts,
runner contracts and pinned submodule paths are unchanged.

| Function | Location |
|---|---|
| Staged inputs | `01_Staging/runs/<run_id>/` for configurable new runners; existing `01_Staging/<run_id>/` contracts retained |
| Reusable references | `01_Staging/references/`; existing EGFR sources remain in `01_Staging/EGFR/` |
| Foundry design runs | `02_Design/foundry/runs/` |
| Filtering and scoring | `03_Filtering/` |
| Campaign navigation and new analysis packages | `03_Filtering/campaigns/<campaign_id>/` |
| Native remote downloads | Existing sync routes under `results/`; historical prediction packages stay in place |
| Software navigation and new installations | `tools/`; existing pinned sources and environments retained |
| Analysis utilities | `04_Analysis/` |
| Per-run metadata | The actual run directory, using `docs/RUN_RECORD_TEMPLATE.md`; existing records retained |
| Model/design homes and sourced statistics | `registry/`, described in `registry/README.md` |
| Operator-supplied provenance | `registry/annotations.json`, with evidence paths |

## Registry refresh

After importing structures, producing statistics, or staging designed sequences,
refresh the additive registry and validate it from the project root:

```bash
python3 04_Analysis/scripts/build_model_registry.py --reuse-coordinate-cache
python3 04_Analysis/scripts/validate_model_registry.py
```

The registry links to existing scientific files and does not replace run records.
It preserves source IDs, checksums, chains and individual statistical observations.
Unknown execution details stay unknown; complete a run record for future work
rather than relying on reconstruction from filenames.

When campaign evidence or available reports change, refresh the additive
navigation indexes and check their references:

```bash
python3 scripts/build_project_navigation.py
python3 scripts/build_project_navigation.py --check
```

Navigation-only changes do not require rebuilding the scientific registry or
running scientific jobs. Recorded status is local source evidence, not live
remote state.

## Reproducibility rule

A result is deliverable only when another operator can identify the exact source structure/sequence, command or script, environment, threshold/version, and output artifact used to produce it.
