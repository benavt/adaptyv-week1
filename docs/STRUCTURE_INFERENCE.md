# Structure inference workflow

Use this guide for a new prediction, input handoff, rerun, or imported prediction.
AIS Stager prepares inputs; an engine-specific runtime executes inference.
The [bundled preparation references](../01_Staging/Mimas_Folding_Skills/input_workflows.md) describes
its supported preparation APIs and known limits. Use the project's
[data contract](DATA_CONTRACT.md) and [run record](RUN_RECORD_TEMPLATE.md)
throughout. For Mimas execution, also use the pinned
[Folding Skills guide](../01_Staging/Mimas_Folding_Skills/README.md) and
[project routing guide](MIMAS_SKILLS.md).

## Resolve the request

Record the following before generating inputs, using supplied artifacts and
existing run records to resolve choices where possible:

- Stable run/design IDs and whether the request is preparation, prediction,
  rerun, comparison, or import.
- Engine and destination: native AF3, AlphaFold Server, AF2 monomer/multimer,
  ColabFold, Boltz2, or Foundry RF3; selected runtime/version/checkpoint.
- Each ordered chain's role, organism/isoform, source accession/version or file,
  actual sequence and hash, inclusive construct range, mutations, tags/linkers,
  and copy count. Distinguish a fused DD1/DD2 sequence from a DD1–DD2 assembly.
- Seeds, sampling/recycles, per-chain template/MSA policy, requested chemical
  state, ligand copies, modifications, glycan sites, bonds, and custom CCDs.
- Output destination and batch policy: ordered/unordered combinations,
  self-pairs, duplicate handling, source/fallback authority, and failed-row policy.

For the established paired EGFR route, stage binder first and receptor second;
keep `Human`/`Mouse` identity tied to source sequences. The README documents a
mouse EGFR sequence discrepancy: choose and record the intended source rather
than copying whichever representation appears first. D3 segment `n` maps to
EGFR `309 + n` only for the documented 310–501 construct. Map other constructs
explicitly, including PDB author numbering and insertion codes.

## Choose preparation and execution routes

| Request/source | Preparation reference | Execution/handoff |
| --- | --- | --- |
| Exact sequence, engineered construct, fusion or assembly | AIS input workflows §4; programmatic `IntakeSpec` with exact sequences | Generate the requested engine's format; preserve all chain copies |
| UniProt accession, natural variant | §2; `builder.build_accession_bundle` | Review annotations/state and runtime keys; multiple accessions create separate bundles |
| Complete local accession bundles with selected ranges/copies | §3; `compose_input_from_specs` | Inclusive bounds, aligned seeds, kept/dropped annotations; not an arbitrary complex importer |
| BMRB/STAR | §5; `resolve_bmrb_sequence` | Save returned resolver source/entity/chain and any fallback; check incomplete assignments |
| RefSeq, PDB chain/entity, gene, isoform, domain | §6; `references` helpers, then exact-sequence preparation | Persist chosen reference and numbering; do not infer full chemistry from polymer sequence |
| Ligands, glycans, PTMs, covalent bonds/custom CCD | §7 plus curated source job | Use a chemistry-aware native workflow and explicit residue/atom remapping; general composer support is incomplete |
| Spreadsheet or partner grid | §9 and blueprint limitations | Save per-row/job outcomes and enumeration policy; historical wrappers are dataset-specific |
| AF2/ColabFold handoff | §10 and validation guide | Expand copies; select monomer/multimer or runner query format; record omitted chemistry |
| Boltz2 | Existing project's engine runner and prior run manifest | AIS AF3 JSON is not a Boltz2 input contract; prepare/validate with the selected runner |
| Foundry RF3 | [RF3 runner reference](../02_Design/foundry/docs/rf3.md) | JSON `components[].seq` or supported coordinate inputs; FASTA is not accepted by that runner |

Locate the actual AIS Stager checkout from source provenance or discover
the relocated dependency; the Mimas guides bundle reference material only. Inspect its `AGENTS.md`, source APIs, and environment
before calling them. Reuse its package rather than copying implementation code.
Source recipes use AIS Stager's `workspace/`; this project does not adopt that
directory layout. Pass explicit output/cache paths where supported, then retain
the prepared payloads and source evidence in the Week1 run. Builder/composer
also write exports in the source repository and can overwrite names; check
those side effects and filesystem permissions before invoking them. The wizard
full-length UniProt route ignores `output_dir`; collect its actual returned paths.
Do not assume those APIs are offline unless the selected route and inputs are.

## Stage files and record provenance

Create a fresh run directory using the
[organization contract](PROJECT_ORGANIZATION.md). New runners with a configurable
input parent should use the layout below. Existing runners retain their original
`01_Staging/<run_id>/` contract; record the actual path rather than relocating
inputs or changing a runner implicitly.

```text
01_Staging/runs/<run_id>/
  RUN_RECORD.md
  intake_spec.json              resolved biological request
  manifest.json                 source choices, transformations, output hashes
  validation.json               checks performed, failures, validation level
  provenance/                   source snapshots, numbering/chain maps
  handoffs/<destination>/        alternative exports when requested
  folding/
    molecules/                  only selected runner-compatible molecular inputs
    assets/                     required MSA/template/CCD dependencies
    references/                 optional geometry references
    config.yaml                 selected runner's configuration, if applicable
    reference_map.json          if applicable
```

Create only needed files. The names of the intake/manifest/validation sidecars
are a project convention, not a claim that every AIS API emits them. A bundle's
source JSONs and manifest stay outside `molecules/`. Keep alternate native/server
JSONs out of the same discovery folder unless the runner explicitly needs both.
Give each expanded FASTA copy a unique header. Record absolute source paths and
hashes while using portable dependency paths in runtime payloads. Materialize
dependencies before transfer; the uploader rejects symlinks.

For provenance capture the resolved sequence, source choice and evidence,
numbering/chain transformation, resolver source/fallback, retained/omitted
chemistry, per-chain template/MSA settings, seeds, package/runtime versions,
exact command, output paths, and hashes. Keep non-runtime metadata in sidecars.

## Validate before execution

Apply [AIS input validation](../01_Staging/Mimas_Folding_Skills/input_validation.md):

1. Match actual sequences and inclusive ranges to the biological request; check
   mutation identities and preserve intentional tags/linkers.
2. Expand native ID lists/server counts and compare ordered sequences and copies
   across the requested exports. Inspect contents, including files named `_afs3.json`.
3. Check the target dialect/version, seeds and residue alphabet. Unknown/ambiguous
   residues need a recorded destination-specific resolution.
4. Check PTM/attachment positions, ligand identities/copies, bond endpoints/atom
   names, and bundled CCD/MSA/template dependencies. Regenerate or correctly remap
   stale templates/MSAs when sequences change; record any information loss.
5. Use the selected installed engine's parser or Server import validation when
   available. AIS emits native v4/server v1 as documented; confirm the selected
   runtime accepts the payload. Enriched `researchAnnotations` belong in a sidecar
   of a separate runtime copy, preserving the original. Changing only `version`
   is not a verified downgrade.

Report syntax, sequence/topology, chemistry/dependency review, runtime parsing,
submission, and completed prediction separately. A FASTA match cannot validate
PTMs or bonds. No new prediction was run as part of this guide's integration.

## Execute and collect

For a remote folding run, use the private configuration and
[sync guide](remote-sync.md). In the trusted private sync config, set
`LOCAL_INPUT_ROOT` to `01_Staging/runs` for the new nested layout, or to
`01_Staging` for existing root-level runs. If the shared config must
retain another input root, copy it to a private, gitignored per-run config and
select that file with `SYNC_ENV_FILE` for every command. The wrapper sources
the config after environment setup, so prefixing `LOCAL_INPUT_ROOT=...` on the
command does not reliably override a configured value.

```bash
./scripts/sync.sh doctor
./scripts/sync.sh prepare <run_id>
./scripts/sync.sh push <run_id> --dry-run
./scripts/sync.sh push <run_id>
```

Transfer does not submit inference. Inspect the chosen runner's preparation and
submission command, verify uploaded dependencies and requested settings, then
submit when execution is requested. Record the job ID, frozen inputs, environment,
checkpoint and logs. The optional `sync.sh run` executes the configured command;
inspect it before use. Do not substitute a historical launcher with hard-coded
paths for a verified runner. If runtime access is missing, preserve validated
inputs and describe the specific execution blocker.

```bash
./scripts/sync.sh pull folding <run_id> --metadata
./scripts/sync.sh pull folding <run_id> --selected --dry-run
./scripts/sync.sh pull folding <run_id> --selected
# Retain the full run for reproducibility when requested/needed:
./scripts/sync.sh pull folding <run_id>
```

Use `foundry` or the outer `campaign` pull route for those runners as documented.
Selected mode requires indexed complete predictions; it does not establish that
all jobs succeeded. Compare expected and completed job counts and retain failed
or missing outcomes. Preserve native confidence/PAE sidecars and rank identities;
do not compare raw confidence scales across engines as one shared ranking.

## Integrate results and close the run

Verify each structure's actual sequences, species, chain roles, and model identity
against the frozen request. Record the native chain mapping before conversion or
scoring. Canonical Boltz2 complexes use binder A/receptor B; other engine outputs
require inspection rather than that assumption. Put derived metrics/selections
under a run-specific `03_Filtering/` directory and retain raw downloads in
`results/`. Use the project's reporting and PyMOL specifications when requested.

After new scientific artifacts are imported, run from the project root:

```bash
python3 04_Analysis/scripts/build_model_registry.py --reuse-coordinate-cache
python3 04_Analysis/scripts/validate_model_registry.py
```

Complete the run record with prepared/submitted/completed/failed counts, exact
artifact paths, checks performed, missing metrics and unresolved issues.
Input preparation alone should be reported as prepared inputs; prediction
completion requires collected, verified engine outputs.
