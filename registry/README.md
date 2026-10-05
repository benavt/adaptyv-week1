# Model and sequence registry

Start with **[designs.csv](designs.csv)** to find a binder, or **[models.csv](models.csv)** to find a specific structure. Paths in these indexes are relative to the Week1 project root. The `home` column opens that record's folder.

## Homes

```text
registry/
  designs/<design_id>/
    metadata.json           # Identity, aliases, sequence, mutation checks, evidence
    sequence.fasta          # Recorded/extracted binder sequence
    models/                 # Relative links to all associated structure homes
  models/<kind>/<run>/<original_name>__<id>/
    <original_filename>     # Relative link to original coordinates
    metadata.json           # Provenance, chains, parent records, unknowns, warnings
    statistics.json         # Individual source-backed observations
    sequences.fasta         # Observed protein sequences from coordinates
    artifacts/              # Inputs, sidecars, manifests, preparation packages
    statistics_sources/     # Original tables and JSON responses
  runs/<run_id>/
    metadata.json
    original_run/           # Link to complete Foundry or MD run
```

Scientific files remain in their original locations so existing tools and scripts keep working. These homes add a single browsing location without duplicating coordinates or changing chains, residue numbering, or scores. They are not an independent backup: copy the whole project, including original source files, when transferring it.

## What the metadata means

- `MOD_…` identifies a **source file**, using its project-relative path. Different runs, preparations and duplicate copies remain distinct. Moving an original to a new relative path requires an explicit identity migration; ordinary refreshes retain IDs.
- `source_sha256` is the hash of the actual source-file bytes, including compression where applicable. It detects changes and byte-identical copies; it does not equate differently encoded or transformed structures.
- `SEQ_…` identifies the amino-acid string by content. Shared sequence does not imply shared structure or method.
- Existing `DD1_001`–`DD1_342` aliases come from the existing staged FASTAs and are checked against the original sequence table. They are **not reassigned by current table order**. Other RFD3 designs use run-qualified IDs to avoid collisions between Human and protonated-Human runs that reuse filenames.
- Known mutant sequences include their parent and a comparison of actual sequence changes against their mutation labels.
- `provenance` records exact per-model JSON and run-manifest facts when present. A run manifest lists intended stages; it is not proof that MPNN or RF3 ran. RFD3 outputs are labeled RFD3, not MPNN-generated sequences.
- `method_basis` distinguishes output evidence, documentation, and inference from preparation scripts or paths. Null values and `unknowns` mean missing evidence. Default commands/checkpoints and filesystem timestamps are not substituted for actual execution history.
- A recorded checkpoint path such as `rfd3_latest.ckpt` identifies the path selected by the run; it does not establish the exact weight bytes or checkpoint revision when those were not retained.
- Target species comes from exact sequence matching to staged EGFR D3. Filename disagreements are retained as warnings. Chain IDs are preserved, including display remappings and receptor-only preparations.
- Coordinate FASTAs contain **observed alpha-carbon residues**, not necessarily a complete intended sequence. Ensembles remain one file record; `coordinate_model_numbers` lists their models, and sequences describe the first model. GRO has no chain identifiers.
- `statistics.json` retains values as recorded, with a source path and CSV row number (header is row 1) or JSON pointer. Raw confidence units are retained; separate metrics tables may express different scales. Scores from different reruns are separate observations rather than a silently selected latest value.
- JustHISpKa observations extracted from preparation packages describe the **isolated EGFR chain prepared from a complex**, not the binder and not a new complex prediction. These relationships are reconstructed from the preparation script and directory identity, with target sequence verification.
- A DeltaForge response's `model_sha256` is a scorer response field, **not assumed to be the input-coordinate checksum**.
- Empty coordinate artifacts receive homes with `structure_status: empty_file`; they are not counted as valid structures. Failed/unparseable artifacts are flagged separately.

## Other indexes

| File | Purpose |
| --- | --- |
| `sequences.csv`, `sequences.jsonl` | Unique sequences and source observations |
| `sequence_inputs.csv` | Every FASTA record and its identity/status |
| `runs.csv` | Available run homes and manifest coverage |
| `source_files.csv` | Source hashes and file-state observations |
| `duplicates.json` | Byte-identical nonempty coordinate files; no copies deleted |
| `issues.json` | Unresolved sequence observations and invalid coordinate artifacts |
| `registry_manifest.json` | Generation status, scope and counts |
| `validation.json` | Latest integrity and coverage validation |
| `campaigns.csv` | Campaign navigation, available outputs, scope and manifest hashes |
| `campaign_runs.csv` | Exact recorded run/attempt IDs, actual paths, status evidence, counts and explicit supersession |
| `navigation_manifest.json` | Navigation generation scope and source-evidence hashes |
| `navigation_compatibility.json` | Organization-change audit of existing paths, model IDs and reference checks |

Model-specific warnings are in `metadata.json` and counted in `models.csv`.
`validation.json` also lists statistics species conflicts with their exact source
rows. Passing integrity checks does not resolve scientific caveats or make an
empty artifact valid.

The 20 refolded proximity-table species conflicts identified in the initial
registry were corrected on 2026-10-02 using actual receptor sequences. The
[correction audit](../03_Filtering/RFD3_filtering/out/history/2026-10-02_species_label_fix/README.md)
preserves the original table, changed rows, source hashes and ranking/slide checks.

Software fixtures, examples, installed environments and model weights are excluded. The user's ProtonPottsMPNN inference outputs inside the software submodule are explicitly included. PDB, CIF (including gzip), GRO and `.pdb.original` artifacts are indexed; Amber MOL2/topologies, arrays, logs and trajectories remain in linked preparation/run packages.

## Refresh and validate

From the project root, using Python 3.10 or newer (standard library only):

```bash
python3 04_Analysis/scripts/build_model_registry.py
python3 04_Analysis/scripts/validate_model_registry.py
```

For a faster refresh, `--reuse-coordinate-cache` reuses coordinate-derived
sequences only after verifying that each original file's SHA-256 still matches.
All provenance, statistics and relationships are rebuilt from current sources.

Refresh writes only the generated registry. It does not submit scientific jobs or recalculate scores. Source values, links, aliases and identities are validated separately. Validation checks discovery coverage, every home and link, binder/design sequence agreement, parent IDs, and source hashes. A source edited while indexing can require a refresh; it is not treated as unchanged.

Generated model/design/run homes and the verbose sequence-observation file are ignored by Git because they are reproducible from the source files and builder. The compact indexes, documentation, annotations and validation report remain available for version control. No existing Git history was rewritten.

## Campaign navigation

[Campaign pages](../03_Filtering/campaigns/README.md) group staging, original
predictions, scoring, selections, renders and reports without moving scientific
sources or changing model IDs. Their curated `manifest.json` files are the
navigation sources; generated pages and CSVs extend this registry rather than
introducing another scientific inventory. Existing `runs.csv` retains its
Foundry/MD scope. `campaign_runs.csv` adds documented preparation, submission and
scoring attempts without changing that schema.

Regenerate navigation separately from the coordinate registry:

```bash
python3 scripts/build_project_navigation.py
python3 scripts/build_project_navigation.py --check
```

The read-only check verifies reference existence, documentation links, generated
freshness and explicit successor evidence. Recorded status is local evidence,
not live remote state. Blank counts are unknown, while zero is recorded as zero.
Missing reports named by older documentation are explicitly listed as unavailable;
they are never selected as current available outputs. See the
[organization contract](../docs/PROJECT_ORGANIZATION.md) for source fields and
future migration requirements. Navigation changes alone do not require a full
scientific registry refresh.

Add operator knowledge to **[annotations.json](annotations.json)**, keyed by a model ID or design ID. It is retained across refreshes under the record's `annotations` field. Include your evidence path and distinguish confirmed facts from recollection. Do not edit generated metadata directly: a refresh replaces it.

```json
{
  "models": {
    "MOD_example": {
      "operator_note": "Explain what is known",
      "evidence_path": "path/to/original/run_record.md"
    }
  },
  "designs": {}
}
```

Older generated homes are not deleted automatically when an original file is removed. The active CSV indexes define current membership; cleanup and physical migration can be reviewed separately.
