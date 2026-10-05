# EGFR campaign workflow

Use one project (`EGFR`), a campaign for each scientific question, and immutable
runs for its attempts. The four stage names match locally and remotely:

| Stage | EGFR responsibility | Execution |
| --- | --- | --- |
| `01_staging` | Validate constructs, topology and hotspot inputs | Mac |
| `02_design` | Generate designs with the selected campaign engine | Mimas |
| `03_filtering` | Selected predictions/refolding, scoring and selection | Mimas, then Mac |
| `04_analysis` | Compare, rank, plot and report | Mac |

New campaign scientific choices are explicit. Commands below show the interface; placeholders do
not choose scientific settings or create evidence.

## Prepare and transfer

```sh
python3 scripts/workspace.py show
python3 scripts/workspace.py show --stage 02_design
python3 scripts/workspace.py campaign CAMPAIGN --title "TITLE" --objective "QUESTION"
python3 scripts/workspace.py new CAMPAIGN LABEL --input-dir drafts/NATIVE_INPUTS --settings drafts/settings.json
python3 scripts/workspace.py verify CAMPAIGN RUN
python3 scripts/workspace.py upload CAMPAIGN RUN
python3 scripts/workspace.py upload CAMPAIGN RUN --apply
python3 scripts/workspace.py collect CAMPAIGN RUN --files-from local/selected-files.txt
python3 scripts/workspace.py collect CAMPAIGN RUN --files-from local/selected-files.txt --apply
python3 scripts/workspace.py overview
python3 scripts/workspace.py standard presentations
```

Use `--parent-run CAMPAIGN/PRIOR_RUN` for successors; repeat for multiple origins.
Materialize exact selected native inputs explicitly. Parent links do not copy
bulky results. Freeze changed inputs/settings in a new date-first run; previous
runs stay intact. Frozen inputs include project, campaign, settings and workflow
pins. Verification checks hashes and identities, not scientific correctness.

`run.json` and `01_staging/` inputs are immutable. `state.json` records progress;
receipts hold separate transport/execution observations. The front door does not
submit jobs or infer native settings. Validate actual runtime paths, checkpoints,
constructs and native inputs before first selected-engine execution in this layout.

Upload checks remote target identity, verifies hashes and retains the frozen
remote preparation record. Collection requires explicit relative files under
`02_design/`, `03_filtering/` or `04_analysis/`, verifies hashes, and refuses
existing local output paths. Bulk predictions stay remote; collect the small
FASTA/score files needed for local trace verification and analysis. Record private
canonical/archive locations in receipts, not portable source.

## Trace candidates across campaigns

After collecting validated native outputs, prepare records in `drafts/` using the
pinned shared `examples/candidate.example.json` and `observation.example.json`.
Templates contain placeholders; replace them with actual source evidence.

```sh
python3 scripts/workspace.py candidate register --record drafts/candidate.json
python3 scripts/workspace.py candidate observe --record drafts/observation.json
python3 scripts/workspace.py candidate trace --id CANDIDATE
```

A candidate records its exact sequence, native FASTA identity/file hash, origin
campaign/run, generation method and parent candidates. A changed sequence creates
a new candidate with its native source and parent links. Missing generation fields
remain explicit; supplied metadata is not independent execution certification.

An observation links one candidate to an exact score artifact and JSON Pointer or
CSV row/column, method/version, native value and scientific scope. Include relevant
species, construct, pH and assembly. Keep predictions distinct from lab results.
Registry records are immutable; trace validates source hashes and recursive parent
records before final ranking.

## Freeze the chosen set

Write candidate IDs one per line in intended rank order. Ranking JSON must contain
`method`, `version` and `rationale`. Include thresholds, manual decisions and the
actual observations used; optional `evidence` entries link `candidate_id` and
`observation_id`. Missing ranking evidence must remain explicit.

```sh
python3 scripts/workspace.py selection --ids-file drafts/ranked-ids.txt --ranking-record drafts/ranking.json
```

The new `submissions/YYYY-MM-DD_selection_NNN/` contains ranked CSV and frozen
candidate/parent/observation traces with ranking rationale and project competition
metadata. Later observations do not change it. Duplicate candidate IDs or sequence
hashes are rejected. Another ranking or selected set needs a new package.

The package is a draft. Researcher review must check eligibility, track limits,
scientific interpretation and final order against docs/CHALLENGE.md. Nothing is
submitted automatically; retain the actual submission receipt separately.

Runs, registry records and selection packages are ignored by current Git policy.
Back them up with referenced artifacts in the recorded temporary archive; code on
GitHub alone does not preserve scientific lineage. Read the shared
`docs/provenance.md` for record details. Old guides under `reference/` remain
historical source material with prior path contracts, not new execution defaults.
