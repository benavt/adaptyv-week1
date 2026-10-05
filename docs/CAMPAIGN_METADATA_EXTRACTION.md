# Campaign metadata extraction and updates

Initial local extraction: 2026-10-04T19:51:44-04:00. This document describes the nine `CAMPAIGN_METADATA.md` companions and how to refresh their evidence. Metadata is curated Markdown; scientific artifacts, scores, run IDs and data paths are unchanged.

## Snapshot inventory

| Campaign | Local recorded status | Metadata |
| --- | --- | --- |
| Paired DD1 evaluation against Human and Mouse EGFR | Paired prediction/metrics and reports available; historical job provenance incomplete | [putative_dd1_egfr_d3](../03_Filtering/campaigns/putative_dd1_egfr_d3/CAMPAIGN_METADATA.md) |
| DD1 ASP placement and paired Boltz2 evaluation | 400 designs / 800 paired Boltz2 complete in saved records; full pKa collection and final reports unfinished | [dd1_asp_20261003](../03_Filtering/campaigns/dd1_asp_20261003/CAMPAIGN_METADATA.md) |
| DD2 N-terminal DD1 extensions: Boltz2 monomers | 1,000 stitched designs and 3,000 Boltz2 predictions complete in saved records; monomer report available | [dd2_pair11_20261003](../03_Filtering/campaigns/dd2_pair11_20261003/CAMPAIGN_METADATA.md) |
| Separate DD1–DD2 complexes: RF3 top-5% reports | Top-5% reports complete for two targets; broader RF3 structure coverage unfinished | [dd2_rf3_top5pct_20261003](../03_Filtering/campaigns/dd2_rf3_top5pct_20261003/CAMPAIGN_METADATA.md) |
| DD1_062 mutant species comparisons | Four mutant/species predictions, DeltaForge and pKa summaries available; ipSAE and gate metadata incomplete | [dd1_062_mutants](../03_Filtering/campaigns/dd1_062_mutants/CAMPAIGN_METADATA.md) |
| EGFR D3 references and historical MD | Reference artifacts and two finished 10 ns MD logs available; provenance and trajectory QA remain to be consolidated | [egfr_d3_references_md](../03_Filtering/campaigns/egfr_d3_references_md/CAMPAIGN_METADATA.md) |
| Stitched DD1/DD2–Human EGFR AFSample3 attempts | Attempt 005 submitted/running in saved observation; 0 completed models recorded; live status needs checking | [dd1_dd2_stitch_afsample3_20261004](../03_Filtering/campaigns/dd1_dd2_stitch_afsample3_20261004/CAMPAIGN_METADATA.md) |
| DD1_062 Human PROPKA comparison | Audited smoke 002 complete for two local PROPKA inputs; no unfinished jobs in recorded scope | [dd1_062_propka_20261004](../03_Filtering/campaigns/dd1_062_propka_20261004/CAMPAIGN_METADATA.md) |
| ProtonPottsMPNN existing outputs and smoke-test evidence | Packaged PD-L1 example outputs available; EGFR/DD2 EV6 smoke blocked before final labels | [protonpottsmpnn_existing_outputs](../03_Filtering/campaigns/protonpottsmpnn_existing_outputs/CAMPAIGN_METADATA.md) |

## Uniform fields and units

Every companion uses the same sections: scope/targets, inputs/parameters, run/job ledger, filtering/selection, available outputs/validation, unfinished work/unknowns, and evidence catalog. The header records stable campaign ID, extraction time, status and observation scope. Each evidence row supplies a project-relative file link, locator (JSON pointer, CSV field, log line or document section) and SHA-256 of the local source bytes.

A **design** is a source sequence/geometry candidate; a **request** is one complete molecular input; a scheduler **task** may process many requests; a **model** is one seed/sample output; a **site row** is a residue observation. Species contexts, seeds, samples, table rows and duplicate source copies are never interchanged. Blank/unknown means unrecorded; N/A means the field does not apply; zero means explicitly observed or recorded zero. Expected minus locally observed is a collection/evidence gap, not an automatic remote failure count.

These nine navigation campaigns overlap in source dependencies and are not nine disjoint experiment batches. The ProtonPottsMPNN campaign includes a PD-L1 packaged example as well as a blocked EGFR/DD2 smoke test; neither is silently merged with completed EGFR design runs.

## Evidence and precedence

1. Start from `03_Filtering/campaigns/*/manifest.json` for the nine identities, paths and authoritative report revisions. Those manifests are navigation sources, not complete scientific run records.
2. Read native intake/FASTA/JSON/YAML, run manifests, frozen commands, submission receipts, native completion summaries, log exits, scoring definitions and selection manifests. Read sequence contents and chemistry; never infer topology or species from a family name or filename. Check exact sequence/chain correspondences before transferring a score or status.
3. Prefer source-backed executed-run records over superseded local preparation, while retaining both. Pair11's unsubmitted separate-binder preparation differs from the executed single-chain fusion. AFSample3 prepared 001/002 and failed 003/004 do not replace latest submitted 005. PROPKA audited 002 supersedes 001. Native final completion summaries outrank stale intermediate counters, without deleting those counters.
4. Separate preparation, transfer, scheduler acceptance, runtime parsing/featurisation, completed predictions, collection, scientific validation and reporting. A script documents configured parameters; it is not an execution receipt. A fold manifest is not a completed fold. A slide deck's presence does not close an upstream scoring campaign.
5. Label counts **source-reported** when copied from manifests/summaries and **locally measured** when counting rows/files. This extraction measured table coverage, selected inputs and completion markers; it did not recalculate scientific scores, validate every trajectory/model, query remote queues or launch jobs. Embedded source timestamps take precedence over extraction date for execution status. MD log timestamps without timezone remain unspecified.
6. Record conflicts and unknowns explicitly. Generic `latest` checkpoint paths do not identify weight bytes. DeltaForge response model hashes identify scorer weights, not coordinate hashes. An unassigned gate is not failed binding. Native AF3 and AlphaFold Server inputs must be identified by contents, not folder label.

## Extraction recipes

Run read-only inspection from the project root. Use `rg --files` and narrowly scoped `rg` searches to locate files; avoid recursively mixing current and superseded output directories. Python standard-library JSON/CSV parsing was used for field reads and coverage; SHA-256 was computed on exact file bytes. Structure/sequence claims come from native input and retained run validation, not freshly executed molecular calculations.

CSV/TSV data rows exclude the header; JSON list entries count records, not distinct molecular sequences unless deduplicated explicitly:

```python
import csv, json, hashlib
from pathlib import Path
csv_rows = list(csv.DictReader(Path("03_Filtering/Boltz2_metrics/out/top20_grid_search_candidates.csv").open()))
assert len(csv_rows) == 20  # 1,066-row grid table is a different artifact
sweep_rows = list(csv.DictReader(Path("03_Filtering/ProtonPottsMPNN/inference/outputs/sweep_designs.tsv").open(), delimiter="\t"))
print(len(csv_rows), len(sweep_rows))
print(hashlib.sha256(Path("results/dd2_rf3_top20_20261003/coverage_summary.json").read_bytes()).hexdigest())
```

For DD1 ASP pKa, use **only** the current output root; join marker identities to the 800 expected coordinate-source hashes. Markers prove locally recorded completion for that preparation, not complete conditional workflows or fresh scientific validation. Do not add priority-site rows or superseded-preparation markers to this job count:

```python
import json, math
from pathlib import Path
root = Path("03_Filtering/JustHIpKA/DD1_ASP_Boltz2_20261003")
jobs = json.loads((root / "manifest.json").read_text())["jobs"]
expected = {job["sha256"] for job in jobs}
markers = [json.loads(p.read_text()) for p in (root / "outputs").glob("*/complete.json")]
matching = {m["sha256"] for m in markers if m.get("sha256") in expected}
valid_values = [m for m in markers if m.get("sha256") in expected
                and len(m.get("pKa", [])) == 2
                and all(isinstance(v, (int, float)) and math.isfinite(v) for v in m["pKa"])]
print(len(jobs), len(expected), len(markers), len(matching), len(valid_values), len(expected - matching))
print(sum(bool(m.get("conditional")) for m in markers))
# Initial snapshot: 800, 800, 59, 59, 59, 741; nonempty conditional=0.
```

RF3 coverage is read from the four branch dictionaries in coverage_summary.json, whose original analyzer joins exact DD1/DD2 sequence pairs. Sum unique_designed_pairs, has_RF3_structure and missing_RF3_structure separately; do not combine source-structure and RF3 gaps. Fraction selection uses ceil(0.05 × 1,040) = 52 per target. Pair11 native summary provides 1,000 predictions/context; AFSample3 expected count is five requests × 200 seeds × five samples = 5,000; the earlier interpretation was 200 total samples/request.

For GROMACS, inspect native mdp dt/nsteps and log engine/version/finished lines; check nonempty native files separately. A finished log does not establish trajectory stability. For PROPKA, retain sentinel/terminus exclusions and bound-versus-isolated context. For mutants, parse species/site rows and raw DeltaForge gate fields; do not fill missing ipSAE from ipTM.

## Shared receptor identity

These values identify the prepared 192-residue EGFR 310–501 construct, not every historical file named EGFR. For this exact construct, local residue n maps to full receptor 309+n. Campaign chain IDs vary and are stated in each companion. Native JSON includes chemistry; sequence hashes alone cannot verify glycans/bonds.

| Construct | Sequence length | SHA-256 of amino-acid sequence only | Native input |
| --- | --- | --- | --- |
| Human / Homo sapiens | 192 | `619f4e93960e69b0d6e37b312fafdddb76975a101bec4934d49f1871f6663002` | [Human native JSON](../01_Staging/EGFR/Human_EGFR_310_501_afs3.json) |
| Mouse / Mus musculus | 192 | `251b388a69a5747da9af623ca62a43d0ae3ed304e472b6c4ee10aa38b0f3b816` | [Mouse native JSON](../01_Staging/EGFR/Mouse_EGFR_310_501_afs3.json) |

No organism isoform accession is asserted from this sequence hash alone; record the exact source accession/isoform when preparing future intake. Exact supplied constructs remain authoritative. Root README's Mouse full-receptor 598 discrepancy does not authorize altering these existing D3 inputs.

## Refresh and verification

Re-extract when new receipts/results arrive. Check existing remote job state when reconciling execution or planning reruns; for this snapshot, all status is local saved evidence. Update the relevant campaign's ledger, counts, source observation time, extraction time, output authority and unfinished actions together. Preserve supersession/failure history and distinguish absent local evidence from remote outcomes. Recompute every changed evidence hash after reading the new source; never update hashes alone to disguise unread changes.

Run `python3 scripts/check_campaign_metadata.py` to check nine companions, uniform headings, evidence-file hashes, E-number references and local links. This detects changed source files and broken references; it does **not** automatically refresh prose or validate molecular predictions. Re-run the count recipes for directory/table-derived facts because directory contents are not covered by a manifest's file hash.

Run `python3 scripts/build_project_navigation.py` then `python3 scripts/build_project_navigation.py --check` after adding/updating companions. The navigation builder exposes these files and tracks their hashes; it does not regenerate curated metadata. Original scientific files and model IDs must retain their paths. After importing new scientific models, separately refresh/validate the model registry using the project's existing workflow.

Unfinished-work checkboxes are a backlog for the next campaign work, not authorization or proof of reruns. First reconcile existing jobs/outputs, then use project structure-inference/run-record instructions for any execution. No additional prediction, scoring or remote submission was performed to create these companions.

