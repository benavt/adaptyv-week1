# Campaign analysis completion and recovery

Policy recorded from the user's instruction on 5 October 2026:
**Complete means the top binding sequences for the run have been defined by
Boltz2 ipSAE, with pKa and DeltaForge predictions for those sequences.**

The working selection for the current Human/Mouse reporting pass is the top ten
distinct binder sequences by mean Human/Mouse EGFR D3 ipSAE, or every eligible
sequence if fewer than ten exist. Retain both species' individual scores. The
existing reports use the native interface ipSAE at 15 Å PAE / 15 Å distance
cutoffs; retain the directional scores, aggregation, engine and settings in the
selection evidence. Deduplicate by actual sequence, retaining a documented
representative model. Preserve run-specific constructs, chemistry and chain
roles. Never average a missing species as zero or substitute DD1–DD2 assembly
ipSAE for binder–EGFR ipSAE.

## Acceptance gates

1. **Binding ranking:** completed, validated Boltz2 complexes and native PAE for
   the declared candidate cohort; Human/Mouse sequence agreement; stable run,
   design and model IDs; exact source hashes; ranked sequence list and explicit
   exclusions. If any requested candidate lacks a prediction, label the ranking
   provisional until it is recovered or its exclusion is recorded and reviewed.
2. **pKa predictions:** results for each selected sequence/species with method,
   molecule, scope, residue mapping and preparation recorded. Use the established
   bound-complex PROPKA/JustHISpKa workflow where supported: PROPKA sidechains and
   JustHISpKa EGFR D3 HIS37/HIS100 (full EGFR H346/H409). Historical receptor-only
   values are useful controls but do not certify a bound-complex calculation.
   Reuse only a result matching the selected model and validated preparation.
   Keep conditional states separate; failed, missing or conditional-only
   unforced sites remain unresolved. A completed prediction with a documented,
   reviewed geometry warning is distinct from a failed calculation. Unsupported
   chemistry needs an appropriate preparation route; do not strip glycans or
   change disulfides to pass a protein-only runner.
3. **DeltaForge predictions:** retained raw responses and successful ΔG/Kd
   predictions for each selected sequence/species, with exact input structure,
   chain roles, source hash, scorer/runtime and units. Existing-structure scoring
   should use the selected native Boltz2 geometry after a validated conversion.
   Sequence refolding changes the scoring input and must be labelled separately.
4. **Final joined selection:** one reviewable table joining sequence, species,
   Boltz2 ipSAE, pKa values/scope and DeltaForge values to their raw evidence.
   Record exclusions, method disagreements and review outcomes. Refresh campaign
   metadata/reports and validate registry integration after new scientific
   results. No new numerical affinity or pKa cutoff is imposed by this policy;
   the existing pKa display range is not automatically an acceptance threshold.

Track **execution status** and **analysis completion** separately. A stage's
historical `complete` record means that stage finished; it does not establish all
four gates. Failed or superseded attempt IDs retain their historical outcomes
and link to a recovery run. A campaign is complete only after the final joined
selection is accepted. Predictions do not establish experimentally confirmed
binding.

## Exact flagged run IDs

The audit of all curated campaign run records and top-level staging status JSONs
finds four indexed attempts with the requested status labels, all in
`dd1_dd2_stitch_afsample3_20261004`.

| Run ID | Recorded status / cause | Recovery action |
| --- | --- | --- |
| `2026-10-04-DD1_DD2_Stitch_EGFR_AFS3_submit_001` | `prepared_not_submitted`; superseded when the requested sample count changed | Preserve as superseded; follow 002 → 003 → corrected 005. Do not submit the obsolete preparation. |
| `2026-10-04-DD1_DD2_Stitch_EGFR_AFS3_submit_002` | `prepared_not_submitted`; saved reason: remote execution approval declined; later explicitly superseded by 003 | Preserve its receipt/history. Continue through the already-executed corrected retry. |
| `2026-10-04-DD1_DD2_Stitch_EGFR_AFS3_submit_003` | `failed_before_prediction`; all five requests failed because EGFR MSA/template specification was partial | Retry 004 supplied the sequence-matched cached EGFR MSA. Corrected 005 incorporates that fix. |
| `2026-10-04-DD1_DD2_Stitch_EGFR_AFS3_submit_004` | `failed_before_prediction`; all five requests failed because the custom template lacked `revision_date` | Retry 005 supplied conversion-date metadata and passed native data-pipeline/featurisation preflight. Continue with its outputs. |

**Mimas check, 5 October 2026:** corrected
`2026-10-04-DD1_DD2_Stitch_EGFR_AFS3_submit_005` has **5,000 native sample
`model.cif` files**, plus five best-model copies. All five task logs finish with
`Done running 1 fold jobs.` and report seed 200. This establishes generated
remote predictions; collection, sample-by-sample molecular/sidecar validation
and registry import are still pending in this audit. The earlier saved
`running`/zero-model snapshot is historical.

For analysis completion, the same five fusion sequences **already have paired
protein-only Boltz2 ranking** in the
[stitched selection](../03_Filtering/campaigns/dd1_dd2_stitch_afsample3_20261004/top10_egfr_20261005/selection.json).
The leading sequence is
`design_DD1_pair11_3_model2_Nterm_KHK_E2_D46_127_model_1`, mean ipSAE **0.809086**.
Score these ten native complexes with pKa and DeltaForge, then publish the joined
five-sequence selection. All ten currently lack those values in the report.
The 5,000 AFSample3 samples are a separate glycosylated Human-only protocol;
their collection does not replace the paired Boltz2 evidence or silently change
the chemistry of either analysis.

## Other unfinished work relevant to completion

- `2026-10-04_human_rfd3_rank1_pka_plan_001` is recorded as
  `prepared_not_executed`, but execution attempts `...pka_batch_001` and
  `...pka_batch_002` subsequently exist. The latter is marked
  `propka_outputs_removed`; its removal receipt records user-requested deletion,
  retained JustHISpKa files and a suspended batch. The former's saved `running`
  flag is not proof of a live process. Do not resume either from its status label
  alone. Audit reusable JustHISpKa checkpoints by exact source hash, and prepare a
  fresh focused top-ten sequence batch from the updated **10,000-pair** ranking.
  The original plan covers only 8,693 pairs. Resuming its combined runner would
  regenerate PROPKA outputs and requires reviewing that intended scope first.
- `DD2_DD1_top10_HK_linker_20261005_001` failed before RFD3 because Bio.PDB was
  unavailable. Attempt 002 generated 100 compressed backbones for each of the
  first four parents; its MPNN handoff failed on `use_auth_fields` versus native
  `use_author_fields`. At 14:24 UTC, job **2057** still has parents 5–10 queued;
  the MPNN-only retry `..._003`, job **2077**, waits on `afterany:2057`. Retain and
  finish those backbones without rediffusing completed parents. Audit 100 unique
  backbones per parent, then exact domains, one-chain DD2–linker–DD1 topology,
  5–20 linker residues and the recorded H/K-rich acceptance. Fold accepted
  fusions as **one binder chain plus one EGFR chain**, separately for Human and
  Mouse, rank by Boltz2 ipSAE, and add pKa/DeltaForge. No three-protein-chain
  Boltz2 systems are requested.
- `Ligand_AI_Human_Mouse_Boltz2_20261005_001` is an additional recorded campaign
  outside the ten PPTX cohorts. Its saved submission is 2,000 complexes, Slurm
  2071 plus collector 2072. This audit does not establish its latest prediction
  completeness. Collect/validate existing outputs before considering retries,
  then establish the top-ten paired ranking and scoring cohort.

## Campaign scoring coverage

Counts below mean **numeric fields present in the frozen report**, not verified
scientific completion or a claim that no other raw scoring output exists.

| Campaign | Selected sequences / paired candidate pool | Selected complexes | DeltaForge ΔG/Kd fields | HIS37/HIS100 pKa fields | Main remaining step |
| --- | ---: | ---: | ---: | ---: | --- |
| Putative DD1 | 10 / 342 | 20 | 20 | 20 | Validate DeltaForge input provenance; historical pKa is receptor-only. Obtain validated bound-complex pKa. |
| DD1 ASP | 10 / 400 | 20 | 0 | 0 | Score both species for the selected ten; existing parent-focused pKa is not the whole top-ten cohort. |
| DD2 pair11 fusions | 10 / 1,000 | 20 | 0 | 0 | pKa and DeltaForge on existing two-protein-chain complexes. |
| DD1_062 mutants | 2 / 2 | 4 | 4 | 4 | Validate DeltaForge and replace/relabel receptor-only pKa for bound-context acceptance. |
| Stitched AFSample3 cohort | 5 / 5 | 10 | 0 | 0 | pKa and DeltaForge; collect completed AFSample3 separately. |
| DD1_062 PROPKA comparison | 1 / 1 | 2 | 2 | 2 | Reuse validated Human bound pKa with warnings; add Mouse bound pKa and join existing DeltaForge. |
| ProtonPotts existing EGFR smoke binders | 3 / 3 | 6 | 6 | 6 | Scope/prepare bound pKa; keep PD-L1 examples and unfinished EV6 labels outside this EGFR ranking. |
| Separate DD1–DD2 RF3 | Pending / 0 EGFR pairs | 0 | 0 | 0 | Finish linker design → paired EGFR Boltz2 → ranking → pKa/DeltaForge. |
| Reference binder | 1 / 1 | 2 | 0 | 0 | Resolve Human glycan geometry flag; use chemistry-aware pKa/DeltaForge preparation. |
| Human RFD3 paired evaluation | 10 / 10,000 | 20 | 0 | 0 | Check retained pKa checkpoints; new top-ten focused scoring batch. |

Every selected sequence is retained with exact Human/Mouse native structures,
full binder/target sequences and input hashes in the
[machine-readable audit](../03_Filtering/workflows/campaign_completion_20261005/audit.json).
Each campaign has an `ANALYSIS_COMPLETION.md` companion with its next action.
Overlap between campaigns can reuse validated results on the **same molecule,
coordinates, method and preparation**; a shared family name alone is insufficient.

## Extraction and refresh

Run `python3 scripts/snapshot_campaign_completion.py` from the project root to
refresh the local audit, recovery CSV and ten campaign companions. It reads each
campaign's curated `manifest.json`, follows the exact status evidence and JSON
pointer, and hashes the original evidence. It also audits finite numeric
DeltaForge ΔG/Kd and both histidine pKa fields in the retained report snapshot.
It launches no scorers, predictions or remote commands. Remote observations are
explicitly dated evidence supplied separately; rerunning the extraction does
not refresh cluster state.

The output lists every exact flagged ID, native failure cause, explicit
supersession where recorded, recovery route, selected input sequences/hashes and
completion gaps. Textual flags and JSON job outcomes are read separately; blank
counts remain unrecorded. Retained pKa pipeline source shows receptor extraction
and distinguishes controls from the newer NRC-corrected bound calculation.
Numeric display coverage alone never certifies complete. Native historical run
records are retained; this audit adds the user's final analysis criterion.

Recovery data: [CSV](../03_Filtering/workflows/campaign_completion_20261005/recovery_runs.csv),
[audit and evidence hashes](../03_Filtering/workflows/campaign_completion_20261005/audit.json),
[Mimas observation](../03_Filtering/workflows/campaign_completion_20261005/mimas_observation.json).
Execution recipes: [complex pKa](BATCH_COMPLEX_PKA.md),
[DeltaForge](../03_Filtering/DeltaForge/README.md),
[structure inference](STRUCTURE_INFERENCE.md), [Mimas routing](MIMAS_SKILLS.md).
