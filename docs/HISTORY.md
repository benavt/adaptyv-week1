# Historical EGFR workspace

The eleven campaigns now live under campaigns/<historical-id>/, each with
01_staging, 02_design, 03_filtering and 04_analysis. Existing scientific IDs
and manifest bytes are preserved. The original manifest for each campaign is
04_analysis/history/manifest.json; its recorded paths are historical evidence.

Use `python3 scripts/workspace.py history [CAMPAIGN]` for the scientific overview.
Use `python3 scripts/artifacts.py list [--campaign CAMPAIGN] [--stage STAGE]`
for compact counts, or add --files for the full per-file catalog.
Use `python3 scripts/artifacts.py locate ORIGINAL_PATH --verify` for the
current canonical location and a checksum of the file present there.
The --verify flag reports current bytes; compare with artifacts.jsonl for the
migration checksum. It does not certify scientific lineage.

Shared tools, targets, registry evidence and legacy configuration remain in
reference/. Two broad folders, RFD3_filtering and PyMol_top20, remain there
because they mix reusable helpers and cross-campaign material. Do not infer
exclusive campaign ownership or rerun archived path-dependent generators.

Scientific files were renamed within one filesystem and checked against inode,
size, timestamp and mode inventories. Each moved file has a SHA-256 catalog.
Navigation symlinks were backed up and replaced by canonical locator records;
actual targets were preserved. Generated runtime build data remains in reference.
Private migration journals, inventories and exact shortcut backups live in
local/migrations/2026-10-05_primary-repository/. Reverse journaled moves only
after checking consumers; restore aliases from their backup only if needed.

Historical campaigns are not retrospectively frozen executions. New attempts
belong under campaigns/<id>/runs/ and freeze inputs/settings with parent links.
Git tracks portable source and reviewed campaign configs; artifact catalogs,
scientific data and private mappings require their own backup. The original
repository is a read-only GitHub archive, with unpublished original refs retained
locally. Never check old archive branches out over the scientific control folder.


## 2026-10-05 first EGFR candidate batch

Assembled the first batch of 55 candidate designs and a combined 70-slide
PowerPoint in the requested order:

1. Three Foundry DD1 designs: DD1_310, DD1_195 and DD1_062.
2. One DD1 Foundry design originating from ASP–HIS interface seeding.
3. One campaign-2 DD1 / Foundry DD2 fusion: D46_130_model_1.
4. Thirty ProtonPottsMPNN redesigns: ten from each of the three Foundry DD1 parents.
5. Ten LigandAI designs.
6. Ten AlphaProteinNovo designs.

Downloaded and verified Human/Mouse EGFR Boltz2 models for the LigandAI and
AlphaProteinNovo campaigns, each containing 1,000 designs and 2,000 complexes.
Selected the AlphaProteinNovo top ten by mean Human/Mouse ipSAE. DeltaForge
predicted affinity for only their 20 selected complexes; bound JustHISpKa
completed for 18 complexes. APNovo_0293 pKa remains unresolved at the user's
request. These scores are predictions, not measured binding affinities.

Corrected the ten LigandAI model slides to follow the common saved EGFR views,
with cyan design chains. Updated the campaign report, corresponding candidate
package and combined master deck, preserving the earlier report evidence.
The master is directly under `candidates/`:
`candidates/2026-10-05_EGFR_combined_candidates.pptx`, with a source-linked
`2026-10-05_EGFR_combined_candidates_manifest.json` alongside it. Candidate
subfolders retain sequences, scores, PSE sessions and supplementary information.

Submitted the top ten LigandAI designs and the top ten AlphaProteinNovo designs
for ProtonPottsMPNN redesign on Mimas. Each campaign requests 1,000 redesigns
per parent (10,000 sequences), followed by dependent Human/Mouse Boltz2 folding
(20,000 complexes per campaign). The pH 6.5 versus 7.4 preference is a design
objective and has not been established experimentally.

| Source | Submitted run | Redesign jobs | Prepare | Boltz2 array | Collector |
| --- | --- | --- | --- | --- | --- |
| LigandAI top ten | `2026-10-05_Ligand_AI_top10_pH65-vs74_1000_human_mouse_001` | 2302–2311 | 2312 | 2313 | 2314 |
| AlphaProteinNovo top ten | `2026-10-05_APNovo_top10_pH65-vs74_1000_human_mouse_001` | 2359–2368 | 2369 | 2370 | 2371 |

Submission evidence was checked in the Mimas task “Redesign EGFR binders with
MPNN” (`01a10ed8-7079-7ea2-bb7c-6b49a8a7b6d4`). The job IDs document
submission, not completion or current scheduler status. Native input validation
checked exact binder sequences, chain roles and paired species; the EGFR receptor
was fixed during redesign.

This commit records portable campaign configurations and the batch summary.
Bulk models, scientific run records, PowerPoints, PSE files and candidate
supplements remain outside Git under the existing artifact storage policy.
