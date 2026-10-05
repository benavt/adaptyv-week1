# EGFR design data contract

## Identity and naming

- Use a stable design ID such as `DD1_062`, `DD1_062_Y32D`, or `DD1_062_Y32D_T3H`.
- Encode mutations as `<original><position><new>` and join multiple mutations with `_`.
- Use `Human` and `Mouse` exactly for species labels.
- Use `DD1` and `DD2` exactly for design families.
- Do not use directory order or display rank as the primary identity.

Recommended artifact names:

```text
<design_id>_<species>_<artifact>.<ext>
```

## Sequence contract

- Store the complete binder sequence, not only the mutation string.
- Record the target sequence source and version/date.
- For paired EGFR complexes, binder sequence precedes EGFR D3 in FASTA inputs.
- Preserve the sequence used for scoring alongside the result table.
- A mutation label must be checked against the actual sequence before analysis.

## Chain contract

For Boltz2 DD1/EGFR complexes, the canonical mapping is:

| Role | Chain |
|---|---|
| Binder/DD1 or DD2 | A |
| EGFR D3/receptor | B |

DeltaForge must therefore receive receptor `B` and peptide `A` for these complexes.

Some display-preparation artifacts intentionally reverse the chains for a PyMOL-only display. Such files must be labelled as display remappings and must not be used for metric extraction without an explicit remapping step.

## Residue contract

- Preserve source chain IDs and residue numbers during preparation unless a tool requires a documented remapping.
- Record both segment-local and full-protein numbering when relevant. For the current EGFR D3 staging, segment residue `n` corresponds to EGFR residue `309 + n`.
- The project histidine landmarks are D3 `HIS37` and `HIS100` (full-protein H346 and H409 in the documented staging structures).
- Distances must state whether they use heavy atoms, all atoms, or representative atoms and must include units.

## File and provenance contract

- Model/sequence homes and indexes are maintained in `registry/`; see
  `registry/README.md` for scope, schema meanings and refresh commands.
- Record structure identity separately from binder identity and ranking. One
  binder can have several prediction, preparation and species-specific structures.
- Keep source-file checksums, exact input evidence and statistics row/JSON locators.
- Add operator-supplied provenance in `registry/annotations.json`, with evidence;
  generated metadata is refreshed from original artifacts.

- Never overwrite an original PDB, CIF, FASTA, or external response.
- Save converted intermediates beside the source or in the run-specific directory.
- Keep raw JSON responses and normalized CSVs together.
- Record tool, version/checkpoint, command, environment, date, and source paths.
- Label placeholders, proxies, failed runs, and manually corrected structures explicitly.

## Missing values

- Use `N/A` in human-facing reports when a metric is unavailable.
- Use empty CSV fields only when the schema requires a numeric field to remain blank.
- Never infer one metric from another or substitute a different model/score without recording the substitution.
