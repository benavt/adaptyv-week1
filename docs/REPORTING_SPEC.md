# EGFR–DD1/DD2 reporting specification

## Scope

Use this specification for a one-design summary or a ranked multi-design deck. The current DD1 top-20 implementation is `03_Filtering/PyMol_top20/build_dd1_top20_deck.py`, with the reference deck at `03_Filtering/PyMol_top20/DD1_top20_EGFR_interface_renderings_v1.pptx`.

## Required score table

Use consistent rows and units:

- ipSAE
- ipTM, when available
- complex pLDDT, when available
- DeltaForge ΔG (kcal/mol)
- DeltaForge Kd (nM)
- HIS37 pKa
- HIS100 pKa
- Mouse-versus-Human whole-complex TM-score and RMSD (Å)

Use `N/A` when unavailable. State the USalign reference orientation. Do not silently substitute a different model or score.

## Slide layout

- 16:9, light/white background.
- Title at top with design ID and Human/Mouse EGFR D3.
- Human panel left; Mouse panel right.
- Equal displayed panel dimensions with preserved aspect ratio.
- Complete binder sequence above the color key/table.
- Color key and method note outside the table.
- No explanatory text over molecular images.

## Deck structure

For a top-20 deck, use:

1. Selection-method/grid-search summary.
2. Aggregate pKa/DeltaForge summary.
3. One slide per ranked design with paired Human/Mouse renderings and score table.

The reference builder creates 21 slides for a top-20 deck.

## QA checks

- [ ] Correct design ID and rank.
- [ ] Human and Mouse panels are present and not stretched.
- [ ] Sequence matches the source CSV.
- [ ] Chain mapping is consistent with the data contract.
- [ ] Required rows, units, and missing-value labels are present.
- [ ] pKa values use the documented precision and 6.5–7.4 heatmap range.
- [ ] TM-score/RMSD reference orientation is stated.
- [ ] Color key is present and accurate.
- [ ] No table, caption, or legend overlap.
- [ ] Deck opens and contains the expected slide count.
