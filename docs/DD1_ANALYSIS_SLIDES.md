# Creating DD1 analysis slides

This guide documents the repeatable workflow behind `03_Filtering/PyMol_top20/DD1_top20_EGFR_interface_renderings_v1.pptx`.

## Inputs

Required normalized inputs include:

- `03_Filtering/Boltz2_metrics/out/top20_grid_search_candidates.csv`
- `03_Filtering/Boltz2_metrics/out/boltz2_dd1_metrics_long.csv`
- Human and Mouse DeltaForge results
- JustHISpKa results
- Mouse-versus-Human USalign results
- Human/Mouse rendered PNGs and fixed-view JSON files

The builder also converts SVG summary assets to PNG for PowerPoint insertion.

## Build sequence

1. Run the Boltz2 metric analysis and top-20 selection.
2. Run Human/Mouse structure rendering with `render_top20_egfr_interfaces.py`.
3. Generate the pKa/DeltaForge summary assets.
4. Confirm all per-design score keys exist for both species.
5. Run `build_dd1_top20_deck.py`.
6. Open the resulting PPTX and inspect representative slides, including the first summary slide, one high-ranked design, and the final design.

Typical commands from the repository root:

```bash
python3 03_Filtering/PyMol_top20/render_top20_egfr_interfaces.py
python3 03_Filtering/PyMol_top20/build_dd1_top20_deck.py
```

## Slide template

- 16:9 canvas.
- Light background.
- Title and rank at top.
- Human EGFR D3 image on the left; Mouse EGFR D3 image on the right.
- Sequence line below images.
- Color key below sequence.
- Compact three-column score table at bottom.

The reference implementation uses 1400×1000 rendered PNGs, places each at equal dimensions, and creates 22 slides: two summary slides plus one slide per top-20 design.

## Per-design data rows

The reference table includes ipSAE, TM score, DeltaForge ΔG, DeltaForge Kd, cross-species TM score, HIS37 pKa, and HIS100 pKa. Extend the table only when the reporting specification is updated; do not add ad hoc rows to individual slides.

## Common failure modes

- Missing one species’ metric or render.
- Using directory names instead of actual sequence-based design IDs.
- Passing display-remapped chain IDs to DeltaForge.
- Stretching one molecular panel to fit the slide.
- Manually typing values that already exist in normalized outputs.
- Forgetting to regenerate PNG summary assets from changed SVG inputs.

## Final inspection

- [ ] PPTX opens without repair warnings.
- [ ] Slide count is correct.
- [ ] Both species are present on every design slide.
- [ ] Human is `lightorange`; Mouse is `skyblue`; DD1 is `cyan`; DD2, when present, is `magenta`.
- [ ] HIS/NEG/POS/NAG are purple/blue/red/yellow sticks with standard non-carbon element colors; hotspot precedence matches `PYMOL_VISUALIZATION.md`.
- [ ] Sequence and score values match source outputs.
- [ ] No images, tables, or captions overlap.
- [ ] The deck and its source assets are stored together.
