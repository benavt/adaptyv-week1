# Uniform single-slide summaries for EGFR–DD1 complexes

## Purpose

Use this specification for every one-slide summary of a DD1 design or mutant evaluated against Human and Mouse EGFR D3. Each slide should make the structural model and the corresponding scores readable at a glance and should be reproducible from the transferred Boltz2 outputs.

## Required slide layout

Use a 16:9 slide with a light background.

1. Title at the top: `DD1_<design or mutant> | Human and Mouse EGFR D3`.
2. Two structure panels directly below the title:
   - left: Human EGFR D3 complex
   - right: Mouse EGFR D3 complex
3. A compact score table beneath the structures, with columns:
   - Metric
   - Human EGFR D3
   - Mouse EGFR D3
4. The full single-letter DD1 binder sequence directly above the color key and score table.
5. The color key below the sequence and above the score table.
6. A short method note below the table.
7. Keep captions and color keys outside the table boundaries. Do not place explanatory text over the molecular images or table.

The Human and Mouse PyMOL renderings must use identical displayed dimensions and must preserve the source image height-to-width ratio. Resize by setting only one dimension and deriving the other from the source ratio, or use a crop/fit operation that does not stretch the image. Never independently assign width and height values that distort the molecular view.

The structure panels are required. A summary slide is incomplete if it contains only the score table.

The sequence line should identify the design or mutant and show the complete DD1 sequence in a left-aligned, readable text box. Wrap the sequence only if necessary, keeping both the sequence and color key above the score table.

## Structure source and chain mapping

Use the model-0 CIF transferred from the Boltz2 prediction directory:

```text
01_Staging/DD1_<design_or_mutant>_EGFR_species_intersection/
  predictions/DD1_<design_or_mutant>_Human_EGFR/*_model_0.cif
  predictions/DD1_<design_or_mutant>_Mouse_EGFR/*_model_0.cif
```

Boltz2 chain convention:

- chain A: DD1 binder
- chain B: EGFR D3

Do not silently change this mapping when extracting metrics or rendering structures. DeltaForge should score receptor chain B and peptide chain A.

## PyMOL rendering convention

Render one PNG per species and place them above the table. Use the fixed views already captured for the project:

```text
03_Filtering/PyMol_top20/view_capture/human_view.json
03_Filtering/PyMol_top20/view_capture/mouse_view.json
```

Before applying the fixed view, align Boltz2 chain B to the matching EGFR medoid in `01_Staging/EGFR/`:

- Human reference: `Human_EGFR_seg_medoid.pdb`
- Mouse reference: `Mouse_EGFR_seg_medoid.pdb`

Place the two resulting PNGs at the same width and height on the slide. If the sequence and color key require more vertical space, move the score table downward and reduce unused margins rather than shrinking or stretching one molecular panel independently.

Use these colors and representations:

| Object | Selection | Display | Color |
|---|---|---|---|
| Human EGFR D3 | EGFR chain B | cartoon | `lightorange` |
| Mouse EGFR D3 | EGFR chain B | cartoon | `skyblue` |
| DD1 | DD1 chain A | cartoon | `cyan` |
| DD2 (when present) | recorded DD2 chain | cartoon | `magenta` |
| HIS | all molecules, `resn HIS` | sticks | `purple` |
| NEG | all molecules, `resn ASP+GLU` | sticks | `blue` |
| POS | all molecules, `resn ARG+LYS` | sticks | `red` |
| NAG (when present) | `resn NAG` | sticks | `yellow` |
| Original EGFR hotspots | EGFR residues `35+37+67+96+100+122+125` | sticks | `green` except residue-class overrides |

Use the shared `03_Filtering/PyMol_top20/molecular_colors.py` implementation and `docs/PYMOL_VISUALIZATION.md`. These rules reproduce the original EGFR hotspot script. Apply residue colors after hotspot colors, with `util.cnc` for standard non-carbon element colors. DD1 ASP is blue NEG sticks; EGFR HIS37/HIS100 are purple HIS sticks.

Use a white background, orthoscopic projection, and the fixed species-specific view. Colors follow molecule identity: EGFR is not DD2. Set cartoon colors independently of atom/stick colors and apply them before ray tracing or PNG export. Save the PML and rendered PNG beside the report or in a run-specific rendering directory.

## Required score table

Report the following rows whenever the source output provides them:

- ipSAE
- ipTM
- complex pLDDT
- DeltaForge ΔG in kcal/mol
- DeltaForge Kd in nM
- HIS37 pKa
- HIS100 pKa
- Mouse-versus-Human whole-complex TM-score and RMSD

Use `N/A` when a metric is not present in the transferred output. Do not infer ipSAE from a different model or silently substitute another interface score.

For pKa values, retain three decimals in the CSV and three decimals or two decimals on the slide. Use the project heatmap scale of 6.5 to 7.4 when color-coding pKa cells. Flag values outside that range in red or with an `out of range` note.

For inter-model comparison, run USalign on the Mouse and Human model-0 structures and report the TM-score normalized by the reference structure together with RMSD in Å. State which structure was used as the reference.

## Recommended output layout

For each design or mutant, create a unique directory under `03_Filtering`, for example:

```text
03_Filtering/DeltaForge/out/DD1_062_Y32D/
  Human.json
  Mouse.json
  results.csv
  summary.csv
  DD1_062_Y32D_summary.pptx
  renders/
    DD1_062_Y32D_Human.png
    DD1_062_Y32D_Mouse.png
  pml/
    DD1_062_Y32D_Human.pml
    DD1_062_Y32D_Mouse.pml
```

Keep JustHISpKa outputs separately under:

```text
03_Filtering/JustHIpKA/DD1_062_Y32D/
```

The slide builder should read the score CSVs and the rendered PNGs rather than embedding ad hoc values or manually edited screenshots.

## Reproducibility checklist

Before delivering a slide, verify:

- Human and Mouse model-0 CIFs are both present.
- Both structure panels are visible and use the fixed views.
- Human EGFR cartoon is lightorange; Mouse EGFR cartoon is skyblue; DD1 is cyan; DD2 is magenta when present.
- All present HIS/NEG/POS/NAG use purple/blue/red/yellow sticks with standard non-carbon element colors; cartoon colors remain independent.
- EGFR HIS37/HIS100 and original EGFR hotspots are shown as sticks with the documented colors.
- Chain A/B mapping is consistent across PyMOL, DeltaForge, and metric extraction.
- The table contains the same metric rows and units across runs.
- Missing metrics are labeled `N/A`.
- The inter-model TM-score and RMSD specify the comparison and reference orientation.
- No caption, legend, or color key overlaps the table or structure panels.
- The slide opens with exactly one summary slide for the requested design or mutant.

## Implementation note

Existing summary builders under `03_Filtering/DeltaForge/` can be used as data and table templates, but they must be extended to generate and insert the Human and Mouse PyMOL PNGs before the score table. Use the established rendering logic in `03_Filtering/PyMol_top20/` as the reference implementation for alignment, views, and color conventions.
