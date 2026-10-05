# PyMOL visualization conventions

## Source and molecule identity

The EGFR and residue conventions come from `01_Staging/EGFR/EGFR_seg_medoid_hotspots.pml` and `01_Staging/EGFR/render_egfr_targeting_pptx.py`. The updated design-backbone convention is DD1 cyan and DD2 magenta.

Colors follow molecule identity. In the canonical Boltz2 DD1–EGFR complex, chain A is DD1 and chain B is EGFR D3. Chain B is **not DD2**. In DD1–DD2 complexes, assign DD1 and DD2 using the recorded chain mapping. Display-only remapped structures must state their actual roles; the DD1_062 hotspot display uses Human EGFR chain A and DD1 chain B. Never change scoring inputs to suit a display convention.

## Colors and representations

| Component | Representation | Color |
|---|---|---|
| Human EGFR D3 | cartoon | `lightorange` |
| Mouse EGFR D3 | cartoon | `skyblue` |
| DD1 | cartoon | `cyan` |
| DD2 | cartoon | `magenta` |
| HIS (`resn HIS`, including EGFR HIS37/HIS100) | sticks | `purple` |
| NEG (`resn ASP+GLU`) | sticks | `blue` |
| POS (`resn ARG+LYS`) | sticks | `red` |
| NAG (`resn NAG`, when present) | sticks | `yellow` |
| Original EGFR hotspots | sticks | `green` |
| Editable interface hotspots | sticks | `tv_orange` |

Apply the residue classes to all molecules present, including DD1 and DD2. Their stick carbon colors take precedence over hotspot/mismatch highlights. Use `util.cnc` to retain standard non-carbon element colors, as in the source EGFR script. Histidines are purple, acidic residues are blue, and ARG/LYS are red. The old filtering-specific overrides (blue HIS, red EGFR acidic residues, cyan DD1 ASP) are superseded.

Set **cartoon-specific** colors independently of atom colors: residue highlighting must never change a molecule's backbone color. NAG is displayed only when it exists in the source coordinates; do not invent glycans. Remove solvent and hide unrelated reference objects. Record any other scientific ligand and its display convention separately.

## Shared filtering implementation

`03_Filtering/PyMol_top20/molecular_colors.py` generates the common PML block. Supply verified molecule-to-chain roles and EGFR species. Unknown species or overlapping chain roles must stop rendering. All filtering renderers use this block, and saved PMLs contain a self-contained copy so they remain replayable.

Apply the block after hotspot colors and **before `ray` or `png`**, so cached ray images cannot retain stale colors. Compact ipSAE color aliases also apply the block when invoked.

## Camera and rendering

- Use white background and orthoscopic projection for report figures.
- Align the EGFR chain to the matching Human/Mouse EGFR medoid before applying the fixed species view.
- Use `human_view.json` and `mouse_view.json` for standardized DD1 interface figures.
- Preserve matching panel dimensions and source aspect ratio.
- Save the generated PML beside the PNG and update embedded PowerPoint images when regenerating assets.

## Minimal command pattern

```text
load <complex.cif>, complex
hide everything
show cartoon, complex
set cartoon_color, cyan, complex and chain A
set cartoon_color, lightorange, complex and chain B  # skyblue for Mouse
select HIS, complex and resn HIS
select NEG, complex and resn ASP+GLU
select POS, complex and resn ARG+LYS
select NAG, complex and resn NAG
show sticks, HIS or NEG or POS or NAG
color purple, HIS
color blue, NEG
color red, POS
color yellow, NAG
util.cnc HIS
util.cnc NEG
util.cnc POS
util.cnc NAG
bg_color white
set orthoscopic, on
```

For a DD1–DD2 complex, use the recorded DD1 and DD2 chains and assign cyan and magenta cartoons respectively. Every figure must state any explicit deviation from these defaults.
