# PyMOL visualization snippets

Portable, well-commented code to reproduce **case-study PyMOL assets** from the CSP_UBQ pipeline:

| Pipeline output | Snippet output |
|-----------------|----------------|
| `color_csp_mask.pml` | Same (written to `--pml-dir`) |
| `color_occlusion.pml` | Same |
| `csp_classification_original.pml` | Same |
| `case_study_assets/case_study_color_csp_mask.png` | Same basename |
| `case_study_assets/case_study_color_occlusion.png` | Same |
| `case_study_assets/case_study_csp_classification_original.png` | Same |

Also supports `case_study2_*` naming via `--variant v2`.

## Prerequisites

- **PyMOL** on `PATH` as `pymol` (headless: `pymol -c -q`)
- Python 3.10+ (stdlib only; see [requirements.txt](requirements.txt))
- Input **`master_alignment.csv`** per target (see [DATA_CONTRACT.md](DATA_CONTRACT.md))
- Structure PDB files (typically `{pdb}_csp.pdb` and `{pdb}_delta_sasa.pdb` from pipeline outputs)

## Quick start

From this directory:

```bash
python examples/render_case_study_assets.py \
  --master-csv ../outputs/1CF4_18251/master_alignment.csv \
  --structure-pdb ../outputs/1CF4_18251/1cf4_csp.pdb \
  --occlusion-pdb ../outputs/1CF4_18251/1cf4_delta_sasa.pdb \
  --output-dir ../outputs/1CF4_18251/case_study_assets \
  --view-json ../pymol_views/1cf4_case_study_view.json
```

See [examples/sample_config.yaml](examples/sample_config.yaml) for path templates.

### Centroid gallery (BRD3ET / BRD4ET)

Loop each `centroid_*.pdb` under `BRD3ET/` and `BRD4ET/`, capture a camera in the PyMOL GUI (**F5**), then ray-trace a PNG beside the model:

```bash
python examples/render_centroid_views.py
```

Outputs next to each PDB: `centroid_CN.png` and `centroid_CN_view.json`. Re-render without reopening the GUI with `--reuse-views`; skip models that already have a PNG with `--skip-existing`. Requires a display (same F5 workflow as `--interactive-view` below).

## Three-step workflow

1. **Prepare inputs** — `master_alignment.csv` + holo structure PDB(s)
2. **Generate `.pml` scripts** — automatic inside the render CLI (or call `snippets.pymol_script_writers.write_all_case_study_scripts`)
3. **Render PNGs** — headless PyMOL ray-tracing with a saved or interactive camera view

## Module map

| File | Role |
|------|------|
| [snippets/colors.py](snippets/colors.py) | TP/FP/TN/FN hex colors |
| [snippets/binding_site.py](snippets/binding_site.py) | Binding-site union + classification logic |
| [snippets/csv_adapter.py](snippets/csv_adapter.py) | Load `master_alignment.csv` → residue sets |
| [snippets/pymol_script_writers.py](snippets/pymol_script_writers.py) | Write the three `.pml` files |
| [snippets/centroid_display.py](snippets/centroid_display.py) | Minimal cartoon `.pml` for centroid PDBs |
| [snippets/pymol_render.py](snippets/pymol_render.py) | Headless render + optional F5 view capture |
| [snippets/render_assets.py](snippets/render_assets.py) | End-to-end orchestrator |

## Camera / view

- **Reuse view:** pass `--view-json pymol_views/<slug>_case_study_view.json` (same format as pipeline)
- **Capture new view:** `--interactive-view` opens PyMOL GUI on `color_csp_mask.pml`; press **F5** to save and quit

## Differences from full pipeline

These snippets intentionally **do not** include:

- `CSPResult` / full CSP computation ([scripts/csp.py](../scripts/csp.py))
- Sequence alignment ([scripts/align.py](../scripts/align.py))
- PDB B-factor rewriting ([create_pdb_with_csp_bfactors](../scripts/visualize.py))
- Matplotlib case-study figure composition ([scripts/case_study.py](../scripts/case_study.py))

For full fidelity with live pipeline objects, use the original [scripts/visualize.py](../scripts/visualize.py) writers directly.

## Porting to another project

Copy the entire `pymol_visualization_snippets/` folder. Provide:

1. A CSV with the columns listed in [DATA_CONTRACT.md](DATA_CONTRACT.md)
2. A holo PDB with correct chain IDs
3. PyMOL installed

No dependency on the rest of CSP_UBQ is required.

## Original source references

| Snippet | Adapted from |
|---------|--------------|
| `pymol_script_writers.py` | `scripts/visualize.py` lines 700–800, 1521–1662, 2063–2273 |
| `pymol_render.py` | `scripts/case_study.py` lines 83–184 |
| `binding_site.py` | `scripts/merge_csv.py` lines 173–248 |
| `colors.py` | `scripts/config.py` lines 130–136 |
