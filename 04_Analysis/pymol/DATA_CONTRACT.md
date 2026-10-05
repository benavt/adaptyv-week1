# Data contract: master_alignment.csv

The snippet bundle uses **`master_alignment.csv`** as its primary input. This file is produced by the CSP_UBQ pipeline (`scripts/merge_csv.py`) and already aligns CSP sequential indices with PDB residue numbers.

## Required columns

### All three panels

| Column | Type | Usage |
|--------|------|-------|
| `holo_resi` | int | PyMOL `resi` for receptor coloring |
| `chain` | str | Receptor chain ID (e.g. `A`) |

Optional override: pass `--receptor-chain` / `--ligand-chain` on the CLI (ligand defaults to `B`).

### CSP mask panel (`color_csp_mask.pml`)

| Column | Type | Usage |
|--------|------|-------|
| `significant` | 0/1 or bool | Residue has significant combined N/H CSP |

Rows with `significant == 1` are colored **red** on the receptor chain.

### Binding-site panel (`color_occlusion.pml`)

A residue is in the binding-site union if **any** of the following is true (matches pipeline ground truth):

| Column | Meaning |
|--------|---------|
| `is_occluded_occlusion` | SASA backbone occlusion |
| `passes_filter_distance` | Within CA-distance threshold |
| `has_hbond_interaction` | H-bond to ligand |
| `has_charge_complement_interaction` | Charge complementarity |
| `has_pi_contact_interaction` | Pi contact |
| `passes_sub_2A_filter_any_atom` | Sub-2 Å inter-chain contact flag |
| `min_any_atom_distance` | Used as fallback: `< 2.0` Å counts as binding |

Binding-site residues are colored **red**; other receptor residues **gray30**; ligand **cyan**.

### Classification panel (`csp_classification_original.pml`)

| Column | Type | Usage |
|--------|------|-------|
| `classification` | TP / FP / TN / FN | Precomputed label per residue |

If `classification` is missing, the snippet derives it from `significant` + binding-site flags (same logic as `scripts/merge_csv.py`).

## Example row

From `outputs/1DDM_4683/master_alignment.csv`:

```
holo_resi=3, significant=1, classification=FP,
is_occluded_occlusion=false, passes_filter_distance=false,
has_hbond_interaction=false, ...
```

## PDB file expectations

| Panel | Typical pipeline PDB | Snippet argument |
|-------|---------------------|------------------|
| CSP mask | `{pdb}_csp.pdb` | `--structure-pdb` |
| Binding site | `{pdb}_delta_sasa.pdb` | `--occlusion-pdb` |
| Classification | `{pdb}_csp.pdb` | `--structure-pdb` |

You may pass the **same PDB** for both structure and occlusion if a delta-SASA PDB is unavailable; coloring is by residue number, not B-factor.

PDB requirements:

- Cartoon-compatible (standard ATOM records)
- Receptor and ligand on distinct chains matching CSV / CLI chain IDs
- `holo_resi` values must exist as `resi` in the receptor chain

## Classification colors

| Label | Hex | PyMOL custom color |
|-------|-----|-------------------|
| TP | `#2ecc71` | `tp_color` |
| FP | `#9b59b6` | `fp_color` |
| TN | `#3498db` | `tn_color` |
| FN | `#f39c12` | `fn_color` |

Ligand chain: **cyan**. Receptor without CSP/classification data: **gray**.

## View JSON format

Saved camera files (e.g. `pymol_views/1CF4_18251_case_study_view.json`) are JSON arrays of **18 floats** — the PyMOL `get_view()` tuple. Produced by interactive F5 capture in the pipeline or `--interactive-view` in this bundle.

## Minimal CSV for a new project

If building CSV from scratch, each row needs at least:

```csv
holo_resi,chain,significant,classification,is_occluded_occlusion,passes_filter_distance,has_hbond_interaction,has_charge_complement_interaction,has_pi_contact_interaction,passes_sub_2A_filter_any_atom,min_any_atom_distance
42,A,1,TP,true,false,false,false,false,false,
57,A,0,TN,false,false,false,false,false,false,12.5
```

Only include rows you want considered; prolines / gaps are typically omitted in pipeline output.
