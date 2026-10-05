# ProtonPottsMPNN: DD2 protonation smoke-test plan

Prepared 2026-10-02. Submodule: `03_Filtering/ProtonPottsMPNN`, upstream
[christian-creator/ProtonPottsMPNN](https://github.com/christian-creator/ProtonPottsMPNN),
pinned at `09682abfa7d20e0abcdeea0490b7a4b1c190aee3`.

## Scope and interpretation

Label the existing `top_02`, `top_07`, and `top_11` DD2 seed binders in their
staged Human EGFR complexes. These are the Boltz2 model-0 structures staged under
Foundry's DD2 hotspot-display runs, rather than newly generated RFD3 outputs.
Chain A is Human EGFR D3 (192 residues); chain B is the binder.

Use the deployed **EV6 FLAML ensemble** to predict HIS/ASP/GLU labels and retain
`p_protonated`, `sd`, `ambiguous`, and `metal_adjacent`. Cross-check the discrete
labels through `prepare_potts_input(..., extended_vocab="v6")` on the first
complex. This labeling task requires HBPLUS and the Python dependencies, but
does not require a Potts checkpoint forward pass, RF3 weights, a GPU, or redesign.

EV6 has **no pH argument**. Its scores describe the learned structural classifier;
they are not calculated equilibrium occupancy at pH 5.0 or 7.4 and must not be
converted directly into pKa. If pH-specific states are needed, follow this smoke
test with a separate pKa calculation (e.g. the project's JustHISpKa workflow for
His, plus a suitable acid pKa method). Use the resulting pKa values to estimate
protonated fractions at the chosen assay pH values. Those estimates and the EV6
scores must remain separately identified.

## Verified input inventory

All paths below are relative to the workspace root. Preserve chain IDs, residue
numbers, insertion codes, sequence, and heavy-atom coordinates.

| Binder | Complex PDB under `02_Design/foundry/runs/` | Binder length | Binder HIS | Binder ASP | Binder GLU | Bound / binder-only titratable rows |
| --- | --- | ---: | --- | --- | --- | ---: |
| top_02 | `DD2_top02_Human_EGFR_hotspot_display/DD2_top02_Human_EGFR_complex.pdb` | 89 | none | 30, 55, 66 | 70 | 29 / 4 |
| top_07 | `DD2_top07_Human_EGFR_hotspot_display/DD2_top07_Human_EGFR_complex.pdb` | 107 | 59, 98 | 29, 40, 81, 82, 90 | 28, 67, 91, 102, 103 | 37 / 12 |
| top_11 | `DD2_top11_Human_EGFR_hotspot_display/DD2_top11_Human_EGFR_complex.pdb` | 102 | 54, 70 | 11, 15, 21, 25, 28, 47, 83 | 73, 89, 90, 100, 101 | 39 / 14 |

Each target chain has 25 titratable residues: six His, nine Asp, ten Glu.
The inventory is an input check, **not predicted protonation results**.
`top_02` is a useful acid-only binder case; target histidines still exercise the
His classifier in its complex. There are 30 binder titratable sites total.

Source SHA-256 values:

```text
top_02  5a2d20a5fed622f760dd04cef2d2b19274cc206ecaa628ec9222d9878ebd242d
top_07  06439e92566e6ac71f38471415d4616ff15ff8c41bb237228e69f655d74f3a8b
top_11  c0643e9279bc5b661f7a02fd37a3c9c7bc4b3edeac8fed2757588c77fcd8c53f
```

Derive every binder-only control from chain B of its corresponding complex,
including `top_02`, which has no pre-existing binder-only PDB in its display run.
This isolates the effect of removing EGFR at fixed binder coordinates. Call the
control `binder_only_fixed_coordinates`; it is not a relaxed apo structure.

## Documentation and source review

Reviewed the upstream README, `install.sh`, `requirements-extra.txt`,
`labeller/label_pdb.py`, and the bundled `mpnn` implementation:

- `transforms/ev6/predictor.py`: five-fold ensemble mean and population SD;
  deterministic per-residue calls.
- `transforms/ev6/weights/thresholds.json`: default His probability threshold
  0.40 / SD cutoff 0.20; acid probability threshold 0.06 / SD cutoff 0.30.
- `transforms/extended_vocab_v6.py`: the pipeline can override probability
  thresholds with `EV6_HIS_PROB_THR` / `EV6_ACID_PROB_THR`.
- `transforms/ev6/features.py`: heavy-atom geometry, SASA, and HBPLUS features;
  H···A cutoff 3.2 Å and donor–acceptor cutoff 4.0 Å. HBPLUS pools must use these
  settings. SASA failure can silently substitute zeros, so inspect warnings and
  explicitly verify SASA works during preflight.
- `potts_inference.py`: `prepare_potts_input` defaults to vocabulary v4;
  explicitly request **v6**, deterministic labeling, and zero structure noise.

At the default operating point, a metal-adjacent site or SD **greater than** its
cutoff produces `-A`; otherwise mean probability **at least** its threshold
produces `-P`. Other sites are HIS-S, ASP-D, or GLU-D. HIS-S is neutral and does
not distinguish HID/HIE tautomers; ambiguous calls must remain ambiguous.

The default cutoffs are prevalence-matched to the upstream neutron training
set, rather than calibrated for these predicted EGFR complexes. Report raw
scores alongside tokens; do not assume that the acid threshold 0.06 implies
most sites called `-P` have a probability above 0.5. The v6 Potts design
checkpoint's directory uses `his0.3_acid0.06`; that training setting is separate
from the EV6 default His threshold 0.40. Use the documented EV6 defaults for this
labeling smoke test and record them explicitly.

## Planned execution

1. **Environment preflight.** Use a dedicated Python 3.12 environment for this
   submodule's bundled Foundry, avoiding import collisions with `02_Design/foundry`.
   Verify `mpnn.__file__` resolves inside the new submodule. Confirm HBPLUS is an
   executable compatible with the execution host; use an explicit absolute
   `HBPLUS_PATH` to avoid the author's hard-coded fallback. Record the submodule
   commit, Python/platform information, package versions, HBPLUS path/version or
   binary hash, input hashes, and thresholds.
2. **Structure preflight.** Check one model, A/B chains, one CA per residue,
   canonical residue names, finite coordinates, unique atom/residue keys, absence
   of alternate conformers and insertion-code collisions, and complete backbone
   and titratable side-chain heavy atoms. Require HIS CG/ND1/CD2/CE1/NE2,
   ASP CG/OD1/OD2, and GLU CD/OE1/OE2. Reject missing atoms rather than silently
   impute them. Confirm the inventory above and that SASA computation succeeds.
3. **First inference: top_07 bound.** This exercises both binder histidines and
   acids. Load the complex using `MPNNInferenceInput`, remove H/D if present,
   and call `EV6Predictor().predict(atom_array)`. Save all 37 target/binder rows
   and the 12 chain-B rows. Repeat once and check identical keys/tokens plus
   scores/SD within `1e-6`. Check v6 pipeline token parity by residue key.
4. **Complete the matrix.** Label the other two complexes, then label chain-B
   subsets of all three complexes separately: six condition runs in total,
   plus the first-case repeat and pipeline cross-check. Preserve the original
   PDBs and keep generated outputs outside the upstream submodule.
5. **Compare.** Join bound and binder-only rows on binder/residue identity;
   report tokens, raw probability/SD, and `bound minus binder_only` probability.
   Flag state changes and ambiguous calls for review, rather than treating them
   as a measured binding or pH-switch effect.

The configured local environment is `03_Filtering/ProtonPottsMPNN/.venv` and
uses the compiled HBPLUS executable at
`03_Filtering/ProtonPottsMPNN_runtime/build/hbplus/hbplus`. The upstream
requirements leave XGBoost unbounded; this runtime pins `xgboost==2.1.4`
because the shipped serialized EV6 acid models are incompatible with the
newly selected 3.x API.

To recreate the Python environment from the workspace root:

```bash
python3.12 -m venv 03_Filtering/ProtonPottsMPNN/.venv
03_Filtering/ProtonPottsMPNN/.venv/bin/python -m pip install \
  -e ./03_Filtering/ProtonPottsMPNN/foundry \
  -r ./03_Filtering/ProtonPottsMPNN/requirements-extra.txt
03_Filtering/ProtonPottsMPNN/.venv/bin/python -m pip install --force-reinstall \
  'xgboost==2.1.4'
```

Use a new environment directory. The upstream `install.sh` uses `--clear` and
registers a user Jupyter kernel; neither operation is needed for this smoke test.
Before inference, unset `EV6_HIS_PROB_THR`, `EV6_ACID_PROB_THR`, and stray `DEBUG`;
set `CUDA_VISIBLE_DEVICES=""` and the compiled `HBPLUS_PATH` in the smoke-test
process. On this macOS sandbox, Matplotlib's font discovery can block on
`system_profiler`; set a writable `MPLCONFIGDIR` with a prebuilt font cache when
running the model process.

The APIs for the runner are:

```python
from mpnn.utils.inference import MPNNInferenceInput
from mpnn.transforms.ev6 import EV6Predictor
from mpnn.potts_inference import prepare_potts_input

inf = MPNNInferenceInput.from_atom_array_and_dict(
    input_dict={"structure_path": str(complex_pdb), "structure_noise": 0.0}
)
aa = inf.atom_array
aa = aa[(aa.element != "H") & (aa.element != "D")]
predictor = EV6Predictor()
bound = predictor.predict(aa)
binder_only = predictor.predict(aa[aa.chain_id == "B"].copy())

# First-case parity check: compare CA labels to bound.token by residue identity.
pipeline = prepare_potts_input(
    aa.copy(), extended_vocab="v6", structure_noise=0.0,
    deterministic=True, protonation_seed=0,
)
ca = pipeline["atom_array"]
ca = ca[ca.atom_name == "CA"]
tokens = ca.get_annotation("protonation_label")
```

The Python imports, EV6 feature generation, and HBPLUS subprocess path have
been validated locally. Full EV6 model prediction still exits in a native
scikit-learn/XGBoost path on this macOS arm64 host, so the six-row output matrix
remains pending a compatible runtime or Linux execution host. Do not treat the
partial setup as protonation results.

## Outputs and acceptance criteria

Write outputs under `03_Filtering/ProtonPottsMPNN_smoke/DD2/`:

- `manifest.json`: provenance, environment, thresholds, condition definitions,
  and execution status/timing.
- `<binder>/<condition>/protonation_labels.csv`: all HIS/ASP/GLU rows, with
  binder ID, condition, chain, residue ID/name, token, `p_protonated`, `sd`,
  `ambiguous`, and `metal_adjacent`. Record insertion codes if supported; otherwise
  require that input insertion codes are blank during preflight.
- `binder_comparison.csv`: 30 matched binder sites with both conditions' calls,
  probabilities/SD, and probability differences.
- `summary.md` and per-run logs: count calls by residue type/token, highlight
  binder H59/H98 and H54/H70, and list unresolved warnings.

Pass when all six conditions finish without missing HBPLUS/model/import errors;
bound row counts are 29/37/39 and binder-only counts are 4/12/14; binder residue
keys match exactly between conditions; probabilities are finite in [0, 1]; SDs
are finite and nonnegative; every token obeys the recorded threshold rule; the
top_07 repeat is reproducible; and its pipeline tokens agree with direct EV6
calls. A legitimate ambiguous prediction or lack of protonated sites is not a
failure. Do not prescribe a particular biological token as an expected result
before inference.

## Status at handoff

Submodule registration and checkout verified; upstream checkout is clean.
Source hashes, CA-based lengths, and HIS/ASP/GLU inventory verified locally.
The Python 3.12 environment and HBPLUS v3.06 build are configured. HBPLUS
successfully processed the `top_07` complex and reported 253 hydrogen bonds.
EV6 feature generation succeeded for 8 His and 29 Asp/Glu rows, but the full
model call still has a native macOS arm64 compatibility issue as described above.
