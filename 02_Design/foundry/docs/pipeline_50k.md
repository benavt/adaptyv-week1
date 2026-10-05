# 50k pipeline (RFD3 → RF3, no MPNN)

[`pipeline_50k.slurm`](../pipeline_50k.slurm) is the large-library wrapper: **50k RFD3 designs, skip MPNN, fold those RFD3 CIFs with RF3, then USalign**. One GPU allocation. Nested scripts run as bash, not as extra Slurm jobs.

Default walltime is **14 days**. Override with `sbatch --time`.

This is the same machinery as [the design pipeline](pipeline.md) except MPNN is skipped and RF3’s input directory is `runs/<id>/rfd3/outputs` instead of MPNN outputs.

## Submit

Required when starting at RFD3 (the default): `YAML` and `PDB`. Same pattern as the other pipeline:

```bash
sbatch --export=ALL,YAML=runs/Brd4ET_hot_4_50000/Brd4ET.yaml,PDB=runs/Brd4ET_hot_4_50000/Brd4ET_closed_state_fix.pdb,RUN_ID=Brd4ET_hot_4_50000 pipeline_50k.slurm
```

Do **not** point `RUN_ID` at a finished small pipeline folder (`Brd4ET_hot_4`, isoform fingerprint without `_50000`, and so on). Stage a new directory with YAML + PDB first.

## What is different from `pipeline.slurm`

| Topic | `pipeline.slurm` | `pipeline_50k.slurm` |
|---|---|---|
| Default design count | `N_DESIGNS=100` (13 batches / 104 CIFs) | `N_DESIGNS=50000` (6250 batches / 50000 CIFs) |
| MPNN | SolubleMPNN, `N_SEQS=10`, interface skip on by default | **Skipped.** Log line: `skipping MPNN (pipeline_50k folds RFD3 CIFs directly)` |
| What RF3 folds | `runs/<id>/mpnn/outputs` (or `MPNN_OUTPUTS`) | `runs/<id>/rfd3/outputs` |
| Walltime | 48 hours | 14 days |
| Archived script name | `pipeline.slurm` | `pipeline_50k.slurm` |
| `START_STAGE` you would actually use | `rfd3`, `mpnn`, `rf3`, `tm_align` | `rfd3`, `rf3`, or `tm_align` (`mpnn` is accepted by the parser but there is no MPNN stage) |

RFD3 still writes 8 CIFs per batch, so `N_DESIGNS` is rounded up to a multiple of 8 internally via `N_BATCHES = ceil(N_DESIGNS / 8)`.

## Knobs (same names as the design pipeline)

| What you set | Env var | Default here | Meaning |
|---|---|---|---|
| Design spec / target | `YAML`, `PDB` | required at `START_STAGE=rfd3` | See [RFD3](rfd3.md). Binder length, hotspots, loopy vs not, A30–A50 `BKBN` vs `[]` all live in the YAML, not in these env vars. |
| First stage | `START_STAGE` | `rfd3` | `rfd3`, `rf3`, or `tm_align` in practice. |
| Existing RFD3 CIFs | `RFD3_OUTPUTS` | unset | Snapshot into `runs/<id>/rfd3/outputs/` when that folder is empty (mid-start). |
| How many backbones | `N_DESIGNS` | `50000` | Target CIF count (8 per batch). |
| Folder name | `RUN_ID` | `<jobid>_pipe` | Exact name under `runs/`. Nested `rfd3/` and `rf3/` only (no `mpnn/`). |
| RFD3 / RF3 weights | `RFD3_CKPT`, `RF3_CKPT` | same cluster defaults as the other pipeline | Passed through as `CKPT` to each stage. |
| RF3 diffusion steps | `NUM_STEPS` | `50` | Same speed tradeoff as the design pipeline. |
| RF3 samples per input | `DIFFUSION_BATCH_SIZE` | `1` | One fold per RFD3 CIF. |
| USalign / archive | `USALIGN`, `ARCHIVE_DIR`, `RSYNC_REFRESH_SEC` | same as [design pipeline](pipeline.md) | Linux USalign is still required up front. |

`N_SEQS`, `DESIGNED_CHAINS`, `EXCLUDE_INTERFACE`, `INTERFACE_CUTOFF`, and `MPNN_CKPT` are still given defaults in the script and written to `manifest.json`, but **MPNN never runs**, so they do not change the science.

RF3 extras (`N_RECYCLES`, `SEED`, `EARLY_STOPPING_PLDDT`, `SKIP_EXISTING`, `RF3_EXTRA`) still leak through the environment into `rf3.slurm` if you export them. See [RF3](rf3.md).

## Mid-pipeline restart

To fold already-finished 50k RFD3 CIFs:

```bash
sbatch --export=ALL,START_STAGE=rf3,RFD3_OUTPUTS=runs/some_rfd3/outputs,RUN_ID=some_rfd3_rf3_50k pipeline_50k.slurm
```

If `runs/<RUN_ID>/rfd3/outputs` already has CIFs, you can omit `RFD3_OUTPUTS`.

## What you get

- `runs/<id>/rfd3/outputs/` — 50k (or `N_DESIGNS`) RFD3 CIFs
- `runs/<id>/rf3/outputs/` — RF3 folds of those CIFs (same sequences as RFD3; no MPNN redesign)
- `runs/<id>/scores.csv` — TM / RMSD of ranked RF3 `*_model.cif` vs the parent RFD3 CIF, plus RF3 confidence columns

Scoring is the same pairing rule as the design pipeline (strip `_bN_dN` from the RF3 stem). See [TM-align](tm_align.md).
