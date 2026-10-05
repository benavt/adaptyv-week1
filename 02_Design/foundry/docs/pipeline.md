# Design pipeline

[`pipeline.slurm`](../pipeline.slurm) runs **RFD3 → SolubleMPNN (non-interface residues by default) → RF3 → USalign** in **one GPU allocation**. It invokes the existing Slurm scripts as bash (they do not queue separate jobs).

Default walltime is **48 hours**. Override with `sbatch --time`.

Outputs live under `runs/<id>/{rfd3,mpnn,rf3,tm_align}/` plus `runs/<id>/scores.csv`. Nested stage folders reuse `RUN_ID=<id>/rfd3` (and `mpnn`, `rf3`).

A Linux USalign binary is required **before** GPU work starts; see [TM-align](tm_align.md).

## Submit

Required when starting at RFD3 (the default): `YAML` and `PDB`.

```bash
sbatch --export=ALL,YAML=runs/Brd4ET_design_inputs/Brd4ET.yaml,PDB=runs/Brd4ET_design_inputs/Brd4ET_closed_state_fix.pdb,RUN_ID=Brd4ET_design pipeline.slurm
```

## Which stage to start at

| What you set | Env var | Default | Meaning |
|---|---|---|---|
| First stage | `START_STAGE` | `rfd3` | `rfd3`, `mpnn`, `rf3`, or `tm_align`. Earlier stages are skipped. |
| Existing RFD3 CIFs | `RFD3_OUTPUTS` | unset | Directory of `*.cif` / `*.cif.gz` to copy into `runs/<id>/rfd3/outputs/` when that folder is empty. Needed if you start at `mpnn` (or later) without leftover CIFs in the run. The source directory is not modified; copies are then `gzip -dkf` extracted. |
| Existing MPNN CIFs | `MPNN_OUTPUTS` | `runs/<id>/mpnn/outputs` | CIF directory RF3 should fold when you start at `rf3` (or if you want a different MPNN folder). |

YAML and PDB are **not** required unless `START_STAGE=rfd3`. If you pass them anyway on a mid-start, they are recorded in `manifest.json`.

Use a **new** `RUN_ID` when snapshotting an existing gz-only RFD3 folder so the original stays gz-only.

```bash
sbatch --export=ALL,START_STAGE=mpnn,RFD3_OUTPUTS=runs/NS1BctdRNA/outputs,RUN_ID=NS1BctdRNA_pipe,N_SEQS=10,DESIGNED_CHAINS=A pipeline.slurm
```

## How many designs and sequences

| What you set | Env var | Default | Meaning |
|---|---|---|---|
| How many RFD3 backbones | `N_DESIGNS` | `100` | Target design count. RFD3 writes **8 CIFs per batch**, so the job asks for `ceil(N_DESIGNS / 8)` batches. **100 → 13 batches → 104 designs.** |
| Sequences per backbone | `N_SEQS` | `10` | MPNN redesigns per RFD3 CIF. 104 designs × 10 ≈ 1040 RF3 folds. |

`N_BATCHES` is computed from `N_DESIGNS`; do not expect to set `N_BATCHES` yourself on this wrapper (it is passed into `rfd3.slurm` and unset for MPNN).

## MPNN behavior in this wrapper

| What you set | Env var | Default | Meaning |
|---|---|---|---|
| Chains to redesign | `DESIGNED_CHAINS` | `A` | Binder chain. Target stays fixed. |
| Skip interface residues | `EXCLUDE_INTERFACE` | `True` | Redesign only non-interface residues on those chains. Set `False` to sequence **all** of chain A. |
| Interface distance | `INTERFACE_CUTOFF` | `5.0` | Ångströms; see [MPNN](mpnn.md). |

This wrapper **hardcodes** SolubleMPNN:

- `MPNN_CKPT` default `solublempnn_v_48_020.pt` (repo-local)
- `MODEL_TYPE=protein_mpnn`
- `IS_LEGACY_WEIGHTS=True`

You can change the checkpoint with `MPNN_CKPT`. You cannot change `MODEL_TYPE` / `IS_LEGACY_WEIGHTS` via export; the wrapper always sets those two.

## Checkpoints and RF3 speed

| What you set | Env var | Default | Meaning |
|---|---|---|---|
| RFD3 weights | `RFD3_CKPT` | `/opt/databases/foundry/rfd3_latest.ckpt` | Passed to `rfd3.slurm` as `CKPT`. |
| MPNN weights | `MPNN_CKPT` | `<repo>/solublempnn_v_48_020.pt` | Passed to `mpnn.slurm` as `CKPT`. Relative paths are from the repo root. |
| RF3 weights | `RF3_CKPT` | `/opt/databases/foundry/rf3_foundry_01_24_latest_remapped.ckpt` | Passed to `rf3.slurm` as `CKPT`. |
| RF3 diffusion steps | `NUM_STEPS` | `50` | Faster than RF3’s own default of 200. |
| RF3 samples per input | `DIFFUSION_BATCH_SIZE` | `1` | One folded model per MPNN CIF, not RF3’s default of 5. |

RF3 extras that are **not** listed above still work if you export them: `N_RECYCLES`, `SEED`, `EARLY_STOPPING_PLDDT`, `SKIP_EXISTING`, `RF3_EXTRA`. The wrapper does not unset those, so `rf3.slurm` will see them. Details in [RF3](rf3.md).

## Names, archive, USalign

| What you set | Env var | Default | Meaning |
|---|---|---|---|
| Pipeline folder | `RUN_ID` | `<jobid>_pipe` | Exact name under `runs/`. Nested stages are `<id>/rfd3`, `<id>/mpnn`, `<id>/rf3`. |
| YAML / PDB | `YAML`, `PDB` | required at `START_STAGE=rfd3` | Same meaning as [RFD3](rfd3.md). The YAML `input:` field is rewritten inside `runs/<id>/rfd3/`. |
| USalign binary | `USALIGN` | Linux `USalign` on `PATH` | Resolved before any GPU stage. |
| Archive parent | `ARCHIVE_DIR` | `<repo>/archive` | Final copy of `slurm.out`, `manifest.json`, `scores.csv`, `README.md`, and `pipeline.slurm`. Nested stages also archive themselves while they run. |
| Archive refresh | `RSYNC_REFRESH_SEC` | `60` | Honored by the nested `rfd3` / `mpnn` / `rf3` scripts. |

If `runs/<id>/README.md` exists at start, it is copied into `runs/<id>/rfd3/README.md`.

## What you get

- `runs/<id>/rfd3/outputs/` — RFD3 CIFs
- `runs/<id>/mpnn/outputs/` — sequenced CIFs + FASTA
- `runs/<id>/rf3/outputs/` — RF3 models and confidences
- `runs/<id>/scores.csv` — TM / RMSD of ranked RF3 `*_model.cif` vs the parent RFD3 CIF, plus RF3 summary scalars (`overall_plddt`, `ptm`, `iptm`, …)

See [TM-align](tm_align.md) for column meanings.
