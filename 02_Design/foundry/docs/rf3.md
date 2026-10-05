# RF3 (structure prediction)

[`rf3.slurm`](../rf3.slurm) runs Foundry `rf3 fold` on JSON / CIF / PDB inputs. One job call processes a file or a whole directory so the checkpoint is loaded once. Arguments are Hydra `arg=value`, not `--flags`.

Default walltime is **4 hours**. Override with `sbatch --time`.

## What you must pass

| What you set | Env var | Default | Meaning |
|---|---|---|---|
| What to fold | `INPUTS` | (required) | A `.json`, `.cif`, `.cif.gz`, `.pdb`, or `.pdb.gz`, **or** a directory of those. Sequences belong in JSON as `components[].seq`. **FASTA is not accepted**; in a mixed MPNN `outputs/` folder, sibling `.fa` files are ignored and the CIFs are folded. |

```bash
sbatch --export=ALL,INPUTS=runs/Brd4ET_2026-09-08_mpnn_n5/outputs,RUN_ID=Brd4ET_2026-09-08_rf3 rf3.slurm
```

## Weights and output folder

| What you set | Env var | Default | Meaning |
|---|---|---|---|
| Checkpoint | `CKPT` | `/opt/databases/foundry/rf3_foundry_01_24_latest_remapped.ckpt` | RF3 weights (Latest, training cutoff 01/24). Relative paths are from the repo root. |
| Folder name | `RUN_ID` | `<jobid>_rf3` | Exact directory under `runs/`. |
| Archive parent | `ARCHIVE_DIR` | `<repo>/archive` | Live rsync destination. |
| Archive refresh | `RSYNC_REFRESH_SEC` | `60` | Seconds between archive refreshes. |

## Folding knobs (Hydra)

If you **leave these unset**, RF3’s own defaults apply: **10 recycles**, **5 diffusion samples**, **200 steps**, early-stop pLDDT **0.5**.

The design pipelines **do** set two of them: `NUM_STEPS=50` and `DIFFUSION_BATCH_SIZE=1`.

| What you set | Env var | RF3 default if unset | Meaning |
|---|---|---|---|
| Recycle iterations | `N_RECYCLES` | 10 | How many times the network re-reads its own prediction. |
| How many diffusion samples | `DIFFUSION_BATCH_SIZE` | 5 | Structures generated per input in one fold call. Pipelines set this to `1`. |
| Diffusion steps | `NUM_STEPS` | 200 | Length of the diffusion trajectory. Pipelines set this to `50` (faster, coarser). |
| Random seed | `SEED` | (RF3 default) | Makes sampling repeatable when set. |
| Give up on bad folds | `EARLY_STOPPING_PLDDT` | 0.5 | Stop a sample if pLDDT falls below this. |
| Do not redo finished inputs | `SKIP_EXISTING` | (unset) | Hydra `skip_existing`. Useful if you restart a long fold into the same `RUN_ID`. |
| Extra Hydra tokens | `RF3_EXTRA` | (empty) | Appended as-is, for example `template_selection="[B]"` to keep templates on chain B. |

## What a job writes

`runs/<id>/{inputs,outputs,slurm.out,manifest.json}` and `archive/<id>/`. Typical RF3 files:

- ranked `*_model.cif`
- `*_confidences.csv`, `*_ranking_scores.csv`
- `*_summary_confidences.json` (pLDDT, pTM, ipTM, … — used later by TM-align)
- `seed-0_sample-n/` per-sample trees

The pipeline’s scorer uses the ranked `*_model.cif` files and ignores paths whose names contain `seed-` or `sample-`.
