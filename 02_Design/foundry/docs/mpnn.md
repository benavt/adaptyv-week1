# MPNN (sequence design)

[`mpnn.slurm`](../mpnn.slurm) runs Foundry `mpnn` on each RFD3 (or other) structure and writes new sequences plus CIFs. Foundry MPNN takes **one structure per call**; the script loops a file or a directory.

Default walltime is **1 hour**. Override with `sbatch --time` when you have many structures or many sequences per structure.

On Brd4ET designs, **chain A is the binder** (usually what you redesign) and **chain B is the target** (left fixed).

## What you must pass

| What you set | Env var | Default | Meaning |
|---|---|---|---|
| Structures to sequence | `INPUTS` | (required) | A `.cif`, `.cif.gz`, or `.pdb`, **or** a directory of those. Gzipped CIFs are unpacked into the run’s `inputs/` snapshot. |

```bash
sbatch --export=ALL,INPUTS=runs/Brd4ET_2026-09-08/outputs,RUN_ID=Brd4ET_2026-09-08_mpnn,DESIGNED_CHAINS=A,N_SEQS=5 mpnn.slurm
```

## Sequence knobs

| What you set | Env var | Default | Meaning |
|---|---|---|---|
| How many sequences per structure | `N_SEQS` | `1` | Passed to MPNN as `--number_of_batches`. Each structure gets this many redesigns. |
| Same as `N_SEQS` | `N_BATCHES` | (alias) | Used only if `N_SEQS` is unset. |
| Which chains to redesign | `DESIGNED_CHAINS` | `A` | Polymer chain ids. The target chain is not in this list, so it stays the input sequence. |
| Skip the interface | `EXCLUDE_INTERFACE` | unset / off | Set to `True` / `1` / `yes` / `on` to redesign **only non-interface** residues on `DESIGNED_CHAINS`. Anything else leaves the flag off. |
| How close counts as interface | `INTERFACE_CUTOFF` | `5.0` | Ångströms. A residue is interface (and stays fixed) if **any of its atoms** is this close to **any atom on another polymer chain**. Waters are ignored. |

`EXCLUDE_INTERFACE=True` is what the full design pipeline uses. Foundry MPNN cannot combine `--designed_chains` with `--designed_residues`, so the job calls [`interface_residues.py`](../interface_residues.py) per structure, writes `<structure>_interface.json` under `inputs/`, and passes the non-interface list as `--designed_residues`. If every residue on the designed chain is at the interface, that structure fails.

Reusing the same `RUN_ID` **appends** more sequences to existing `*.fa` files. Use a new `RUN_ID` for a clean set.

## Weights

| What you set | Env var | Default | Meaning |
|---|---|---|---|
| Checkpoint | `CKPT` | `/opt/databases/foundry/proteinmpnn_v_48_020.pt` | MPNN weights. Relative paths are from the repo root. |
| Model family | `MODEL_TYPE` | `protein_mpnn` | `protein_mpnn` or `ligand_mpnn`. |
| Old-format weights | `IS_LEGACY_WEIGHTS` | `True` | Must stay `True` for vanilla `proteinmpnn_v_48_*.pt` and for the repo-local SolubleMPNN file. |

### SolubleMPNN (what the design pipeline uses)

Repo-local file `solublempnn_v_48_020.pt`: vanilla ProteinMPNN-format weights (48 nearest neighbors, σ = 0.20 Å noise in training). Use:

| Setting | Value |
|---|---|
| `CKPT` | `solublempnn_v_48_020.pt` |
| `MODEL_TYPE` | `protein_mpnn` |
| `IS_LEGACY_WEIGHTS` | `True` |

```bash
sbatch --export=ALL,INPUTS=runs/Brd4ET_2026-09-08/outputs,RUN_ID=Brd4ET_2026-09-08_mpnn_soluble,DESIGNED_CHAINS=A,N_SEQS=5,CKPT=solublempnn_v_48_020.pt,MODEL_TYPE=protein_mpnn,IS_LEGACY_WEIGHTS=True mpnn.slurm
```

## Job bookkeeping

| What you set | Env var | Default | Meaning |
|---|---|---|---|
| Folder name | `RUN_ID` | `<jobid>_mpnn` | Exact directory under `runs/`. |
| Archive parent | `ARCHIVE_DIR` | `<repo>/archive` | Live rsync destination. |
| Archive refresh | `RSYNC_REFRESH_SEC` | `60` | Seconds between archive refreshes. |

## What a job writes

`runs/<id>/{inputs,outputs,slurm.out,manifest.json}` plus `archive/<id>/`. Outputs typically include `*.fa` and structure CIFs (`--write_fasta True --write_structures True`). Interface JSON sidecars land next to the snapped inputs when `EXCLUDE_INTERFACE` is on.

## Skip-interface example

```bash
sbatch --export=ALL,INPUTS=runs/Brd4ET_2026-09-08/outputs,RUN_ID=Brd4ET_2026-09-08_mpnn_noiface,DESIGNED_CHAINS=A,N_SEQS=5,EXCLUDE_INTERFACE=True mpnn.slurm
```
