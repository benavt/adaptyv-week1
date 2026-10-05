# RFD3 filtering

These scripts perform an initial triage of the non-protonated DD1 RFD3 outputs. They use only the metrics written by RFD3; they do not score EGFR–DD1 interface quality.

From the project root:

```bash
python3 03_Filtering/RFD3_filtering/filter_rfd3.py
```

By default this reads:

- `02_Design/foundry/runs/DD1_Mouse_EGFR_hotspots_09_29/rfd3/outputs`
- `02_Design/foundry/runs/DD1_Human_EGFR_hotspots_09_29/rfd3/outputs`

It writes `summary.csv`, `passing_designs.csv`, and `filter_counts.csv` under `03_Filtering/RFD3_filtering/out/`.

The default hard filters are conservative structural sanity checks: no chain breaks, no backbone or side-chain clashes, and maximum Cα deviation ≤ 2 Å. Thresholds can be overridden on the command line. Passing a design here means “not obviously malformed”; it is not evidence of a favorable interface.

To count designed-chain ASP residues within 5 Å heavy-atom distance of EGFR HIS37 and/or HIS100:

```bash
python3 03_Filtering/RFD3_filtering/DD1_ASP_EGFR_D3_HIS_proximity_filter.py
```

The script writes a per-structure CSV. Use `--cutoff` to change the distance threshold.

To rescreen the 684 Boltz2 refolded DD1/EGFR complexes (342 designs against each species), use:

```bash
python3 03_Filtering/RFD3_filtering/DD1_ASP_EGFR_D3_HIS_proximity_filter.py \
  --refolded-root 03_Filtering/Refolding/Putative_DD1_EGFR_D3 \
  --out 03_Filtering/RFD3_filtering/out/dd1_asp_his_proximity_684_refolded.csv
```

This uses the same 5 Å heavy-atom criterion, with DD1 as chain A and EGFR D3 as chain B in the refolded complexes.

Refolded target species is assigned by exact receptor-sequence matching to
`01_Staging/EGFR/Human_EGFR_310_501_afs3.json` and its Mouse counterpart. Source
design names may contain both species, so they are not used to assign target
identity. An unrecognized or ambiguous receptor sequence stops the calculation
instead of receiving a guessed label. The existing lowercase `species` CSV
values are retained for compatibility.
