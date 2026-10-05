# Boltz2 DD1 refold analysis

Run from the project root:

```bash
python3 03_Filtering/Boltz2_metrics/analyze_boltz2_dd1.py
```

The script indexes the 684 Boltz2 models by their actual chain sequences, so it is robust to the refolding directory-name collisions between Mouse- and Human-derived designs. It writes:

- `out/boltz2_dd1_metrics_long.csv`: one row per design/target pair.
- `out/boltz2_dd1_metrics_paired.csv`: Mouse and Human values side by side for each of the 342 designs.
- `out/scatter_*.svg`: Mouse-EGFR x-axis versus Human-EGFR y-axis plots.

Reported metrics include Boltz2 confidence score, pTM, ipTM, complex and interface pLDDT, complex PDE/IPDE, mean cross-chain PAE (`ipdae_mean_pae`), mean pLDDT, ipSAE, ipTM variants, pDockQ, pDockQ2, LIS, TM-score, and RMSD.

`ipdae_mean_pae` is calculated as the mean PAE across both directions of the binder/EGFR cross-chain blocks. `rmsd_to_rfd3` is a CA RMSD after fitting the complete predicted complex to the original RFD3 complex, using common residue numbers. TM-score and `tm_rmsd` are from Foundry's USalign implementation, using the locally compiled macOS build from the checked-in `USalign.cpp` source.
