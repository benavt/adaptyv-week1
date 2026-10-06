# rank_05_D46_130_model_1_chain_A

Placeholder test-panel position **04 of 20**; not an affinity rank or a confirmed submission.

[Back to the panel](../README.md) · [Exact sequence](sequence.fasta) · [Metrics CSV](metrics.csv) · [Source evidence](provenance.json)

| Species | Predicted complex | ipSAE | EGFR H37 pKa | EGFR H100 pKa | DeltaForge Kd (nM) | Files |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| Human | <img src="Human.png" width="280" alt="Human EGFR complex" /> | 0.468 | 4.400 | 7.029 | 34.543 | [PSE](Human.pse) · [PNG](Human.png) · [CIF](Human.cif) |
| Mouse | <img src="Mouse.png" width="280" alt="Mouse EGFR complex" /> | 0.698 | 4.333 | — | 130.351 | [PSE](Mouse.pse) · [PNG](Mouse.png) · [CIF](Mouse.cif) |

Cyan: binder chain A. Grey: EGFR D3 chain B. Green sticks: EGFR H37/H100 (historical full positions 346/409). Magenta marks the DD2 extension in the D46-130 fusion. Native residue numbering and all native atoms are retained. PSE coordinates are rigidly aligned on EGFR for display; the paired native CIFs are unchanged and were used for scoring.

These are computational predictions. The PNG/PSE view is a presentation of the scored native model, not a pH-conditioned structure.

Affinity workflow: **DeltaForge sequence fold-and-score**. Kd is not an experimental measurement or a prediction at a specified pH.

- Human pKa context: whole bound complex; target EGFR D3 H37/H100. H37 status: `reported`; H100 status: `reported`.
- Mouse pKa context: whole bound complex; target EGFR D3 H37/H100. H37 status: `reported`; H100 status: `unresolved; conditional states available`.
  Conditional H100, acid deprotonated: 8.173. This assumption is kept separate from the displayed unforced value.
  Conditional H100, acid protonated: 2.338. This assumption is kept separate from the displayed unforced value.
