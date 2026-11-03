# DD1_195_phswitch_0591

Placeholder test-panel position **12 of 20**; not an affinity rank or a confirmed submission.

[Back to the panel](../README.md) · [Exact sequence](sequence.fasta) · [Metrics CSV](metrics.csv) · [Source evidence](provenance.json)

| Species | Predicted complex | ipSAE | EGFR H37 pKa | EGFR H100 pKa | DeltaForge Kd (nM) | Files |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| Human | <img src="Human.png" width="280" alt="Human EGFR complex" /> | 0.862 | 6.187 | 1.456 | 136.744 | [PSE](Human.pse) · [PNG](Human.png) · [CIF](Human.cif) |
| Mouse | <img src="Mouse.png" width="280" alt="Mouse EGFR complex" /> | 0.917 | 6.900 | 4.812 | 113.157 | [PSE](Mouse.pse) · [PNG](Mouse.png) · [CIF](Mouse.cif) |

Cyan: binder chain A. Grey: EGFR D3 chain B. Green sticks: EGFR H37/H100 (historical full positions 346/409). Magenta marks the DD2 extension in the D46-130 fusion. Native residue numbering and all native atoms are retained. PSE coordinates are rigidly aligned on EGFR for display; the paired native CIFs are unchanged and were used for scoring.

These are computational predictions. The PNG/PSE view is a presentation of the scored native model, not a pH-conditioned structure.

Affinity workflow: **DeltaForge scoring of supplied structure**. Kd is not an experimental measurement or a prediction at a specified pH.

- Human pKa context: bound protein complex; target EGFR D3 H37/H100. H37 status: `completed`; H100 status: `completed`.
- Mouse pKa context: bound protein complex; target EGFR D3 H37/H100. H37 status: `completed`; H100 status: `completed`.
