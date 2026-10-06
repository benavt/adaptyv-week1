# DD1_310_phswitch_0217

Placeholder test-panel position **10 of 20**; not an affinity rank or a confirmed submission.

[Back to the panel](../README.md) · [Exact sequence](sequence.fasta) · [Metrics CSV](metrics.csv) · [Source evidence](provenance.json)

| Species | Predicted complex | ipSAE | EGFR H37 pKa | EGFR H100 pKa | DeltaForge Kd (nM) | Files |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| Human | <img src="Human.png" width="280" alt="Human EGFR complex" /> | 0.771 | 5.488 | 5.066 | 20.753 | [PSE](Human.pse) · [PNG](Human.png) · [CIF](Human.cif) |
| Mouse | <img src="Mouse.png" width="280" alt="Mouse EGFR complex" /> | 0.803 | 5.484 | 6.178 | 45.446 | [PSE](Mouse.pse) · [PNG](Mouse.png) · [CIF](Mouse.cif) |

Cyan: binder chain A. Grey: EGFR D3 chain B. Green sticks: EGFR H37/H100 (historical full positions 346/409). Magenta marks the DD2 extension in the D46-130 fusion. Native residue numbering and all native atoms are retained. PSE coordinates are rigidly aligned on EGFR for display; the paired native CIFs are unchanged and were used for scoring.

These are computational predictions. The PNG/PSE view is a presentation of the scored native model, not a pH-conditioned structure.

Affinity workflow: **DeltaForge scoring of supplied structure**. Kd is not an experimental measurement or a prediction at a specified pH.

- Human pKa context: bound protein complex; target EGFR D3 H37/H100. H37 status: `completed`; H100 status: `completed`.
- Mouse pKa context: bound protein complex; target EGFR D3 H37/H100. H37 status: `completed`; H100 status: `completed`.
