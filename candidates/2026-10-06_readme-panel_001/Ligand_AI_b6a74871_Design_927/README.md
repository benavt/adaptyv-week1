# Ligand_AI_b6a74871_Design_927

Placeholder test-panel position **02 of 20**; not an affinity rank or a confirmed submission.

[Back to the panel](../README.md) · [Exact sequence](sequence.fasta) · [Metrics CSV](metrics.csv) · [Source evidence](provenance.json)

| Species | Predicted complex | ipSAE | EGFR H37 pKa | EGFR H100 pKa | DeltaForge Kd (nM) | Files |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| Human | <img src="Human.png" width="280" alt="Human EGFR complex" /> | 0.934 | 3.490 | 4.238 | 55.024 | [PSE](Human.pse) · [PNG](Human.png) · [CIF](Human.cif) |
| Mouse | <img src="Mouse.png" width="280" alt="Mouse EGFR complex" /> | 0.876 | 3.595 | 4.614 | 232.734 | [PSE](Mouse.pse) · [PNG](Mouse.png) · [CIF](Mouse.cif) |

Cyan: binder chain A. Grey: EGFR D3 chain B. Green sticks: EGFR H37/H100 (historical full positions 346/409). Magenta marks the DD2 extension in the D46-130 fusion. Native residue numbering and all native atoms are retained. PSE coordinates are rigidly aligned on EGFR for display; the paired native CIFs are unchanged and were used for scoring.

These are computational predictions. The PNG/PSE view is a presentation of the scored native model, not a pH-conditioned structure.

Affinity workflow: **DeltaForge scoring of supplied structure**. Kd is not an experimental measurement or a prediction at a specified pH.

- Human pKa context: packaged target-chain EGFR D3 histidine predictions; see evidence/scoring_methods.json. H37 status: `reported`; H100 status: `reported`.
- Mouse pKa context: packaged target-chain EGFR D3 histidine predictions; see evidence/scoring_methods.json. H37 status: `reported`; H100 status: `reported`.
