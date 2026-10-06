# DD1_062

Placeholder test-panel position **07 of 20**; not an affinity rank or a confirmed submission.

[Back to the panel](../README.md) · [Exact sequence](sequence.fasta) · [Metrics CSV](metrics.csv) · [Source evidence](provenance.json)

| Species | Predicted complex | ipSAE | EGFR H37 pKa | EGFR H100 pKa | DeltaForge Kd (nM) | Files |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| Human | <img src="Human.png" width="280" alt="Human EGFR complex" /> | 0.687 | 6.993 | 6.936 | 212.242 | [PSE](Human.pse) · [PNG](Human.png) · [CIF](Human.cif) |
| Mouse | <img src="Mouse.png" width="280" alt="Mouse EGFR complex" /> | 0.825 | 6.793 | 6.529 | 207.449 | [PSE](Mouse.pse) · [PNG](Mouse.png) · [CIF](Mouse.cif) |

Cyan: binder chain A. Grey: EGFR D3 chain B. Green sticks: EGFR H37/H100 (historical full positions 346/409). Magenta marks the DD2 extension in the D46-130 fusion. Native residue numbering and all native atoms are retained. PSE coordinates are rigidly aligned on EGFR for display; the paired native CIFs are unchanged and were used for scoring.

These are computational predictions. The PNG/PSE view is a presentation of the scored native model, not a pH-conditioned structure.

Affinity workflow: **DeltaForge scoring of supplied structure**. Kd is not an experimental measurement or a prediction at a specified pH.

- Human pKa context: EGFR-only conformation extracted from complex; target EGFR D3 H37/H100. H37 status: `reported`; H100 status: `reported`.
- Mouse pKa context: EGFR-only conformation extracted from complex; target EGFR D3 H37/H100. H37 status: `reported`; H100 status: `reported`.
