# DD1 top-20 selection plan

The current grid search uses the 342 paired Mouse/Human rows.

1. Apply symmetric hard cutoffs to both species: `ipSAE_mouse > cutoff` and `ipSAE_human > cutoff`; likewise for TM-score.
2. Search ipSAE cutoffs from 0.30 to 0.70 in 0.01 increments and TM-score cutoffs from 0.55 to 0.80 in 0.01 increments.
3. Choose the cutoff pair whose retained count is closest to 20. The current pair is ipSAE > 0.59 and TM-score > 0.63, retaining exactly 20 designs.
4. Rank those retained designs using the mean of five min-max normalized components: mean ipSAE, minimum ipSAE across species, mean TM-score, minimum TM-score across species, and inverse mean ipDAE. Higher composite score is better. The minimum terms penalize Mouse/Human imbalance, while inverse ipDAE rewards lower interface PAE.
5. Report the top 20 sequences with both per-species metrics and the composite score.

This is a selection heuristic, not a calibrated probability of binding. A future iteration can change the grid, use a Pareto frontier, or alter the component weights without overwriting the current outputs.
