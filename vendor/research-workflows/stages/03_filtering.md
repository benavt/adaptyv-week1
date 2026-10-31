# Filtering

Follow the campaign recipe for refolding, structural scores and selection.
Resolve native chain roles from exact sequences; preserve confidence/PAE sidecars
and Human/Mouse pairing where requested. Record score definitions, units, cutoffs,
method/version and source file hashes. Missing values remain missing.

GPU prediction and local CPU scoring can be separate substeps of this stage.
Keep raw predictions remote when bulky; collect explicit files for local methods.
Keep retained/rejected counts and selection evidence. A confidence score or pKa
proxy does not establish binding or experimental pH selectivity. Outputs belong
in this run's 03_filtering/, preserving separate prediction/scoring attempts.

For generic new work, read [pipeline defaults](../docs/defaults.md); explicit
user settings override them. Record effective settings before submission.
