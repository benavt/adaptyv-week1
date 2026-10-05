# ProtonPottsMPNN local runtime

The local runtime is kept outside version control under
`03_Filtering/ProtonPottsMPNN_runtime/`.

```text
Python:       03_Filtering/ProtonPottsMPNN/.venv/bin/python
HBPLUS_PATH:  03_Filtering/ProtonPottsMPNN_runtime/build/hbplus/hbplus
HBPLUS:       v3.06 source distribution, compiled locally on macOS arm64
XGBoost:      2.1.4 (required for the shipped serialized EV6 acid models)
```

The environment was installed from the submodule's `foundry` package and
`requirements-extra.txt` using Python 3.12. XGBoost 3.4.1 satisfied the loose
upstream requirement but was incompatible with the serialized acid model
objects, so it was replaced with 2.1.4.

Use this setup before labeling:

```bash
export HBPLUS_PATH="$PWD/03_Filtering/ProtonPottsMPNN_runtime/build/hbplus/hbplus"
export CUDA_VISIBLE_DEVICES=""
unset DEBUG EV6_HIS_PROB_THR EV6_ACID_PROB_THR
```

HBPLUS was checked on the DD2 `top_07` complex. It completed successfully and
reported 253 hydrogen bonds.

## Separate pKa calculation

ProtonPottsMPNN's EV6 labeller predicts a structural token and an ensemble score
(`p_protonated`, `sd`). It has no pH input and does not solve an equilibrium
microstate or return pKa values. To estimate states at pH 5.0 and 7.4, run a
separate pKa calculation for every titratable binder site, preferably in both
the bound complex and the fixed-coordinate binder-only structure. Convert a
site pKa to an approximate protonated fraction with:

```text
fraction_protonated = 1 / (1 + 10**(pH - pKa))
```

The existing JustHISpKa setup is suitable for histidines, but it requires a
properly parameterized Tripos MOL2 plus Amber `prmtop` and is not a general
Asp/Glu calculation. Use PROPKA 3.5, which is installed in the runtime, or
another validated structure-based pKa engine for Asp/Glu. Record the engine,
input structure, chain/residue mapping, and any coupled-group or forced
microstate settings. If the switch analysis depends on coupled titration, use a
multi-site method or enumerate coupled microstates rather than treating sites
as independent.
