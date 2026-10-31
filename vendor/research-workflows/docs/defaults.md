# Pipeline defaults

Explicit user settings win. Otherwise use the selected campaign recipe and this
operational policy, then the audited adapter defaults. Freeze all effective values
with the run before submission. These instructions guide command preparation;
the workspace front door does not execute engines or inject CLI arguments.

| Step | Default for new generic work | Required scientific choices |
| --- | --- | --- |
| Staging | Preserve originals, validate coordinates/chain and residue mapping, create an immutable input snapshot | Target structures, constructs, hotspots and constraints |
| Design with Foundry | RFD3; **1,000 requested backbones**, passed explicitly as the helper count | Validated native inputs and campaign-specific design settings |
| Boltz2 prediction using the audited common adapter | 200 sampling steps, 1 diffusion sample, 3 recycles; common runner seed 0; preserve supplied MSA/templates | Prediction inputs and model/checkpoint identity |
| Filtering | Score the completed records locally and retain score evidence and candidate lineage; preserve all candidates until a selection rule is defined | Score method, thresholds and eligibility criteria |
| Analysis | Local analysis, dated outputs, shared figure/presentation standards, traceable ranked selection | Comparison objective and final panel size |

For Foundry, pass `--count 1000` through the supported common launcher; its legacy
helper default is 8 and must not silently control a generic run. Count is a
requested RFD3 backbone budget, not successful designs or final sequences. MPNN
can generate multiple sequences per backbone; multiple native specifications may
multiply output counts. Record count scope and actual completed counts separately.
Do not invent an MPNN sequence count for a generic request.

The Boltz2 values describe the audited common adapter, not every native/versioned
route. The native Boltz2 seed defaults to unset; pass and record the chosen seed
explicitly when using that route. Alternate engines retain their own validated
presets. Do not copy smoke-test settings or unrelated campaign thresholds.

Submit the reviewed job through the project's supported adapter; it invokes
`sbatch` and retains the accepted job ID. Slurm owns queue order and resource
scheduling. Reuse the selected runner's verified resource profile; explicit
user requests override it. Add stage measurement with the existing recording
helper when validating a native adapter. No separate dispatcher or idle-triggered
launcher is part of this workflow. Resolve host and SSH/Tailscale access details
through Fleet's private inventory and bindings.
