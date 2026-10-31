# Execution contract

This document defines portable expectations for remote execution. Research
projects own their scientific workflow and thin project-specific adapters;
`research-workflows` owns the shared lifecycle and evidence contract. Fleet is
limited to machine inventory and SSH/Tailscale bindings and access. Adapters
must not delegate research policy or execution semantics to Fleet.

Adoption is in progress: some existing project adapters may still use Fleet for
launch or monitoring. Treat that as legacy behavior to migrate deliberately,
not as the target ownership model. Record each project's adapter and explicit
overrides beside its pinned `external/research-workflows` revision.

## Prepare and monitor

- Record resource use and errors across the complete job: setup, CPU preparation,
  workload, arrays and finalization. A stage-only log is incomplete. Use the
  project's supported adapter and pinned checker to validate each running task's
  fresh telemetry and its finished log. Missing samples remain gaps, not zeroes.
- Supported GPU folding requires verified CPU MSA preparation and native CPU
  feature preparation before GPU allocation. GPU inference consumes the prepared
  cache; it must not silently repeat CPU preparation. Unsupported engines or
  adapters stop before submission until their preparation path is validated.
- Preserve the selected chemistry, native runtime and package pins through every
  stage. Freeze inputs and settings in a new immutable run when they change; keep
  submission intent, scheduler response, resource logs and validation receipts.
- A scheduler's finished state does not prove valid telemetry. Valid telemetry
  does not prove scientific success. Record scientific acceptance separately from
  both operational states. High-level status summaries belong in run records or
  project views, not in resource logs or Notion as an authority for job evidence.
- Do not bypass the project adapter with raw `sbatch`, and do not blindly retry
  an ambiguous submission. Inspect its receipt and scheduler state first.

Resolve host and connection details through Fleet's inventory and access
bindings. Keep private host names, paths, credentials and connection settings
out of shared policy. Site-specific launch and monitoring belong in the
project's adapter, under this shared contract. For any credit-consuming
LigandAI operation, first quote the call count,
expected credits, project/run, operation and spend ceiling; wait for explicit
approval of that batch before starting it. Paid language-model inference also
requires the user's explicit request for that model call; never use inference to
probe access. Follow Fleet's credit policy for the full approval procedure.
