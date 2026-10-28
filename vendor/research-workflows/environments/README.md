# Environment roles

Select an environment by task role: `control`, `analysis`, `molecular_views`,
`foundry` or `boltz2`. The role names in roles.json are portable. Fleet's ignored
`local/environments.json` binds each project's roles to actual machine prefixes,
explicit interpreters, ownership, observed versions and inventory hashes. A role
binding records where to execute; its validation status says what was checked.

Use the bound interpreter directly. A project's front door can select its control
interpreter from ignored local configuration and re-execute itself in isolated
Python mode. Do not infer the correct environment from the current shell's `python`
or `conda`: multiple managers can coexist. Local CPU and remote GPU roles may need
different Python, package and CUDA versions. Do not force them into one environment.

Installed package metadata is an inventory, not a complete rebuild lock. It does
not establish working CUDA, checkpoint identity, source checkout compatibility,
renderer availability or the dependencies of a particular analysis method.
Check those dependencies before using the scientific method. Reuse an existing
validated role where possible; avoid installing duplicate engines and models.


For the current project, Foundry and Boltz2 execute on the remote GPU resource.
Do not install or bind either engine on the Mac. The Mac roles are control,
selected CPU analysis and molecular rendering. Existing local testing environments
can be retired only after their consumers and unique contents are audited.

User-owned remote environment prefixes use the existing registered archive
storage, rather than the quota-limited home directory. Institutional/admin-managed
prefixes are an explicit exception: register and preserve their actual locations.
An existing alias or symlink can point to such a shared environment. Do not move
shared software merely to make the visible paths uniform. Audit alias consumers
before replacing an existing link. This phase creates no new environment symlinks.

## Deliberate rebuild

1. Audit launchers, active jobs and scripts consuming the current prefix, including
   existing symlinks. Preserve shared/admin-managed environments without mutation.
2. Record the intended direct dependencies in a portable Conda definition, with
   an explicit Python version and engine/source identities. Separate CPU and GPU
   definitions where requirements differ; keep actual paths and credentials out.
3. Resolve and retain a platform-specific package lock after dependencies and
   channels are reviewed. Record pip packages and source commits separately where
   needed. A generic YAML or copied installed-package list is not a verified lock.
4. Create a new prefix in the machine's registered environment storage; validate
   preparation and scientific or rendering behavior before deliberately rebinding
   the role. Keep the old prefix and rollback receipt until consumers are migrated.
5. Record the resolved environment identity and package-inventory hash in execution
   receipts. A changed environment is a changed scientific setting: freeze a new
   run snapshot and retain the previous run's identity and outputs.

This adoption phase registers existing prefixes. It does not install, delete,
rename or move environments, change shell activation, or remove existing symlinks.
Consolidating managers and rebuilding environments are later audited phases.
