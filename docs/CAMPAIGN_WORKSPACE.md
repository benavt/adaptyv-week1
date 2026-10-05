# Week1 campaign workspace

Start with `python3 scripts/workspace.py history` to browse the 11 existing
campaigns. Add a campaign ID to inspect its original stage locations, designated
outputs and recorded runs. These manifests are historical navigation evidence;
their statuses are recorded local observations. No historical campaign or
scientific file has been moved or converted into an immutable snapshot.

`python3 scripts/workspace.py show` describes the new four-stage workspace.
`show --stage 02_design` loads only the selected shared stage and pinned guide.
New work can use the shared campaign/run commands from the same entry point.
Use `campaign`, then prepare explicit inputs under drafts/ and create a new run
with `new`. Choose each campaign's scientific recipe before freezing it: this
project deliberately has no default design engine, construct, scoring method,
seed, checkpoint, or acceptance threshold. Existing native runners retain their
original input/result path contracts. New native execution adapters must be
validated before scientific compute; transfer commands do not submit jobs.

New campaigns live under campaigns/<campaign>/; date-first ISO run snapshots
live under campaigns/<campaign>/runs/. Changed inputs/settings always create a
successor run. Preserve stable historical design/model/sequence IDs and the
existing registry. Shared candidate records concern new frozen runs; they do not
replace the existing registry or fill missing historical execution evidence.
Inspect registry/README.md and original campaign evidence for historical identity.

The project pins shared code and file hashes in workflow.lock.json. Existing guide
pins from the parent Git index have private read-only releases under
local/guide-releases/. Live guide submodule checkouts remain untouched, including
any existing divergence from their parent pin. The new workspace reads its pinned
release; legacy launchers still use their original guide contracts. Guide updates
and environment changes must be deliberate and recorded in successor runs.

Actual Fleet/shared/interpreter paths are configured only in ignored
local/fleet.json. The front door uses its explicit isolated control Python.
Local CPU analysis and rendering require method-specific dependency checks;
remote GPU environments are selected by campaign and resolved through Fleet.
No GPU engine is installed on the Mac by this adoption.

Use `standard figures`, `standard molecular-views` or `standard presentations`
for new artifacts. The shared docs/provenance.md describes candidate lineage,
source-backed observations and frozen selection packages for new runs. Existing
historical selections keep their original source evidence and scientific caveats.

The original AGENTS.md scientific instructions, sync configurations, registry,
application pointers, four stage directories and scientific artifacts remain in
place. The repository now lives in the standard local Code project home. Historical stage paths inside it remain unchanged; further artifact moves require a separate consumer and backup audit.

The primary GitHub source is adaptyv-week1 under the personal profile. The old
organization repository is retained as a read-only historical archive. Legacy
submodule checkouts remain local; new work resolves the deliberately pinned
guide releases in workflow.lock.json, not their dirty working trees.
