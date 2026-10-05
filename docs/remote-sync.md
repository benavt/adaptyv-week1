# Local input and result synchronization

`scripts/sync.sh` runs rsync over a home-directory SSH config alias. All site
values come from the trusted, gitignored `.env.local` Bash file. The example
configuration contains variable names with empty values. SSH addresses, users,
keys, jump-host rules and agent settings belong in your home SSH configuration.

This checkout has a private configuration populated from the supplied file
contract. For another checkout:

```bash
cp .env.example .env.local
chmod 600 .env.local
# Fill in the values locally, then check SSH and the readable output roots:
./scripts/sync.sh doctor
```

Use an SSH key loaded in your agent or keychain. Commands use `BatchMode=yes`,
so missing authentication fails promptly. Python 3 and rsync must be installed
locally, and rsync must exist on the remote host. The Mac's bundled openrsync
is supported: the wrapper uses `--stats` and shell-quotes remote paths. With
a local rsync supporting `--protect-args`, it uses that option instead; the
remote rsync must also support it. Paths containing spaces, quotes and archive
collision suffixes are supported. No transfer deletes source or destination files.

| Variable | Meaning |
| --- | --- |
| `REMOTE` | Short alias already configured in home SSH config |
| `REMOTE_INPUT_ROOT` | Dedicated writable staging parent; never a run/archive tree |
| `REMOTE_FOLDING_ROOT` | Parent of ordinary run directories; source is `<root>/<id>/folding/` |
| `REMOTE_FOUNDRY_ROOT` | Parent of Foundry runs; source is `<root>/<id>/` |
| `REMOTE_CAMPAIGN_ROOT` | Parent of campaign roots, including their shards |
| `REMOTE_LOG_ROOT` | Parent of ordinary `folding-<jobid>.out` logs |
| `REMOTE_RFD3_ROOT` | Optional standalone RFD3 run parent; source is `<root>/<id>/outputs/` |
| `REMOTE_BUNDLE` | Optional exact discovered managed bundle directory |
| `REMOTE_CMD` | Optional trusted command used only by `sync.sh run` |
| `LOCAL_INPUT_ROOT` | Local staging parent; defaults to project `inputs/` |
| `LOCAL_RESULT_ROOT` | Local download parent; defaults to project `results/` |
| `LOCAL_LOG_ROOT` | Local log parent; defaults to project `logs/` |

Relative local roots resolve against the project, regardless of the current
working directory. `SYNC_ENV_FILE` can select another trusted Bash config
instead of the project `.env.local`; keep that config private and gitignored.
Custom download/staging directories need their own ignore rules.

## Upload one experiment

For this EGFR project's preferred new layout, set
`LOCAL_INPUT_ROOT=01_Staging/runs` in a trusted private per-run configuration.
Existing root-level staging runs use `LOCAL_INPUT_ROOT=01_Staging`. Select the
configuration with `SYNC_ENV_FILE` on each invocation; the shell wrapper sources
it after environment setup, so an inline environment prefix cannot reliably
override a configured value. Keep a single unchanged run ID in commands.
The generic `inputs/` default and all existing download routes remain unchanged.
See [PROJECT_ORGANIZATION.md](PROJECT_ORGANIZATION.md).

Choose a fresh run ID using letters, numbers, underscores, dots or hyphens,
starting with a letter or number. Stage real files as:

```text
inputs/<id>/
  folding/
    molecules/          molecular requests, FASTAs or structures only
    assets/             MSA/template/CCD dependencies when needed
    references/         optional geometry references
    reference_map.json  optional
    config.yaml         optional
  foundry/
    design.yaml
    target.pdb
    mpnn_structures/    optional standalone MPNN inputs
    rf3_inputs/         optional standalone RF3 inputs
```

Create only what the experiment uses. Keep configs and reference maps out of
`molecules/`: folding discovery can mistake them for molecular requests.
Request dependency paths must resolve within the uploaded layout. For example,
a request in `molecules/` can use `../assets/target.a3m`. Laptop absolute paths
are not remote dependencies. Materialize local symlinks before staging them;
the uploader rejects symlinks and empty staging directories.

```bash
# First preview of a new destination: create the staging directory only.
./scripts/sync.sh prepare design_comparison_001
./scripts/sync.sh push design_comparison_001 --dry-run
./scripts/sync.sh push design_comparison_001
```

The trailing source slash copies the input contents to the exact run staging
directory. A real push creates the remote staging destination and rejects
symlinked upload ancestors. Dry runs never create remote directories; a preview
can fail if its remote parent has not been created yet. Use `prepare` to create
the input destination before the first preview. Uploading only
places files. Job preparation/submission is a separate step. An optional
locally configured command executes only when you explicitly call:

```bash
./scripts/sync.sh run
```

That command can submit a job if `REMOTE_CMD` is configured to do so; `push`,
`pull`, `prepare`, and `doctor` never submit jobs. Preserve immutable remote input
snapshots and use a fresh run ID for changed settings or archived parents.

## Download indexes, selected predictions or complete runs

```bash
# Ordinary folding run; preview, then collect indexes and execution metadata.
./scripts/sync.sh pull folding design_comparison_001 --metadata --dry-run
./scripts/sync.sh pull folding design_comparison_001 --metadata

# Complete rank-1 predictions with structure/confidence/PAE sidecars.
./scripts/sync.sh pull folding design_comparison_001 --selected --dry-run
./scripts/sync.sh pull folding design_comparison_001 --selected

# Optional exact index filters; ranks are specific to each backend.
./scripts/sync.sh pull folding design_comparison_001 --selected --rank 1 --engine boltz2 --species human

# Full ordinary folding and Foundry trees, materializing linked files/logs.
./scripts/sync.sh pull folding design_comparison_001
./scripts/sync.sh pull foundry design_comparison_001 --dry-run
./scripts/sync.sh pull foundry design_comparison_001

# Compressed RFD3 designs plus JSON metadata, without extracted CIF duplicates.
./scripts/sync.sh pull rfd3 design_comparison_001

# An ordinary folding Slurm log, using the actual job ID.
./scripts/sync.sh log 12345
```

Downloads land in `results/folding/<id>/`, `results/foundry/<id>/`, or
`results/foundry/<id>/rfd3_compressed/`; logs land in `logs/`.
Metadata mode includes optional score/campaign tables when present and omits
native structures, PAE, and source snapshots. It traverses the run tree.

Selected mode uses the already downloaded `results.jsonl` and requires
`status=complete` and the requested rank (default 1). Missing ranks do not match.
It preserves every relative subdirectory and fails if a selected artifact lies
outside the configured source root. Artifact lists are temporary private files.
The selector does not compare confidence across engines, download references,
or establish that a whole experiment succeeded. Pull full runs for reproducibility.

For a sharded campaign, always use the outer campaign root:

```bash
./scripts/sync.sh pull campaign campaign_id --metadata
./scripts/sync.sh pull campaign campaign_id --selected --dry-run
./scripts/sync.sh pull campaign campaign_id --selected
./scripts/sync.sh pull campaign campaign_id
```

Campaign downloads land in `results/campaigns/<id>/`; selection reads
`folding/results.jsonl` and preserves `shards/NNN/folding/...` paths. The full
campaign includes frozen inputs, shard indexes, native artifacts and campaign
metadata. Combined tables alone can be partial while work is active.

## Migrated archives and provenance

```bash
./scripts/sync.sh resolve design_comparison_001
# Inspect the resolved run and archive catalog remotely. Set REMOTE_BUNDLE
# locally to the exact containing managed bundle; do not guess stable IDs.
./scripts/sync.sh pull bundle archive_label --dry-run
./scripts/sync.sh pull bundle archive_label
```

Bundle pulls download the whole configured bundle into
`results/bundles/<archive_label>/`, retaining launcher/inputs/logs/provenance.
`resolve` can point at its `output/` child: verify the containing bundle before
setting `REMOTE_BUNDLE`. Do not upload into managed bundles, retired runs, or
live archive mirrors. Inspect selected trees before pulling: downloads use `-L`
and can materialize externally linked dependencies. Avoid syncing whole
convenience shortcut trees.

Indexes and snapshots retain their original absolute remote paths. Resolve
an artifact relative to the configured remote root and join that relative path
to its local destination; never match only basenames. Sharded runs need the
campaign root mapping. Historical archive paths may need additional mappings
from the archive catalog. Remote collection/scoring/resume should run on the
original tree unless deliberately adapted.

The default `inputs/`, `results/`, `logs/`, and local env files are gitignored
because they can contain private paths or provenance. Do not copy the supplied
site-specific contract or downloaded manifests into public tracked documentation.

## Repeated transfers and verification

For a project-specific subset, `pull ... --files-from <private-list>` downloads
only listed paths relative to that run's source root. Lists can use newline or
NUL separators; absolute paths and traversal components are rejected. Keep
selection lists and comparison receipts under the gitignored results tree.
A subset or delta download is not a complete standalone run: retain a mapping
to the original local files used for comparison.

Rsync is incremental and retains partial transfers. A nonzero SSH/rsync exit
status propagates to the caller. During live jobs, files can change mid-copy;
repeat after completion. Use `--checksum` on a final pull for a deliberate
content comparison, and verify important recorded artifact checksums locally.
It does not itself validate a run's scientific completeness or recorded hashes.

Local verification, including actual rsync through a local SSH stand-in:

```bash
bash -n scripts/sync.sh
python3 -m unittest discover -s scripts/tests -v
```
