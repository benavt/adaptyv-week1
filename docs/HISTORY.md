# Historical EGFR workspace

Use `python3 scripts/workspace.py history [CAMPAIGN]` for the 11 existing
campaigns. Historical stage folders and helpers currently remain in their
original top-level locations; physical relocation into reference/ is pending
a reviewed consumer audit. Do not rename scientific IDs or infer complete
lineage from historical manifests.

The original repository is retained as a read-only GitHub archive. Local
archive branches preserve unpublished Git history; inspect them with git show
and avoid checking them out over the scientific control folder.
Private source/index backups and rollback records are in local/migrations/.
Legacy source omitted from the lightweight primary Git tree remains on disk.
New work belongs in campaigns/ and uses the same skeleton as Week2.
