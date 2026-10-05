# Historical EGFR workspace

The eleven campaigns now live under campaigns/<historical-id>/, each with
01_staging, 02_design, 03_filtering and 04_analysis. Existing scientific IDs
and manifest bytes are preserved. The original manifest for each campaign is
04_analysis/history/manifest.json; its recorded paths are historical evidence.

Use `python3 scripts/workspace.py history [CAMPAIGN]` for the scientific overview.
Use `python3 scripts/artifacts.py list [--campaign CAMPAIGN] [--stage STAGE]`
for compact counts, or add --files for the full per-file catalog.
Use `python3 scripts/artifacts.py locate ORIGINAL_PATH --verify` for the
current canonical location and a checksum of the file present there.
The --verify flag reports current bytes; compare with artifacts.jsonl for the
migration checksum. It does not certify scientific lineage.

Shared tools, targets, registry evidence and legacy configuration remain in
reference/. Two broad folders, RFD3_filtering and PyMol_top20, remain there
because they mix reusable helpers and cross-campaign material. Do not infer
exclusive campaign ownership or rerun archived path-dependent generators.

Scientific files were renamed within one filesystem and checked against inode,
size, timestamp and mode inventories. Each moved file has a SHA-256 catalog.
Navigation symlinks were backed up and replaced by canonical locator records;
actual targets were preserved. Generated runtime build data remains in reference.
Private migration journals, inventories and exact shortcut backups live in
local/migrations/2026-10-05_primary-repository/. Reverse journaled moves only
after checking consumers; restore aliases from their backup only if needed.

Historical campaigns are not retrospectively frozen executions. New attempts
belong under campaigns/<id>/runs/ and freeze inputs/settings with parent links.
Git tracks portable source and reviewed campaign configs; artifact catalogs,
scientific data and private mappings require their own backup. The original
repository is a read-only GitHub archive, with unpublished original refs retained
locally. Never check old archive branches out over the scientific control folder.
