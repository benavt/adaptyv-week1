# Vendored Research Workflows

This directory contains a sanitized snapshot of selected files from the private `benavt/research-workflows` repository.

## Source

- **Repository**: `benavt/research-workflows` (private)
- **Commit**: `e2e2f3ac9d2f6dcd5965c16da957c93d37df4d4c`
- **Date**: October 16, 2026
- **Branch**: main

## Vendoring Rationale

The research-workflows repository is private and not distributed publicly. This snapshot contains only the workflow infrastructure files needed by Week1 EGFR project, with sensitive information removed or neutralized.

## Sanitization

All files have been scanned and cleaned of:
- Private infrastructure hostnames and embedded references
- User-specific paths
- SSH configuration
- Private IP addresses and Tailscale references
- Internal git hosts
- Institution-specific identifiers

Test fixtures that referenced private infrastructure have been updated to use neutral values (e.g., `remote-gpu` instead of actual hostnames).

## Files Not Vendored

Files containing inherently private infrastructure details or not required by Week1 were deliberately excluded. The .git history is not vendored.

## Usage

Week1 scripts reference these vendored files directly. The workflow binding layer automatically resolves to this vendor directory when the repository is distributed publicly.

## Updates

This is a point-in-time snapshot. Updates require manual re-vendoring from the private repository with the same sanitization process.
