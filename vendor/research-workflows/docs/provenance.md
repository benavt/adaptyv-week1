# Trace a final sequence

Keep three identities distinct: a campaign is a scientific question, a run is an
immutable input/settings snapshot, and a candidate is an exact sequence from one
generation origin. Two campaigns may produce the same sequence; keep both origins.
A final selection refuses duplicate sequences and preserves its ranked order.

Candidate and observation records live in the project's `registry/`. They are
private working records unless the project explicitly reviews and allows them in
Git. Retain and back them up with their referenced scientific artifacts. Git code
alone does not restore ignored provenance, runs or results.

## Register the native result

After native design validation, collect the small original FASTA into its run's
`02_design/` stage. Fill the candidate template in editable `drafts/` using the
exact FASTA identity and file SHA-256. The source path is relative to the run and
must name a pipeline stage; physical machine paths never belong in the record.

```sh
python3 scripts/workspace.py candidate register --record drafts/candidate.json
python3 scripts/workspace.py candidate trace --id CANDIDATE
```

Registration verifies the run snapshot, source hash and exact FASTA sequence. It
writes an immutable `registry/candidates/CANDIDATE.json` containing sequence hash,
origin run hash, generation metadata and parent record hashes. Missing generation
fields remain explicit in `generation_unknowns`; supplied metadata is recorded,
not independently certified. Use original native execution receipts, tool versions,
checkpoint identities and settings as generation evidence when available.

A mutation, fusion or other sequence-changing step creates a new candidate. Name
its existing candidate parents in `parents`, and preserve the new native FASTA and
generation method. Changed inputs/settings also create a new child run with
`new ... --parent-run CAMPAIGN/PRIOR_RUN`; run lineage and candidate lineage serve
different purposes. Never overwrite a candidate to attach a changed sequence.

## Record measured or predicted evidence

Collect the exact small score artifact into its run's stage. Fill the observation
template with candidate identity, method/version, scientific scope, file hash,
source locator and native value. Register a new observation for every result.

```sh
python3 scripts/workspace.py candidate observe --record drafts/observation.json
python3 scripts/workspace.py candidate trace --id CANDIDATE
```

JSON locators use JSON Pointer syntax. CSV locators use a column name and physical
row number, starting at 2 after the header; CSV values remain strings. The CLI
checks that the recorded value equals the located source value.

When the score artifact exposes a sequence or sequence hash, add
`candidate_binding` with `kind: "sequence"` or `"sequence_sha256"` and a source
`locator`. The locator reads the same artifact as the score; the expected sequence
or hash comes from the registered candidate. This verifies that the source score
is attached to that exact candidate. A `"native_fasta_id"` binding can check the
native record name, but is weaker because the score source does not expose the
sequence/hash. Without a binding, the association remains researcher-declared;
`score_candidate_binding_unknowns` records the missing identity evidence. Do not
treat a correctly copied score value as proof of its candidate association.

Record species,
construct, pH, assembly and other relevant conditions in `scope`; prediction scores
must remain distinct from experimental measurements. Method version identifies
the actual scoring implementation, not merely a tool's display name.

Observations are immutable files under
`registry/observations/CANDIDATE/OBSERVATION.json`. Trace recursively verifies
source artifacts, parent records and origin run metadata. Keep the referenced
small FASTA and score artifacts available locally for trace verification, even
when bulk predictions remain remote. Transport receipts separately identify
verified remote copies.

## Freeze the ranked lab selection

Prepare a text file with one candidate ID per line in intended rank order. Prepare
ranking JSON containing `method`, `version` and `rationale`; include the actual
criteria, thresholds and manual decisions in that rationale. Optional `evidence`
entries link `candidate_id` and `observation_id` for the observations used to rank.
Missing evidence links remain explicit rather than implying a verified ranking.

```sh
python3 scripts/workspace.py selection --ids-file drafts/ranked-ids.txt --ranking-record drafts/ranking.json
```

This creates `submissions/YYYY-MM-DD_selection_NNN/` with `ranked.csv`, immutable
`selection.json` and a review note. The selection embeds full candidate/parent
traces, observations, ranking rationale and project competition metadata. Later
observations do not change an existing package. It rejects repeated candidate IDs
and repeated sequence hashes.

Review eligibility, track limits, scientific interpretation and final ordering
before submission. The package has not been submitted; its generation metadata
and computational scores do not certify scientific or competition eligibility.
Make another dated package when the chosen set or ranking changes. Preserve an
actual external submission receipt separately when submission occurs.
