# New design/mutant checklist

## 1. Define and stage

- [ ] Assign a stable ID and mutation notation.
- [ ] Confirm the complete binder sequence.
- [ ] Confirm Human and Mouse EGFR source sequences.
- [ ] Create the paired FASTA inputs with binder first.
- [ ] Copy source structures without overwriting them.
- [ ] Check atom names, chain IDs, residue numbering, ligands, and glycans.
- [ ] Start a run record from `RUN_RECORD_TEMPLATE.md`.

## 2. Design or import structures

- [ ] Validate Foundry YAML/PDB compatibility and hotspot atom names.
- [ ] Use a unique run ID.
- [ ] Record checkpoint, environment, job ID, and manifest.
- [ ] If starting mid-pipeline, record the imported output directory and stage.

## 3. Filter

- [ ] Run structural sanity filtering.
- [ ] Run interface/proximity filtering when relevant.
- [ ] Save threshold values and retained/rejected counts.
- [ ] Treat filters as triage, not proof of binding.

## 4. Refold and score

- [ ] For structure inference, apply `STRUCTURE_INFERENCE.md`: record the selected
      engine, resolved intake, validation level, execution and collected outcomes.

- [ ] Run Boltz2 for both species.
- [ ] Pair results by stable design ID and actual sequence.
- [ ] Verify chain A=binder and chain B=EGFR in scoring inputs.
- [ ] Run DeltaForge with receptor B and peptide A for canonical complexes.
- [ ] Prepare/run JustHISpKa inputs and label any proxy or placeholder.
- [ ] Run USalign with the reference orientation recorded.
- [ ] Preserve raw responses and normalized outputs.

## 5. Rank

- [ ] Apply the versioned ranking/filter specification.
- [ ] Record all cutoffs, weights, and candidate counts.
- [ ] Do not overwrite prior ranking outputs; create a new version.

## 6. Render and report

- [ ] Render Human and Mouse structures using the fixed species views.
- [ ] Check colors, hotspot selections, residue labels, and aspect ratios.
- [ ] Build the score table from CSV/JSON outputs, not typed values.
- [ ] Include the complete binder sequence, color key, and method note.
- [ ] Run the QA checklist in `REPORTING_SPEC.md`.

## 7. Archive

- [ ] Complete the run record.
- [ ] Confirm all commands, versions, inputs, and outputs are traceable.
- [ ] Keep raw and derived files together.
- [ ] Record unresolved caveats before sharing the result.
