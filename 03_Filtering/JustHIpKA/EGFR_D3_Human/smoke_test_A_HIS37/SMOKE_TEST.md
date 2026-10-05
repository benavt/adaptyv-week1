# justHISpKa smoke test: EGFR D3 human A:HIS37

Date: 2026-09-30

## Preparation

Input: `../EGFR_D3_Human_justhispka.pdb`

AmberTools 26 successfully completed:

1. `pdb4amber` validation/conversion;
2. `tleap` with `leaprc.protein.ff14SB`;
3. `ambpdb -mol2 -sybyl` conversion to MOL2.

Generated artifacts:

- `egfr_d3_human_A_HIS37.prmtop`
- `egfr_d3_human_A_HIS37.inpcrd`
- `egfr_d3_human_A_HIS37.raw.mol2`

The raw MOL2 contains histidine substructure ID `37`, corresponding to PDB residue A:HIS37. The package's `push_renum.sh` helper was also attempted, but it rejected the current AmberTools 26 `pdb4amber` renumber-file format (`Incorrect file`). This is a separate helper-script compatibility issue; it did not prevent AmberTools from producing the parameterized files.

## Executable test

Command attempted:

```bash
arch -x86_64 /.../software/justhispka-1.0.2/justhispka-1.0.2/mac_justHisPKa \
  MOL=egfr_d3_human_A_HIS37.raw.mol2 \
  PRM=egfr_d3_human_A_HIS37.prmtop \
  SID=37 \
  LICENSE=/.../JustHIpKA/justHisPKa.lic
```

Initial result: the process terminated with exit status **137** and emitted no pKa output. Investigation showed that macOS had attached a quarantine attribute to the downloaded executable. Rosetta itself was confirmed functional with translated system binaries. After removing the quarantine attribute, the same command succeeded:

```text
PKA=7.067 SID=37 RES=******:HIS
```

The executable reports the correct SID and pKa, but does not preserve the PDB chain/residue label in this raw MOL2 conversion (`RES=******:HIS`). The PDB-to-MOL2 renumbering helper should be repaired or replaced before relying on chain-aware labels in production.

## Required software/runtime to complete the test

One of the following is needed:

- an Intel macOS machine capable of running the supplied `mac_justHisPKa` binary;
- an x86-64 Linux machine/VM/container capable of running the supplied `justHisPKa` binary; or
- a native Apple Silicon/arm64 build of justHISpKa from NRC.

AmberTools 18+ and PyMOL 2+ are not blockers here: AmberTools 26 and PyMOL 3.1.0 are installed, Rosetta 2 is installed and functional, and the AmberTools preparation and justHISpKa prediction both succeeded after clearing the executable quarantine attribute.
