#!/usr/bin/env python3
"""Select non-interface residues on designed chains for Foundry MPNN."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import gzip

import biotite.structure as struc
import biotite.structure.io as bsio
from biotite.structure.io.pdbx import CIFFile
import numpy as np

WATER_NAMES = {"HOH", "WAT", "DOD", "H2O"}


def residue_mpnn_id(chain_id: str, res_id: int, ins_code: str) -> str:
    ins = (ins_code or "").strip()
    if ins in {".", "?"}:
        ins = ""
    return f"{chain_id}{int(res_id)}{ins}"


def _read_cif_file(path: Path) -> CIFFile:
    name = path.name.lower()
    if name.endswith(".gz"):
        with gzip.open(path, "rt") as handle:
            return CIFFile.read(handle)
    return CIFFile.read(str(path))


def load_atoms_from_cif(path: Path):
    """Read coords from atom_site.Cartn_* so extra MPNN CIF columns do not NaN biotite."""
    cif = _read_cif_file(path)
    site = cif.block["atom_site"]
    if "pdbx_PDB_model_num" in site:
        models = site["pdbx_PDB_model_num"].as_array(str)
        keep = models == models[0]
    else:
        keep = slice(None)
    chain = site["auth_asym_id"].as_array(str)[keep] if "auth_asym_id" in site else site["label_asym_id"].as_array(str)[keep]
    res_id = site["auth_seq_id"].as_array(int)[keep] if "auth_seq_id" in site else site["label_seq_id"].as_array(int)[keep]
    res_name = site["auth_comp_id"].as_array(str)[keep] if "auth_comp_id" in site else site["label_comp_id"].as_array(str)[keep]
    atom_name = site["auth_atom_id"].as_array(str)[keep] if "auth_atom_id" in site else site["label_atom_id"].as_array(str)[keep]
    if "pdbx_PDB_ins_code" in site:
        ins = site["pdbx_PDB_ins_code"].as_array(str)[keep]
        ins = np.array(["" if x in {".", "?", ""} else x for x in ins], dtype="U4")
    else:
        ins = np.full(len(chain), "", dtype="U4")
    coord = np.column_stack(
        [
            site["Cartn_x"].as_array(float)[keep],
            site["Cartn_y"].as_array(float)[keep],
            site["Cartn_z"].as_array(float)[keep],
        ]
    )
    finite = np.isfinite(coord).all(axis=1)
    if not np.any(finite):
        raise ValueError(f"no finite coordinates in {path}")
    array = struc.AtomArray(int(finite.sum()))
    array.coord = coord[finite]
    array.chain_id = chain[finite]
    array.res_id = res_id[finite]
    array.res_name = res_name[finite]
    array.atom_name = atom_name[finite]
    array.ins_code = ins[finite]
    return array


def load_atoms(path: Path):
    name = path.name.lower()
    if name.endswith(".pdb") or name.endswith(".pdb.gz"):
        array = bsio.load_structure(str(path), model=1)
        if isinstance(array, struc.AtomArrayStack):
            array = array[0]
        if not np.isfinite(array.coord).all():
            raise ValueError(f"non-finite coordinates in {path}")
        return array
    return load_atoms_from_cif(path)


def is_water(array) -> np.ndarray:
    names = np.char.upper(np.asarray(array.res_name, dtype=str))
    return np.isin(names, list(WATER_NAMES))


def select_designed_residues(
    path: Path,
    chains: list[str],
    cutoff: float = 5.0,
) -> dict:
    array = load_atoms(path)
    water = is_water(array)
    polymer = array[~water]
    if polymer.array_length() == 0:
        raise ValueError(f"no non-water atoms in {path}")

    chain_set = set(chains)
    designed_mask = np.isin(polymer.chain_id.astype(str), list(chain_set))
    other_mask = ~designed_mask
    if not np.any(designed_mask):
        raise ValueError(f"no atoms for designed chains {chains} in {path}")
    if not np.any(other_mask):
        raise ValueError(f"no other-chain atoms to define an interface in {path}")

    other_coord = polymer.coord[other_mask]
    designed = polymer[designed_mask]
    starts = struc.get_residue_starts(designed, add_exclusive_stop=True)

    designed_ids: list[str] = []
    interface_ids: list[str] = []
    min_dists: dict[str, float] = {}

    for i in range(len(starts) - 1):
        sl = slice(int(starts[i]), int(starts[i + 1]))
        res = designed[sl]
        rid = residue_mpnn_id(str(res.chain_id[0]), int(res.res_id[0]), str(res.ins_code[0]))
        delta = res.coord[:, None, :] - other_coord[None, :, :]
        dist = float(np.sqrt(np.sum(delta * delta, axis=-1)).min())
        min_dists[rid] = dist
        if dist <= cutoff:
            interface_ids.append(rid)
        else:
            designed_ids.append(rid)

    return {
        "structure": str(path.resolve()),
        "chains": chains,
        "cutoff": cutoff,
        "designed_residues": designed_ids,
        "interface_residues": interface_ids,
        "min_interchain_distance": min_dists,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="List designed-chain residues that are not at a 5 A interchain interface."
    )
    parser.add_argument("structure", help="CIF or PDB path")
    parser.add_argument(
        "--chains",
        default="A",
        help="designed chain IDs, comma-separated (default: A)",
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        default=5.0,
        help="interface cutoff in Angstroms (default: 5.0)",
    )
    parser.add_argument(
        "--json",
        dest="json_path",
        help="write selection JSON here (default: <structure>_interface.json)",
    )
    args = parser.parse_args(argv)
    path = Path(args.structure)
    chains = [c.strip() for c in args.chains.split(",") if c.strip()]
    if not chains:
        print("ERROR: --chains is empty", file=sys.stderr)
        return 1
    try:
        result = select_designed_residues(path, chains, cutoff=args.cutoff)
    except (ValueError, FileNotFoundError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    json_path = (
        Path(args.json_path)
        if args.json_path
        else path.with_name(path.name + "_interface.json")
    )
    json_path.parent.mkdir(parents=True, exist_ok=True)
    Path(json_path).write_text(json.dumps(result, indent=2) + "\n")
    if not result["designed_residues"]:
        print(
            f"ERROR: no non-interface residues on chains {chains} "
            f"(cutoff={args.cutoff} A) in {path}",
            file=sys.stderr,
        )
        return 1
    print(",".join(result["designed_residues"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
