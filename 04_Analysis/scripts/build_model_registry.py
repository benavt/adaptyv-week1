#!/usr/bin/env python3
"""Build additive, source-backed homes for project structures and sequences.

Original scientific files are read only. Generated files live under registry/;
relative symlinks retain the original layout used by existing analysis scripts.
No scientific jobs, scoring services, or external software are run.
"""
from __future__ import annotations

import argparse
import collections
import csv
import datetime as dt
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import sys

SCHEMA = 1
ROOT = Path(__file__).resolve().parents[2]
AA = dict(zip("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(), "ARNDCQEGHILKMFPSTWYV"))
AA.update({"HID": "H", "HIE": "H", "HIP": "H", "HSD": "H", "HSE": "H", "HSP": "H", "ASH": "D", "GLH": "E", "CYX": "C", "MSE": "M"})
EXCLUDED = {".git", ".venv", "venv", "node_modules", "__pycache__", ".vendor", "software", "DockQ", "ProtonPottsMPNN", "ProtonPottsMPNN_runtime", "registry"}
SEQUENCE_SUFFIXES = {".fasta", ".fa", ".faa", ".fna"}


def digest(data):
    return hashlib.sha256(data).hexdigest()


def safe(text):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text).strip("._") or "unknown"


def structure_stem(path):
    name = path.name
    for suffix in (".pdb.original", ".cif.gz", ".pdb.gz", ".cif", ".pdb", ".gro"):
        if name.lower().endswith(suffix):
            return name[:-len(suffix)]
    return path.stem


def is_structure(path):
    return path.name.lower().endswith((".cif", ".pdb", ".cif.gz", ".pdb.gz", ".pdb.original", ".gro"))


def read_json(path):
    return json.loads(path.read_text())


def table(path):
    if not path.is_file():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t" if path.suffix == ".tsv" else ","))


def fasta_records(path):
    header = None
    sequence = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(sequence)
            header, sequence = line[1:].strip(), []
        elif line.strip() and header is not None:
            sequence.append(line.strip())
    if header is not None:
        yield header, "".join(sequence)


def cif_tokens(lines):
    """CIF lexical tokens, including quoted and semicolon-delimited values."""
    iterator = iter(lines)
    for line in iterator:
        if line.startswith(";"):
            value = [line[1:].rstrip("\n")]
            for continuation in iterator:
                if continuation.startswith(";"):
                    break
                value.append(continuation.rstrip("\n"))
            yield "\n".join(value)
        elif line.strip() and not line.lstrip().startswith("#"):
            # Atom tables dominate the large RFD3 collection and contain no quotes.
            if "'" not in line and '"' not in line:
                yield from line.split()
            else:
                lexer = shlex.shlex(line, posix=True)
                lexer.whitespace_split = True
                lexer.commenters = "#"
                yield from lexer


def _parse_structure_full(path):
    raw = path.read_bytes()
    text = (gzip.decompress(raw) if path.name.endswith(".gz") else raw).decode("utf-8", errors="replace")
    residues = collections.defaultdict(dict)
    model_numbers = set()
    first_model = None
    header = {}
    if ".cif" in path.name.lower():
        tokens = iter(cif_tokens(text.splitlines(keepends=True)))
        pending = None
        while True:
            token = pending if pending is not None else next(tokens, None)
            pending = None
            if token is None:
                break
            if token == "loop_":
                columns = []
                token = next(tokens, None)
                while token and token.startswith("_"):
                    columns.append(token)
                    token = next(tokens, None)
                if not columns:
                    raise ValueError("empty CIF loop")
                atom = columns[0].startswith("_atom_site.")
                idx = {k.split(".", 1)[1]: i for i, k in enumerate(columns)}
                row = []
                while token is not None:
                    if not row and (token == "loop_" or token.startswith(("_", "data_", "save_"))):
                        pending = token
                        break
                    row.append(token)
                    if len(row) == len(columns):
                        if atom:
                            def value(*keys, default=""):
                                for key in keys:
                                    if key in idx and row[idx[key]] not in {"?", "."}:
                                        return row[idx[key]]
                                return default
                            model = value("pdbx_PDB_model_num", default="1")
                            model_numbers.add(model)
                            if first_model is None:
                                first_model = model
                            name = value("auth_comp_id", "label_comp_id")
                            if model == first_model and value("auth_atom_id", "label_atom_id") == "CA" and name in AA:
                                chain = value("auth_asym_id", "label_asym_id", default="_")
                                number = value("auth_seq_id", "label_seq_id")
                                insertion = value("pdbx_PDB_ins_code")
                                residues[chain].setdefault((number, insertion), AA[name])
                        row = []
                    token = next(tokens, None)
                if row:
                    raise ValueError("incomplete CIF loop row")
            elif token.startswith("_"):
                value = next(tokens, None)
                if token.startswith(("_entry.", "_struct.pdbx_structure_determination_methodology")):
                    header[token] = value
    else:
        model = "1"
        for line in text.splitlines():
            if line.startswith("MODEL "):
                model = line[10:14].strip() or "1"
            elif line.startswith(("ATOM  ", "HETATM")):
                model_numbers.add(model)
                if first_model is None:
                    first_model = model
                name = line[17:20].strip()
                if model == first_model and line[12:16].strip() == "CA" and name in AA:
                    residues[line[21:22].strip() or "_"].setdefault((line[22:26].strip(), line[26:27].strip()), AA[name])
            elif line.startswith(("TITLE ", "HEADER", "EXPDTA")):
                header.setdefault(line[:6].strip(), []).append(line[6:].strip())
    if not residues:
        raise ValueError("no protein alpha-carbon residues found")
    chains = []
    for chain, values in residues.items():
        chains.append({"chain_id": chain, "sequence": "".join(values.values()), "residue_ids": [n + i for n, i in values], "sequence_basis": "observed alpha-carbon residues of first coordinate model; not assumed complete"})
    return {"sha256": digest(raw), "size_bytes": len(raw), "chains": chains, "coordinate_model_numbers": sorted(model_numbers), "embedded_header": header}


def parse_structure(path):
    """Fast path for ordinary one-row-per-line atom loops; full CIF fallback."""
    if path.suffix.lower() == ".gro":
        raw = path.read_bytes()
        lines = raw.decode().splitlines()
        count = int(lines[1].strip())
        residues = {}
        for line in lines[2:count + 2]:
            name, atom = line[5:10].strip(), line[10:15].strip()
            if atom == "CA" and name in AA:
                residues.setdefault(line[:5].strip(), AA[name])
        if not residues:
            raise ValueError("no protein alpha-carbon residues in GRO")
        return {"sha256": digest(raw), "size_bytes": len(raw), "chains": [{"chain_id": "unrecorded", "sequence": "".join(residues.values()), "residue_ids": list(residues), "sequence_basis": "observed protein alpha-carbon residues; GRO has no chain IDs"}], "coordinate_model_numbers": ["1"], "embedded_header": {"title": lines[0]}}
    if ".cif" not in path.name.lower():
        return _parse_structure_full(path)
    raw = path.read_bytes()
    text = (gzip.decompress(raw) if path.name.endswith(".gz") else raw).decode("utf-8", errors="replace")
    lines = iter(text.splitlines())
    header = {}
    residues = collections.defaultdict(dict)
    numbers = set()
    first = None
    try:
        for line in lines:
            if line.startswith("_entry."):
                bits = shlex.split(line)
                if len(bits) == 2:
                    header[bits[0]] = bits[1]
            if line.strip() != "loop_":
                continue
            columns = []
            current = next(lines, "")
            while current.strip().startswith("_"):
                columns.append(current.strip())
                current = next(lines, "")
            if not columns or not columns[0].startswith("_atom_site."):
                continue
            idx = {k.split(".", 1)[1]: i for i, k in enumerate(columns)}
            def index(*keys):
                return next(idx[k] for k in keys if k in idx)
            atom = index("auth_atom_id", "label_atom_id")
            comp = index("auth_comp_id", "label_comp_id")
            chain = index("auth_asym_id", "label_asym_id")
            residue = index("auth_seq_id", "label_seq_id")
            model = idx.get("pdbx_PDB_model_num")
            insertion = idx.get("pdbx_PDB_ins_code")
            while current and not current.lstrip().startswith(("#", "_", "loop_", "data_")):
                values = shlex.split(current) if "'" in current or '"' in current else current.split()
                if len(values) != len(columns):
                    raise ValueError("wrapped atom row; full parser required")
                number = values[model] if model is not None else "1"
                numbers.add(number)
                if first is None:
                    first = number
                if number == first and values[atom] == "CA" and values[comp] in AA:
                    chn = values[chain] if values[chain] not in {"?", "."} else values[idx["label_asym_id"]]
                    rid = values[residue] if values[residue] not in {"?", "."} else values[idx["label_seq_id"]]
                    ins = values[insertion] if insertion is not None and values[insertion] not in {"?", "."} else ""
                    residues[chn].setdefault((rid, ins), AA[values[comp]])
                current = next(lines, "")
        if not residues:
            raise ValueError("no readable protein atom loop")
        chains = [{"chain_id": c, "sequence": "".join(v.values()), "residue_ids": [n + i for n, i in v], "sequence_basis": "observed alpha-carbon residues of first coordinate model; not assumed complete"} for c, v in residues.items()]
        return {"sha256": digest(raw), "size_bytes": len(raw), "chains": chains, "coordinate_model_numbers": sorted(numbers), "embedded_header": header}
    except (ValueError, KeyError, StopIteration):
        return _parse_structure_full(path)


class Registry:
    def __init__(self, root, output, reuse_coordinates=False):
        self.root, self.output = root.resolve(), output.resolve()
        if not self.output.is_relative_to(self.root) or self.output == self.root:
            raise ValueError("registry output must be a directory inside the project")
        self.now = dt.datetime.now(dt.timezone.utc).isoformat()
        self.models = {}
        self.designs = {}
        self.sequences = {}
        self.sequence_observation_keys = collections.defaultdict(set)
        self.aliases = collections.defaultdict(set)
        self.issues = []
        self.sources = {}
        self.path_model = {}
        self.metrics = collections.defaultdict(list)
        self.targets = collections.defaultdict(set)
        self.refolds = {}
        self.sequence_inputs = []
        self.source_designs = {}
        self.generator_sha256 = digest(Path(__file__).read_bytes())
        self.coordinate_cache = {}
        if reuse_coordinates and (self.output / "models.csv").is_file():
            for row in table(self.output / "models.csv"):
                self.coordinate_cache[row["source_path"]] = self.root / row["home"] / "metadata.json"

    def relative(self, path):
        path = Path(path)
        if path.is_absolute():
            try:
                return path.relative_to(self.root).as_posix()
            except ValueError:
                if path.exists():
                    try:
                        return path.resolve().relative_to(self.root).as_posix()
                    except ValueError:
                        pass
                # Historic root paths are recognizable without guessing by basename.
                marker = "/Week1/"
                if marker in path.as_posix():
                    candidate = path.as_posix().split(marker, 1)[1]
                    if (self.root / candidate).exists():
                        return candidate
                marker = "/foundry/runs/"
                if marker in path.as_posix():
                    candidate = "02_Design/foundry/runs/" + path.as_posix().split(marker, 1)[1]
                    if (self.root / candidate).exists():
                        return candidate
                return path.as_posix()
        return path.as_posix()

    def issue(self, code, path, detail):
        self.issues.append({"code": code, "source_path": self.relative(path), "detail": detail})

    def source(self, path):
        rel = self.relative(path)
        if rel not in self.sources:
            p = self.root / rel
            data = p.read_bytes()
            s = p.stat()
            self.sources[rel] = {"source_path": rel, "sha256": digest(data), "size_bytes": len(data), "mtime_ns": s.st_mtime_ns}
        return self.sources[rel]

    def sequence(self, value, source, role, header=None, basis="recorded sequence"):
        if not value:
            return None
        canonical = value.strip().upper()
        if not re.fullmatch(r"[A-Z]+", canonical):
            self.issue("noncanonical_sequence", source, "Sequence retained as a source artifact, not converted from extended-state tokens.")
            return None
        sid = "SEQ_" + digest(canonical.encode())[:20]
        entry = self.sequences.setdefault(sid, {"sequence_id": sid, "sequence_sha256": digest(canonical.encode()), "sequence": canonical, "length": len(canonical), "observations": []})
        observation = {"source_path": self.relative(source), "role": role, "header": header, "basis": basis}
        observation_key = tuple(sorted(observation.items()))
        if observation_key not in self.sequence_observation_keys[sid]:
            self.sequence_observation_keys[sid].add(observation_key)
            entry["observations"].append(observation)
        return sid

    def design(self, identifier, sequence, source, aliases=(), provenance=None):
        sid = self.sequence(sequence, source, "binder")
        entry = self.designs.setdefault(identifier, {"schema_version": SCHEMA, "design_id": identifier, "sequence_id": sid, "aliases": [], "evidence": [], "model_ids": [], "statistics": [], "provenance": provenance or {}, "warnings": []})
        if sid and entry["sequence_id"] and sid != entry["sequence_id"]:
            self.issue("design_sequence_conflict", source, f"{identifier} has conflicting sequences; conflicting observation not merged.")
            return None
        if sid:
            entry["sequence_id"] = sid
        entry["evidence"].append({"source_path": self.relative(source), "basis": "recorded sequence or sequence extracted from this structure"})
        for alias in (identifier, *aliases):
            self.aliases[alias].add(identifier)
            if alias not in entry["aliases"]:
                entry["aliases"].append(alias)
        return identifier

    def resolve_design(self, alias):
        candidates = self.aliases.get(alias, set())
        return next(iter(candidates)) if len(candidates) == 1 else None

    def load_designs(self):
        # Imported LigandAI designs have run-qualified identities and exact FASTA evidence.
        ligand_manifest = self.root / "02_Design/Ligand_AI/designs.json"
        if ligand_manifest.is_file():
            data = read_json(ligand_manifest)
            fasta = self.root / data["source_fasta"]
            observed = list(fasta_records(fasta))
            self.source(ligand_manifest)
            self.source(fasta)
            if digest(fasta.read_bytes()) != data["source_fasta_sha256"]:
                raise ValueError("LigandAI source FASTA checksum disagrees with its manifest")
            if len(observed) != len(data["designs"]):
                raise ValueError("LigandAI source FASTA record count disagrees with its manifest")
            for row, (header, sequence) in zip(data["designs"], observed):
                if header != row["source_header"] or sequence != row["binder_sequence"] or digest(sequence.encode()) != row["sequence_sha256"]:
                    raise ValueError("LigandAI design does not match its source FASTA record")
                self.design(row["design_id"], sequence, fasta,
                            provenance={"method": data["method"], "generation_run_id": data["generation_run_id"],
                                        "checkpoint": data["checkpoint"], "evidence_path": self.relative(ligand_manifest),
                                        "source_record_number": row["source_record_number"], "source_header": header,
                                        "unknowns": ["exact generation checkpoint and generation command"]})
        report = self.root / "03_Filtering/RFD3_filtering/out/dd1_asp_near_both_his_sequences.csv"
        report_rows = {row["design_id"]: row for row in table(report)}
        if report.exists():
            self.source(report)
        # Adopt aliases from files already staged, never from current CSV ordering.
        for path in sorted((self.root / "01_Staging/Putative_DD1/FASTA").glob("*.fasta")):
            records = list(fasta_records(path))
            match = re.match(r"(DD1_\d+)_(Human|Mouse)_(.+)", path.stem)
            if not match or len(records) != 2:
                self.issue("unresolved_staged_design", path, "Expected existing DD1 alias and paired FASTA.")
                continue
            alias, species, legacy = match.groups()
            row = report_rows.get(legacy)
            if not row or records[0][1] != row["binder_sequence"]:
                self.issue("staged_sequence_conflict", path, "FASTA binder does not match recorded source design; alias not adopted.")
                continue
            self.design(alias, records[0][1], path, aliases=[legacy], provenance={"method": "Foundry RFD3", "evidence_path": self.relative(report), "source_structure": self.relative(row["structure_path"]), "mutation_parent": None})
            self.source_designs[self.relative(row["structure_path"])] = alias
            self.targets[records[1][1]].add(species)
        for path in sorted((self.root / "01_Staging").glob("*species_intersection/*.fasta")):
            match = re.match(r"(DD1_\d+_.+)_(Human|Mouse)_EGFR", path.stem)
            records = list(fasta_records(path))
            if match and len(records) == 2:
                identifier, species = match.groups()
                self.design(identifier, records[0][1], path, provenance={"method": "staged sequence mutation", "evidence_path": self.relative(path.parent / "README.md"), "mutation_parent": identifier[:7], "mutation_labels": identifier[8:].split("_")})
                self.targets[records[1][1]].add(species)
        for species in ("Human", "Mouse"):
            path = self.root / f"01_Staging/EGFR/{species}_EGFR_310_501_afs3.json"
            if path.is_file():
                for item in read_json(path).get("sequences", []):
                    if "protein" in item:
                        self.targets[item["protein"]["sequence"]].add(species)
        # Sequence designs inside a software submodule are explicitly allowlisted.
        potts = self.root / "03_Filtering/ProtonPottsMPNN/inference/outputs"
        for path in (potts / "designs.tsv", potts / "sweep_designs.tsv"):
            for rownum, row in enumerate(table(path), 2):
                sequence = row.get("canonical_sequence")
                if sequence and row.get("design_id"):
                    did = self.design(row["design_id"], sequence, path, provenance={"method": "ProtonPottsMPNN sequence design", "checkpoint": None, "seed": None, "evidence_path": self.relative(path), "unknowns": ["executed command", "exact checkpoint used", "parent structure"]})
                    if did:
                        self.designs[did]["statistics"].append({"source_path": self.relative(path), "row_number": rownum, "values": row})
            if path.is_file():
                self.source(path)
        manifest = potts / "pareto_fold_manifest.json"
        if manifest.is_file():
            self.source(manifest)
            for i, row in enumerate(read_json(manifest)):
                did = self.design(row["design_id"], row["canonical_sequence"], manifest, provenance={"method": "ProtonPottsMPNN sequence design", "fold_status": "fold input manifest exists; successful fold not established", "evidence_path": self.relative(manifest)})
                if did:
                    self.designs[did]["statistics"].append({"source_path": self.relative(manifest), "json_pointer": f"/{i}", "values": row})
        full = potts / "designs.json"
        if full.is_file():
            self.source(full)
            for i, row in enumerate(read_json(full)):
                did = self.resolve_design(row.get("design_id"))
                if did:
                    entry = self.designs[did]
                    entry["provenance"]["design_configuration"] = {k: v for k, v in row.items() if k not in {"energy_trajectory", "extended_tokens"}}
                    entry["provenance"]["raw_json"] = {"source_path": self.relative(full), "json_pointer": f"/{i}"}
        for did, entry in list(self.designs.items()):
            parent = entry["provenance"].get("mutation_parent")
            if parent and parent in self.designs:
                a = self.sequences[self.designs[parent]["sequence_id"]]["sequence"]
                b = self.sequences[entry["sequence_id"]]["sequence"]
                observed = [f"{x}{i}{y}" for i, (x, y) in enumerate(zip(a, b), 1) if x != y]
                expected = entry["provenance"]["mutation_labels"]
                verified = len(a) == len(b) and set(observed) == set(expected)
                entry["provenance"].update({"observed_mutations": observed, "mutations_verified": verified})
                if not verified:
                    self.issue("mutation_label_conflict", entry["evidence"][0]["source_path"], f"{did}: labels {expected}; observed {observed}")

    def load_statistics(self):
        # Keep raw values and row locators: do not silently rescale or aggregate.
        sources = [
            "03_Filtering/RFD3_filtering/out/summary.csv",
            "03_Filtering/RFD3_filtering/out/dd1_asp_his_proximity.csv",
            "03_Filtering/RFD3_filtering/out/dd1_asp_his_proximity_684_refolded.csv",
            "03_Filtering/Boltz2_metrics/out/boltz2_dd1_metrics_long.csv",
            "03_Filtering/PyMol_top20/out/dockq_tm_all_designs.csv",
            "03_Filtering/PyMol_top20/out/mouse_vs_human_tm_scores.csv",
        ]
        for name in sources:
            path = self.root / name
            if not path.is_file():
                continue
            self.source(path)
            for number, row in enumerate(table(path), 2):
                paths = [self.relative(row[k]) for k in ("structure_path", "source_cif", "human_cif", "mouse_cif") if row.get(k)]
                if row.get("boltz_dir"):
                    folder = self.root / self.relative(row["boltz_dir"])
                    paths += [self.relative(p) for p in folder.glob("*_model_0.cif")]
                    for rel in paths:
                        self.refolds[rel] = row
                for rel in paths:
                    self.metrics[rel].append({"source_path": name, "row_number": number, "association": "explicit structure path in source table", "values": row})
        for path in sorted((self.root / "03_Filtering/DeltaForge/out").glob("*/*results.json")):
            data = read_json(path)
            if not isinstance(data, list):
                continue
            self.source(path)
            for i, row in enumerate(data):
                if not isinstance(row, dict):
                    continue
                rel = self.relative(row["source_cif"]) if row.get("source_cif") else None
                species = row.get("species")
                if not rel and path.parent.name.startswith("DD1_") and species:
                    folder = self.root / f"01_Staging/{path.parent.name}_EGFR_species_intersection/predictions/{path.parent.name}_{species}_EGFR"
                    matches = sorted(folder.glob("*.cif"))
                    if len(matches) == 1:
                        rel = self.relative(matches[0])
                if not rel and row.get("design_id") and species:
                    did = self.resolve_design(row["design_id"])
                    matches = [p for p, r in self.refolds.items() if self.resolve_design(r["design_id"]) == did and r.get("target") == species] if did else []
                    if len(matches) == 1:
                        rel = matches[0]
                item = {"source_path": self.relative(path), "json_pointer": f"/{i}", "association": "recorded source_cif" if row.get("source_cif") else "source script / design and species lookup; check model metadata warnings", "values": row}
                if rel:
                    self.metrics[rel].append(item)
                else:
                    self.issue("unassigned_statistics", path, f"JSON row {i}: exact source structure unresolved")
        # Selections are design observations, never model identities.
        for name in ("03_Filtering/Boltz2_metrics/out/top20_grid_search_candidates.csv", "03_Filtering/Boltz2_metrics/out/boltz2_dd1_metrics_paired.csv", "03_Filtering/JustHIpKA/top20_DD1_EGFR_D3/top20_hispka_results.csv"):
            path = self.root / name
            if path.is_file():
                self.source(path)
                for i, row in enumerate(table(path), 2):
                    did = self.resolve_design(row.get("design_id"))
                    if did:
                        self.designs[did]["statistics"].append({"source_path": name, "row_number": i, "scope": "design / selection; not automatically a structure score", "values": row})

    def discover(self):
        files = []
        for base, directories, names in os.walk(self.root, followlinks=False):
            directories[:] = sorted(d for d in directories if d not in EXCLUDED and not (Path(base) / d).is_symlink() and not (Path(base) / d).resolve().is_relative_to(self.output) and (Path(base) / d).relative_to(self.root).as_posix() not in {"tools/environments", "tools/cache"})
            for name in sorted(names):
                path = Path(base) / name
                if not path.is_symlink() and (is_structure(path) or path.suffix.lower() in SEQUENCE_SUFFIXES):
                    files.append(path)
        potts = self.root / "03_Filtering/ProtonPottsMPNN/inference/outputs"
        files.extend(p for p in sorted(potts.iterdir()) if p.is_file() and (is_structure(p) or p.suffix in SEQUENCE_SUFFIXES)) if potts.exists() else None
        return sorted(set(files))

    def classify(self, rel):
        if rel.startswith("02_Design/foundry/runs/"):
            run = rel.split("/")[3]
            if "/rfd3/outputs/" in rel:
                return "foundry", run, "Foundry RFD3", "per-model JSON when available"
            return "foundry_inputs", run, None, "staged/prepared input, not proof of new design generation"
        if rel.startswith("03_Filtering/Refolding/"):
            return "refolds", "putative_DD1", "Boltz2", "documented by Boltz2_metrics/README.md and analyzer; exact execution unavailable"
        if "species_intersection/predictions/" in rel:
            return "refolds", "mutants", "Boltz2", "confidence schema and scoring script; inferred tool, exact execution unavailable"
        if rel.startswith("03_Filtering/JustHIpKA/"):
            return "preparations", "justhispka", "structure preparation for JustHISpKa", "source preparation scripts / path; individual invocation unrecorded"
        if rel.startswith("03_Filtering/DeltaForge/"):
            return "preparations", "deltaforge", "CIF-to-PDB scoring preparation", "source runner / path; individual invocation unrecorded"
        if rel.startswith("03_Filtering/PyMol_top20/"):
            return "preparations", "comparison", "structural comparison preparation", "analysis scripts / path; individual invocation unrecorded"
        if rel.startswith("03_Filtering/OpenMM_mini/"):
            return "preparations", "openmm", "PDBFixer preparation / OpenMM energy minimization", "OpenMM_mini/README.md and matching output naming; exact invocation unrecorded"
        if rel.startswith("01_Staging/MD/"):
            return "simulations", "EGFR_D3", "GROMACS preparation / trajectory extraction", "MD README and source scripts; frame timing not inferred from filename"
        if rel.startswith("01_Staging/frozen-tiger-ice/"):
            return "imports", "frozen-tiger-ice", "AlphaFold 3", "imported package data JSON and CIF attribution"
        if rel.startswith("01_Staging/EGFR/"):
            return "references", "EGFR", None, "reference / prepared structure; exact derivation requires source evidence"
        return "unresolved", "other", None, "origin unresolved"

    def sidecars(self, path, group):
        stem = structure_stem(path)
        files = []
        if group == "foundry":
            p = path.parent / (stem + ".json")
            if p.is_file():
                files.append(p)
        elif group in {"refolds", "imports"}:
            for p in sorted(path.parent.iterdir()):
                if p.is_file() and p != path and not is_structure(p) and p.suffix in {".json", ".npz", ".txt"}:
                    if group == "refolds" or p.name.startswith(stem.replace("_model", "")) or path.parent.name.startswith("seed-"):
                        files.append(p)
        return files

    def build_model(self, path):
        rel = self.relative(path)
        group, run, method, basis = self.classify(rel)
        structure_status = "indexed"
        try:
            cached_path = self.coordinate_cache.get(rel)
            cached = read_json(cached_path) if cached_path and cached_path.is_file() else None
            raw = path.read_bytes() if cached else None
            if cached and cached.get("structure_status") == "indexed" and digest(raw) == cached["source_sha256"]:
                record = {"sha256": cached["source_sha256"], "size_bytes": len(raw), "chains": [{k: c[k] for k in ("chain_id", "sequence", "residue_ids", "sequence_basis")} for c in cached["chains"]], "coordinate_model_numbers": cached["coordinate_model_numbers"], "embedded_header": cached["embedded_header"]}
            else:
                record = parse_structure(path)
        except (ValueError, IndexError, KeyError) as exc:
            raw = path.read_bytes()
            structure_status = "empty_file" if not raw.strip() else "unparsed_coordinates"
            pdb_header = [line for line in raw.decode(errors="replace").splitlines() if line.startswith(("HEADER", "TITLE", "EXPDTA", "REMARK"))]
            record = {"sha256": digest(raw), "size_bytes": len(raw), "chains": [], "coordinate_model_numbers": [], "embedded_header": {"pdb_header_lines": pdb_header}}
            self.issue("empty_structure_artifact" if structure_status == "empty_file" else "structure_parse_error", path, str(exc))
        mid = "MOD_" + digest(rel.encode())[:20]
        home = self.output / "models" / group / safe(run) / (safe(structure_stem(path)) + "__" + mid[4:12])
        provenance = {"method": method, "method_basis": basis, "executed_command": None, "tool_version": None, "checkpoint": None, "seed": None, "created_at": None, "evidence": [], "input_paths": [], "parent_model_ids": []}
        meta = {"schema_version": SCHEMA, "model_id": mid, "home": self.relative(home), "source_path": rel, "source_sha256": record["sha256"], "source_size_bytes": record["size_bytes"], "structure_status": structure_status, "artifact_kind": group, "run_id": run, "design_id": None, "target_species": None, "chains": record["chains"], "coordinate_model_numbers": record["coordinate_model_numbers"], "embedded_header": record["embedded_header"], "provenance": provenance, "statistics": self.metrics.get(rel, []).copy(), "related_artifacts": [], "warnings": [] if structure_status == "indexed" else ["No protein coordinates indexed; this is an empty or unreadable artifact, not a valid predicted model."]}
        header_lines = record["embedded_header"].get("pdb_header_lines", [])
        if any("MODEL GENERATED BY ROSETTA" in line for line in header_lines):
            provenance["method"] = "Rosetta"
            provenance["method_basis"] = "explicit source PDB REMARK header"
            version_line = next((line for line in header_lines if "VERSION" in line), None)
            provenance["tool_version"] = version_line.split("VERSION", 1)[1].strip() if version_line else None
        self.sources[rel] = {"source_path": rel, "sha256": record["sha256"], "size_bytes": record["size_bytes"], "mtime_ns": path.stat().st_mtime_ns}
        sidecars = self.sidecars(path, group)
        for p in sidecars:
            meta["related_artifacts"].append({"source_path": self.relative(p), "relationship": "same output bundle"})
            if p.suffix == ".json":
                self.source(p)
                payload = read_json(p)
                if group == "foundry":
                    provenance.update({"checkpoint": payload.get("ckpt_path"), "seed": payload.get("seed"), "configuration": payload.get("specification"), "evidence": [{"source_path": self.relative(p), "basis": "per-model output JSON"}]})
                    if payload.get("specification", {}).get("input"):
                        provenance["input_paths"].append(payload["specification"]["input"])
                    meta["statistics"].append({"source_path": self.relative(p), "json_pointer": "/metrics", "association": "same-stem model output", "values": payload.get("metrics", {})})
                elif p.name.startswith("confidence") or "summary_confidences" in p.name:
                    meta["statistics"].append({"source_path": self.relative(p), "json_pointer": "", "association": "same prediction bundle; original confidence units retained", "values": payload})
        if group == "foundry":
            manifest = self.root / f"02_Design/foundry/runs/{run}/manifest.json"
            if manifest.is_file():
                self.source(manifest)
                provenance["run_manifest"] = self.relative(manifest)
                meta["related_artifacts"].append({"source_path": self.relative(manifest), "relationship": "run manifest"})
                provenance["evidence"].append({"source_path": self.relative(manifest), "basis": "run manifest; lists intended pipeline stages, not proof they all ran"})
            local_inputs = sorted(path.parents[1].glob("*.pdb")) + sorted(path.parents[1].glob("*.yaml"))
            meta["related_artifacts"].extend({"source_path": self.relative(p), "relationship": "run input"} for p in local_inputs)
            header = record["embedded_header"]
            if header.get("_entry.date"):
                provenance["created_at"] = {"date": header.get("_entry.date"), "time": header.get("_entry.time"), "timezone": None, "basis": "embedded CIF header; timezone unrecorded"}
        targets = []
        binders = []
        for chain in meta["chains"]:
            species = self.targets.get(chain["sequence"], set())
            chain["role"] = "EGFR_D3" if species else "unassigned"
            chain["role_basis"] = "exact sequence match to staged EGFR D3" if species else "unknown"
            chain["sequence_id"] = self.sequence(chain["sequence"], path, chain["role"], header=f"chain {chain['chain_id']}", basis=chain["sequence_basis"])
            if len(species) == 1:
                targets.append(next(iter(species)))
            if not species:
                binders.append(chain)
        if len(set(targets)) == 1:
            meta["target_species"] = targets[0]
        if group == "foundry" and len(binders) == 1 and targets:
            binder = binders[0]
            legacy = structure_stem(path)
            did = self.source_designs.get(rel)
            if did and self.sequences[self.designs[did]["sequence_id"]]["sequence"] != binder["sequence"]:
                self.issue("structure_design_sequence_conflict", path, f"{did} does not match actual binder sequence")
                did = None
            if not did:
                did = self.design(run + "--" + legacy, binder["sequence"], path, aliases=[run + "::" + legacy], provenance={"method": "Foundry RFD3", "source_model_id": mid, "original_model_name": legacy, "evidence_path": rel})
            meta["design_id"] = did
            binder.update({"role": "binder", "role_basis": "only non-target protein chain in RFD3 output"})
        elif group == "refolds" and len(binders) == 1 and targets:
            binder = binders[0]
            metric_row = self.refolds.get(rel)
            if metric_row:
                did = self.resolve_design(metric_row["design_id"])
            else:
                match = re.match(r"(DD1_\d+(?:_[A-Z]\d+[A-Z])*)_(Human|Mouse)_", path.name)
                did = self.resolve_design(match[1]) if match else None
            if did and self.sequences[self.designs[did]["sequence_id"]]["sequence"] == binder["sequence"]:
                meta["design_id"] = did
                binder.update({"role": "binder", "role_basis": "exact sequence match to staged design"})
                if metric_row and metric_row["target"] != meta["target_species"]:
                    meta["warnings"].append("Target in metric table disagrees with actual structure sequence.")
                    self.issue("statistics_species_conflict", path, "Recorded target differs from sequence-derived target")
            else:
                self.issue("unresolved_refold_design", path, "Actual binder sequence does not uniquely validate the recorded design alias")
            match = re.match(r"DD1_\d+_(Human|Mouse)_", path.name)
            if match and match[1] != meta["target_species"]:
                meta["warnings"].append(f"Filename says {match[1]}; actual target sequence is {meta['target_species']}. Filename is preserved, not used as target identity.")
        elif group == "foundry_inputs" and len(binders) == 1:
            binder = binders[0]
            candidates = [d for d, v in self.designs.items() if v["sequence_id"] == binder["sequence_id"]]
            if len(candidates) == 1:
                meta["design_id"] = candidates[0]
                binder.update({"role": "binder", "role_basis": "unique sequence match; original design identity retained"})
            elif "DD2_top" in run:
                did = self.design(run, binder["sequence"], path, provenance={"method": "imported/staged DD2 seed binder", "original_design_method": None, "evidence_path": "03_Filtering/ProtonPottsMPNN_DD2_smoke_test.md", "note": "Hotspot display inputs are not newly generated RFD3 designs."})
                meta["design_id"] = did
                binder.update({"role": "binder", "role_basis": "DD2 seed input documented in smoke-test plan"})
        elif group == "imports" and run == "frozen-tiger-ice" and len(binders) == 1 and targets:
            binder = binders[0]
            meta["design_id"] = self.design("frozen-tiger-ice", binder["sequence"], path, provenance={"method": "imported binder sequence", "sequence_design_method": None, "fold_method": "AlphaFold 3", "evidence_path": "01_Staging/frozen-tiger-ice/frozen-tiger-ice_data.json"})
            binder.update({"role": "binder", "role_basis": "only non-EGFR protein chain in imported complex"})
            input_json = self.root / "01_Staging/frozen-tiger-ice/frozen-tiger-ice_data.json"
            if input_json.is_file():
                self.source(input_json)
                provenance["input_paths"].append(self.relative(input_json))
                provenance["evidence"].append({"source_path": self.relative(input_json), "basis": "imported AlphaFold input; includes requested seeds"})
            seed = re.search(r"seed-(\d+)_sample-(\d+)", rel)
            if seed:
                provenance["seed"] = {"seed": int(seed[1]), "sample": int(seed[2]), "basis": "imported output directory label"}
        if len(meta["coordinate_model_numbers"]) > 1:
            meta["warnings"].append("Ensemble file: sequences indexed from first coordinate model; all coordinate models remain in original file.")
        provenance["unknowns"] = [k for k in ("method", "executed_command", "tool_version", "checkpoint", "seed", "created_at") if provenance.get(k) is None]
        for chain in meta["chains"]:
            if chain["role"] == "binder":
                self.sequence(chain["sequence"], path, "binder", header=f"chain {chain['chain_id']}", basis=chain["role_basis"])
        self.models[mid] = meta
        self.path_model[rel] = mid
        if meta["design_id"]:
            self.designs[meta["design_id"]]["model_ids"].append(mid)

    def attach_relationships(self):
        # Sequence equality alone is not derivation: retain candidate relationships.
        by_pair = collections.defaultdict(list)
        by_sequence = collections.defaultdict(list)
        for did, design in self.designs.items():
            if design["sequence_id"]:
                by_sequence[design["sequence_id"]].append(did)
        for mid, model in self.models.items():
            if model["artifact_kind"] == "refolds" and model["design_id"]:
                by_pair[(model["design_id"], model["target_species"])].append(mid)
        for mid, model in self.models.items():
            if model["artifact_kind"] == "refolds" and model["design_id"]:
                design = self.designs[model["design_id"]]
                source = design["provenance"].get("source_structure")
                if source in self.path_model:
                    model["provenance"]["parent_model_ids"].append({"model_id": self.path_model[source], "relationship": "binder design source; refolded geometry independently predicted", "basis": "staged design table and verified binder sequence"})
                for observation in design["evidence"]:
                    p = self.root / observation["source_path"]
                    if p.suffix == ".fasta":
                        pairs = list(fasta_records(p))
                        if len(pairs) == 2 and model["target_species"] in self.targets.get(pairs[1][1], set()):
                            model["related_artifacts"].append({"source_path": self.relative(p), "relationship": "sequence-matched staged input; execution of this exact file not established"})
            if model["artifact_kind"] == "preparations":
                rel = model["source_path"]
                binder_chains = [c for c in model["chains"] if c["role"] != "EGFR_D3"]
                if model["run_id"] in {"deltaforge", "comparison"} and len(binder_chains) == 1 and model["target_species"]:
                    chain = binder_chains[0]
                    candidates = by_sequence.get(chain["sequence_id"], [])
                    if len(candidates) == 1:
                        did = candidates[0]
                        model["design_id"] = did
                        chain.update({"role": "binder", "role_basis": "unique exact sequence match to a registered design"})
                        self.designs[did]["model_ids"].append(mid)
                        parents = by_pair.get((did, model["target_species"]), [])
                        if len(parents) == 1:
                            model["provenance"]["parent_model_ids"].append({"model_id": parents[0], "relationship": "candidate source for scoring/comparison conversion", "basis": "preparation script conventions and exact binder/target sequence match; individual conversion invocation unrecorded"})
                            model["warnings"].append("Conversion source relationship reconstructed; exact invocation not available. Parent statistics remain in the parent's home.")
                if "/JustHIpKA/" in rel:
                    parts = Path(rel).parts
                    if len(parts) < 4:
                        continue
                    runfolder = parts[3] if parts[2] in {"putative_DD1_EGFR_D3", "top20_DD1_EGFR_D3", "intersection_both_species"} else parts[2]
                    name = re.sub(r"^\d+_", "", runfolder)
                    species = model["target_species"]
                    name = re.sub(r"_(Human|Mouse)$", "", name)
                    did = self.resolve_design(name)
                    candidates = by_pair.get((did, species), [])
                    if did and len(candidates) == 1:
                        parent = candidates[0]
                        model["design_id"] = did
                        model["provenance"]["parent_model_ids"].append({"model_id": parent, "relationship": "EGFR chain extraction / preparation", "basis": "preparation script and directory identity; target sequence verified, binder is absent"})
                        self.designs[did]["model_ids"].append(mid)
                        directory = self.root / ("/".join(parts[:4]) if parts[2] in {"putative_DD1_EGFR_D3", "top20_DD1_EGFR_D3", "intersection_both_species"} else "/".join(parts[:3]) + "/" + str(species))
                        for p in sorted(directory.rglob("HIS*.txt")):
                            raw = p.read_text(errors="replace")
                            values = re.findall(r"PKA=([-+\d.eE]+).*?SID=(\d+)", raw)
                            if values:
                                self.source(p)
                                item = {"source_path": self.relative(p), "association": "target-only preparation package; not a binder pKa", "values": [{"pKa": float(v), "residue": int(n), "chain": "A"} for v, n in values]}
                                model["statistics"].append(item)
                                if item not in self.models[parent]["statistics"]:
                                    self.models[parent]["statistics"].append(item)
                                    self.models[parent]["related_artifacts"].append({"source_path": self.relative(directory), "relationship": "JustHISpKa preparation and results package"})
            for p in model["provenance"]["input_paths"]:
                local = self.relative(p)
                if local in self.path_model:
                    model["provenance"]["parent_model_ids"].append({"model_id": self.path_model[local], "relationship": "recorded structure input", "basis": "per-model specification"})
        for model in self.models.values():
            for observation in model["statistics"]:
                declared = observation.get("values", {}).get("species") if isinstance(observation.get("values"), dict) else None
                if declared and model["target_species"] and declared.capitalize() != model["target_species"]:
                    model["warnings"].append("A linked statistics row has a conflicting species label; original observation retained for review.")

    def ingest_sequences(self, paths):
        for path in paths:
            if path.suffix.lower() not in SEQUENCE_SUFFIXES:
                continue
            records = list(fasta_records(path))
            self.source(path)
            for i, (header, value) in enumerate(records, 1):
                sid = self.sequence(value, path, "unassigned", header)
                self.sequence_inputs.append({"source_path": self.relative(path), "record_number": i, "header": header, "sequence_id": sid, "status": "indexed" if sid else "extended_or_unresolved"})

    def write_json(self, path, data):
        path.parent.mkdir(parents=True, exist_ok=True)
        value = json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + "\n"
        if path.exists() and path.read_text() == value:
            return
        temporary = path.with_name(path.name + ".tmp")
        temporary.write_text(value)
        temporary.replace(path)

    def link(self, directory, name, source):
        source = self.root / source
        if not source.exists():
            self.issue("missing_related_artifact", source, "Referenced artifact unavailable")
            return
        directory.mkdir(parents=True, exist_ok=True)
        dest = directory / name
        target = os.path.relpath(source, directory)
        if dest.is_symlink():
            if os.readlink(dest) == target:
                return
            # Never silently retarget user links.
            raise ValueError(f"Existing link has a different target: {dest}")
        if dest.exists():
            raise ValueError(f"Refusing to replace existing file: {dest}")
        dest.symlink_to(target, target_is_directory=source.is_dir())

    def write_csv(self, name, rows, columns):
        path = self.output / name
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns)
            writer.writeheader()
            writer.writerows(rows)

    def write(self):
        marker = self.output / "registry_manifest.json"
        if self.output.exists() and any(self.output.iterdir()) and not marker.exists():
            raise ValueError("Nonempty output is not a generated registry; refusing to overwrite it")
        self.output.mkdir(parents=True, exist_ok=True)
        self.write_json(marker, {"schema_version": SCHEMA, "status": "building", "generated_at_utc": self.now})
        annotation_path = self.output / "annotations.json"
        annotations = read_json(annotation_path) if annotation_path.is_file() else {"models": {}, "designs": {}}
        duplicate = collections.defaultdict(list)
        for mid, model in self.models.items():
            if model["structure_status"] == "indexed":
                duplicate[model["source_sha256"]].append(mid)
        rows = []
        for mid, model in self.models.items():
            home = self.root / model["home"]
            model["annotations"] = annotations.get("models", {}).get(mid, {})
            self.link(home, Path(model["source_path"]).name, model["source_path"])
            model["byte_identical_model_ids"] = [x for x in duplicate[model["source_sha256"]] if x != mid]
            for artifact in model["related_artifacts"]:
                rel = artifact["source_path"]
                self.link(home / "artifacts", digest(rel.encode())[:8] + "__" + Path(rel).name, rel)
            for source in sorted({s["source_path"] for s in model["statistics"]}):
                self.link(home / "statistics_sources", digest(source.encode())[:8] + "__" + Path(source).name, source)
            self.write_json(home / "metadata.json", {k: v for k, v in model.items() if k != "statistics"})
            self.write_json(home / "statistics.json", {"model_id": mid, "observations": model["statistics"], "policy": "Source values retained unchanged; units and scope belong to each source. Different reruns are separate observations."})
            sequence_text = "".join(f">{mid}|chain={c['chain_id']}|role={c['role']}|sequence_id={c['sequence_id']}|observed_coordinates\n{c['sequence']}\n" for c in model["chains"])
            (home / "sequences.fasta").write_text(sequence_text)
            rows.append({"model_id": mid, "design_id": model["design_id"], "artifact_kind": model["artifact_kind"], "run_id": model["run_id"], "target_species": model["target_species"], "method": model["provenance"]["method"], "structure_status": model["structure_status"], "source_path": model["source_path"], "source_sha256": model["source_sha256"], "home": model["home"], "statistics_observations": len(model["statistics"]), "warning_count": len(model["warnings"]), "unknowns": ";".join(model["provenance"]["unknowns"])})
        self.write_csv("models.csv", rows, list(rows[0]) if rows else ["model_id"])
        design_rows = []
        for did, design in sorted(self.designs.items()):
            home = self.output / "designs" / safe(did)
            design["home"] = self.relative(home)
            design["annotations"] = annotations.get("designs", {}).get(did, {})
            design["model_ids"] = sorted(set(design["model_ids"]))
            self.write_json(home / "metadata.json", design)
            if design["sequence_id"]:
                sequence = self.sequences[design["sequence_id"]]["sequence"]
                (home / "sequence.fasta").write_text(f">{did}|{design['sequence_id']}\n{sequence}\n")
            for mid in design["model_ids"]:
                self.link(home / "models", mid, self.models[mid]["home"])
            design_rows.append({"design_id": did, "sequence_id": design["sequence_id"], "aliases": ";".join(design["aliases"]), "model_count": len(design["model_ids"]), "method": design["provenance"].get("method"), "home": design["home"]})
        self.write_csv("designs.csv", design_rows, ["design_id", "sequence_id", "aliases", "model_count", "method", "home"])
        with (self.output / "sequences.jsonl").open("w") as handle:
            for sid, entry in sorted(self.sequences.items()):
                handle.write(json.dumps(entry, sort_keys=True) + "\n")
        self.write_csv("sequences.csv", [{"sequence_id": k, "sequence_sha256": v["sequence_sha256"], "length": v["length"], "sequence": v["sequence"], "observation_count": len(v["observations"])} for k, v in sorted(self.sequences.items())], ["sequence_id", "sequence_sha256", "length", "sequence", "observation_count"])
        self.write_csv("sequence_inputs.csv", self.sequence_inputs, ["source_path", "record_number", "header", "sequence_id", "status"])
        self.write_csv("source_files.csv", list(self.sources.values()), ["source_path", "sha256", "size_bytes", "mtime_ns"])
        self.write_json(self.output / "duplicates.json", {sha: ids for sha, ids in duplicate.items() if len(ids) > 1})
        # Run metadata is separate from structure-specific facts.
        run_rows = []
        for folder in sorted((self.root / "02_Design/foundry/runs").iterdir()):
            if not folder.is_dir():
                continue
            manifest = folder / "manifest.json"
            info = {"run_id": folder.name, "manifest_source": self.relative(manifest) if manifest.is_file() else None, "manifest": read_json(manifest) if manifest.is_file() else None, "unknowns": [] if manifest.is_file() else ["run manifest unavailable"], "models": [m for m, v in self.models.items() if v["run_id"] == folder.name]}
            home = self.output / "runs" / safe(folder.name)
            self.write_json(home / "metadata.json", info)
            self.link(home, "original_run", self.relative(folder))
            run_rows.append({"run_id": folder.name, "manifest_available": manifest.is_file(), "model_count": len(info["models"]), "home": self.relative(home)})
        for species in ("Human", "Mouse"):
            folder = self.root / f"01_Staging/MD/EGFR_D3/{species}"
            if folder.is_dir():
                rid = f"MD_EGFR_D3_{species}"
                home = self.output / "runs" / rid
                artifacts = [self.relative(p) for p in sorted(folder.rglob("*")) if p.is_file() and p.suffix in {".xtc", ".tpr", ".mdp", ".log", ".xvg", ".csv"}]
                self.write_json(home / "metadata.json", {"run_id": rid, "method": "GROMACS", "method_basis": "01_Staging/MD/EGFR_D3/README.md", "target_species": species, "source_directory": self.relative(folder), "artifacts": artifacts, "model_ids": [m for m, v in self.models.items() if v["source_path"].startswith(self.relative(folder) + "/")], "unknowns": ["exact execution command", "absolute generation time; filesystem time not treated as creation time"], "note": "Original logs and MDP parameters linked as evidence; simulation condition and frame times must be read from those sources."})
                self.link(home, "original_run", self.relative(folder))
                run_rows.append({"run_id": rid, "manifest_available": False, "model_count": sum(v["source_path"].startswith(self.relative(folder) + "/") for v in self.models.values()), "home": self.relative(home)})
        self.write_csv("runs.csv", run_rows, ["run_id", "manifest_available", "model_count", "home"])
        self.write_json(self.output / "issues.json", self.issues)
        summary = {"schema_version": SCHEMA, "status": "complete", "generated_at_utc": self.now, "generator": self.relative(Path(__file__)), "generator_sha256": self.generator_sha256, "scope": "Project PDB/CIF/GRO structures and FASTA files, plus allowlisted ProtonPottsMPNN inference outputs. Software fixtures, installed dependencies, weights, and other coordinate formats excluded; preparation packages remain linked.", "source_policy": "Originals retained in place; model/design homes use project-relative links. Rebuild refreshes generated metadata, not scientific data.", "models": len(self.models), "designs": len(self.designs), "unique_sequences": len(self.sequences), "models_with_statistics": sum(bool(v["statistics"]) for v in self.models.values()), "models_by_kind": dict(collections.Counter(v["artifact_kind"] for v in self.models.values())), "structure_status_counts": dict(collections.Counter(v["structure_status"] for v in self.models.values())), "issues_by_code": dict(collections.Counter(v["code"] for v in self.issues)), "duplicate_groups": sum(len(v) > 1 for v in duplicate.values())}
        self.write_json(marker, summary)
        print(json.dumps(summary, indent=2), flush=True)

    def build(self):
        self.load_designs()
        self.load_statistics()
        files = self.discover()
        self.ingest_sequences(files)
        structures = [p for p in files if is_structure(p)]
        # Generation records precede input preparations when assigning aliases.
        structures.sort(key=lambda p: (0 if "/rfd3/outputs/" in p.as_posix() else 1, p.as_posix()))
        for i, path in enumerate(structures, 1):
            try:
                self.build_model(path)
            except (OSError, ValueError, KeyError, IndexError) as exc:
                self.issue("structure_parse_error", path, f"{type(exc).__name__}: {exc}")
            if i % 1000 == 0:
                print(f"Indexed {i}/{len(structures)} structure files", flush=True)
        self.attach_relationships()
        self.write()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--output", type=Path, help="Output directory inside project (default registry/)")
    parser.add_argument("--reuse-coordinate-cache", action="store_true", help="Reuse coordinate sequences only when the original file SHA-256 still matches; provenance/statistics are rebuilt")
    args = parser.parse_args()
    Registry(args.root, args.output or args.root / "registry", args.reuse_coordinate_cache).build()


if __name__ == "__main__":
    main()
