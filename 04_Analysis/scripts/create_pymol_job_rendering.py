#!/usr/bin/env python3
"""Render a unique PyMOL session and annotated PNG for each RFD3 run YAML.

Example:
  python create_pymol_job_rendering.py foundry/runs \\
    --pdb foundry/runs/Brd4ET_hot_7/Brd4ET_closed_state_fix.pdb
"""

from __future__ import annotations

import argparse
import os
import re
import sys
import tempfile
import textwrap
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import yaml

SNIPPET_ROOT = Path(__file__).resolve().parents[1] / "pymol"
if str(SNIPPET_ROOT) not in sys.path:
    sys.path.insert(0, str(SNIPPET_ROOT))

from snippets.pymol_render import (  # noqa: E402
    render_pymol_panel_with_view,
    resolve_view,
    run_pymol_command,
    view_command,
)

HOTSPOT_KEY_RE = re.compile(r"^([A-Za-z]+)(-?\d+)$")
CONTIG_CHAIN_RE = re.compile(r"^([A-Za-z]+)\d+")
CONTIG_BINDER_RE = re.compile(r"^(\d+)-(\d+)$")
CONTIG_TARGET_RE = re.compile(r"^([A-Za-z]+\d+(?:-\d+)?)$")
YAML_SUFFIXES = {".yaml", ".yml"}
NS1A_RUN_RE = re.compile(r"^NS1A_([A-Za-z]+)_20000$")
# PyMOL named-color RGB, scaled to 0-255 for the legend swatch.
PYMOL_RGB = {
    "gray": (128, 128, 128),
    "red": (255, 0, 0),
    "chartreuse": (128, 255, 0),
    "cyan": (0, 255, 255),
    "magenta": (255, 0, 255),
    "marine": (0, 128, 255),
    "orange": (255, 128, 0),
    "yelloworange": (255, 192, 64),
}
KV_RE = re.compile(r"^([A-Za-z_][A-Za-z0-9_]*)=(.*)$")
N_DESIGNS_RE = re.compile(
    r"^n_designs=(\d+)\s*->\s*n_batches=(\d+)(?:\s*\((.+)\))?\s*$"
)
RF3_FOLD_KV_RE = re.compile(r"([A-Za-z_][A-Za-z0-9_]*)=([^\s]+)")
SLURM_HEADER_STOP_PREFIXES = (
    "=====",
    "Saving script",
    "Streaming ",
    "mpnn --",
    "rf3 ",
)


@dataclass(frozen=True)
class ResidueRef:
    chain: str
    resi: int

    def label(self) -> str:
        return f"{self.chain}{self.resi}"


@dataclass
class JobSpec:
    yaml_path: Path
    run_dir: Path
    run_name: str
    target: str
    contig: str | None
    is_non_loopy: Any
    hotspots: dict[str, Any]
    fixed_atoms: dict[str, Any] | None

    @property
    def stem(self) -> str:
        return f"{self.run_name}_{self.target}"

    @property
    def display_dir(self) -> Path:
        return self.run_dir / "run_display"

    @property
    def pml_path(self) -> Path:
        return self.display_dir / f"{self.stem}.pml"

    @property
    def pse_path(self) -> Path:
        return self.display_dir / f"{self.stem}.pse"

    @property
    def pymol_png_path(self) -> Path:
        return self.display_dir / f"{self.stem}_pymol.png"

    @property
    def composite_png_path(self) -> Path:
        return self.display_dir / f"{self.stem}_job_rendering.png"

    @property
    def annotation_txt_path(self) -> Path:
        return self.display_dir / f"{self.stem}_job_rendering.txt"

    @property
    def plain_language_txt_path(self) -> Path:
        return self.display_dir / f"{self.stem}_job_rendering_plain_language.txt"


@dataclass
class SlurmStepParams:
    kind: str
    path: Path
    params: dict[str, str]


@dataclass
class RunSlurmParams:
    steps: dict[str, SlurmStepParams]

    def get(self, kind: str) -> SlurmStepParams | None:
        return self.steps.get(kind)


def read_slurm_text(source: str | Path) -> str:
    if isinstance(source, Path):
        return source.read_text(encoding="utf-8", errors="replace")
    path = Path(source)
    if "\n" not in source and path.is_file():
        return path.read_text(encoding="utf-8", errors="replace")
    return str(source)


def slurm_header_lines(text: str) -> list[str]:
    lines: list[str] = []
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith(SLURM_HEADER_STOP_PREFIXES):
            break
        lines.append(line.rstrip("\n"))
    return lines


def parse_kv_header(text: str) -> dict[str, str]:
    params: dict[str, str] = {}
    for line in slurm_header_lines(text):
        stripped = line.strip()
        designs = N_DESIGNS_RE.match(stripped)
        if designs:
            params["n_designs"] = designs.group(1)
            params["n_batches"] = designs.group(2)
            if designs.group(3):
                params["batch_note"] = designs.group(3)
            continue
        match = KV_RE.match(stripped)
        if match:
            params[match.group(1)] = match.group(2)
    return params


def classify_slurm_out(source: str | Path) -> str | None:
    text = read_slurm_text(source)
    keys = set(parse_kv_header(text))
    first = next((line.strip() for line in text.splitlines() if line.strip()), "")
    if "start_stage" in keys or "n_designs" in keys or first.startswith("Checking Linux USalign"):
        return "pipeline"
    if "n_structures" in keys and "designed_chains" in keys:
        return "mpnn"
    if "n_inputs" in keys:
        return "rf3"
    if "n_batches" in keys and "yaml_src" in keys:
        return "rfd3"
    return None


def parse_pipeline_slurm(source: str | Path) -> dict[str, str]:
    return parse_kv_header(read_slurm_text(source))


def parse_rfd3_slurm(source: str | Path) -> dict[str, str]:
    return parse_kv_header(read_slurm_text(source))


def parse_mpnn_slurm(source: str | Path) -> dict[str, str]:
    return parse_kv_header(read_slurm_text(source))


def parse_rf3_slurm(source: str | Path) -> dict[str, str]:
    text = read_slurm_text(source)
    params = parse_kv_header(text)
    for line in text.splitlines()[:30]:
        stripped = line.strip()
        if stripped.startswith("rf3 ") or "num_steps=" in stripped:
            for key, value in RF3_FOLD_KV_RE.findall(stripped):
                params.setdefault(key, value)
            break
    return params


SLURM_PARSERS = {
    "pipeline": parse_pipeline_slurm,
    "rfd3": parse_rfd3_slurm,
    "mpnn": parse_mpnn_slurm,
    "rf3": parse_rf3_slurm,
}


def collect_run_slurm_params(run_dir: Path) -> RunSlurmParams:
    expected = {
        "pipeline": run_dir / "slurm.out",
        "rfd3": run_dir / "rfd3" / "slurm.out",
        "mpnn": run_dir / "mpnn" / "slurm.out",
        "rf3": run_dir / "rf3" / "slurm.out",
    }
    steps: dict[str, SlurmStepParams] = {}
    for expected_kind, path in expected.items():
        if not path.is_file():
            continue
        text = read_slurm_text(path)
        classified = classify_slurm_out(text) or expected_kind
        if classified != expected_kind:
            print(
                f"[SLURM] {path} expected {expected_kind}, classified as {classified}; "
                "using classified parser"
            )
        parser = SLURM_PARSERS[classified]
        steps[classified] = SlurmStepParams(
            kind=classified,
            path=path,
            params=parser(text),
        )
    return RunSlurmParams(steps)


def ckpt_basename(value: str | None) -> str | None:
    if not value:
        return None
    return Path(value).name


def parse_residue_key(key: str) -> ResidueRef:
    match = HOTSPOT_KEY_RE.match(str(key).strip())
    if not match:
        raise ValueError(f"Unrecognized residue key {key!r}; expected e.g. A30")
    return ResidueRef(chain=match.group(1), resi=int(match.group(2)))


def chains_from_contig(contig: str | None) -> list[str]:
    if not contig:
        return []
    chains: list[str] = []
    for part in str(contig).split(","):
        match = CONTIG_CHAIN_RE.match(part.strip())
        if match:
            chain = match.group(1)
            if chain not in chains:
                chains.append(chain)
    return chains


def receptor_chains(job: JobSpec) -> list[str]:
    chains: list[str] = []
    for key in job.hotspots:
        chain = parse_residue_key(key).chain
        if chain not in chains:
            chains.append(chain)
    for chain in chains_from_contig(job.contig):
        if chain not in chains:
            chains.append(chain)
    return chains or ["A"]


def residue_selection(keys: Iterable[str]) -> str:
    by_chain: dict[str, list[int]] = defaultdict(list)
    for key in keys:
        ref = parse_residue_key(key)
        by_chain[ref.chain].append(ref.resi)
    parts = []
    for chain, residues in by_chain.items():
        resi = "+".join(str(r) for r in sorted(set(residues)))
        parts.append(f"(chain {chain} and resi {resi})")
    return " or ".join(parts)


def hotspot_atom_names(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, list):
        return [str(item).strip() for item in value if str(item).strip()]
    return [part.strip() for part in str(value).split(",") if part.strip() and part.strip() != "[]"]


def ns1a_hotspot_color(job: JobSpec) -> str | None:
    match = NS1A_RUN_RE.fullmatch(job.run_name)
    if match is None:
        return None
    return match.group(1)


def hotspot_residue_keys(job: JobSpec) -> list[str]:
    return [key for key, value in job.hotspots.items() if hotspot_atom_names(value)]


def hotspot_selection(job: JobSpec) -> str:
    parts = []
    for key in hotspot_residue_keys(job):
        ref = parse_residue_key(key)
        names = "+".join(hotspot_atom_names(job.hotspots[key]))
        parts.append(f"(chain {ref.chain} and resi {ref.resi} and name {names})")
    return " or ".join(parts)


def is_bkbn(value: Any) -> bool:
    if value is None or isinstance(value, list):
        return False
    return str(value).strip().upper() == "BKBN"


def is_fully_flexible(value: Any) -> bool:
    if value is None:
        return True
    if isinstance(value, list):
        return len(value) == 0
    text = str(value).strip()
    return text == "" or text == "[]"


def flexible_keys(fixed_atoms: dict[str, Any] | None) -> tuple[list[str], list[str]]:
    if not fixed_atoms:
        return [], []
    bkbn: list[str] = []
    full: list[str] = []
    for key, value in fixed_atoms.items():
        if is_bkbn(value):
            bkbn.append(key)
        elif is_fully_flexible(value):
            full.append(key)
    return bkbn, full


def flexibility_mode(value: Any) -> str:
    if value is None:
        return "fully flexible ([])"
    if isinstance(value, list):
        if len(value) == 0:
            return "fully flexible ([])"
        atoms = ",".join(str(atom) for atom in value)
        return f"partial (fixed: {atoms})"
    text = str(value).strip()
    if text == "" or text == "[]":
        return "fully flexible ([])"
    if text.upper() == "BKBN":
        return "side chain (BKBN)"
    return f"partial (fixed: {text})"


def collapse_flexible_residues(fixed_atoms: dict[str, Any]) -> list[str]:
    items: list[tuple[ResidueRef, str]] = []
    for key, value in fixed_atoms.items():
        items.append((parse_residue_key(key), flexibility_mode(value)))
    items.sort(key=lambda item: (item[0].chain, item[0].resi, item[1]))

    if not items:
        return []

    lines: list[str] = []
    start, mode = items[0]
    prev = start

    def emit(first: ResidueRef, last: ResidueRef, label: str) -> str:
        if first.chain == last.chain and first.resi == last.resi:
            return f"  {first.label()}: {label}"
        return f"  {first.label()}–{last.label()}: {label}"

    for ref, label in items[1:]:
        consecutive = (
            ref.chain == prev.chain
            and ref.resi == prev.resi + 1
            and label == mode
        )
        if consecutive:
            prev = ref
            continue
        lines.append(emit(start, prev, mode))
        start = prev = ref
        mode = label
    lines.append(emit(start, prev, mode))
    return lines


def parse_contig_parts(contig: str | None) -> tuple[str | None, str | None]:
    if not contig:
        return None, None
    binder = None
    target = None
    for part in str(contig).split(","):
        part = part.strip()
        if part.startswith("/"):
            continue
        if CONTIG_BINDER_RE.match(part):
            low, high = part.split("-", 1)
            binder = f"{low}–{high}"
        elif CONTIG_TARGET_RE.match(part):
            target = part.replace("-", "–", 1)
    return binder, target


def hotspot_residue_labels(job: JobSpec) -> list[str]:
    refs = [parse_residue_key(key) for key in job.hotspots]
    refs.sort(key=lambda ref: (ref.chain, ref.resi))
    return [ref.label() for ref in refs]


def flag_is_on(value: str | None) -> bool | None:
    if value is None:
        return None
    text = str(value).strip().lower()
    if text in {"1", "true", "yes", "on"}:
        return True
    if text in {"0", "false", "no", "off"}:
        return False
    return None


def plain_language_text(job: JobSpec, slurm: RunSlurmParams | None = None) -> str:
    """Docs-based explanation of YAML + pipeline knobs (foundry/docs)."""
    paragraphs: list[str] = []
    binder, target = parse_contig_parts(job.contig)
    if job.is_non_loopy is False:
        style = "loopy binders"
    elif job.is_non_loopy is True:
        style = "non-loopy binders"
    else:
        style = "binders"

    recipe = f"This job designed {style}"
    if binder:
        recipe += f" of {binder} residues"
    if target:
        recipe += f" against the input target stretch {target}"
    recipe += "."
    hotspots = hotspot_residue_labels(job)
    if hotspots:
        recipe += (
            " The binder is aimed at isoform-fingerprint hotspot sidechains: "
            + ", ".join(hotspots)
            + "."
        )
    if job.fixed_atoms:
        flex = "; ".join(line.strip() for line in collapse_flexible_residues(job.fixed_atoms))
        recipe += (
            f" On the target, {flex}. BKBN means the backbone stays in place "
            "while the sidechain can move; [] would leave the whole residue unconstrained."
        )
    paragraphs.append(recipe)

    slurm = slurm or RunSlurmParams({})
    pipeline = slurm.get("pipeline")
    rfd3 = slurm.get("rfd3")
    mpnn = slurm.get("mpnn")
    rf3 = slurm.get("rf3")

    if pipeline or rfd3:
        bits = ["The design pipeline started"]
        start = _param(pipeline, "start_stage") if pipeline else None
        if start:
            bits.append(f"at {start}")
        else:
            bits.append("at RFD3")
        bits.append("and then ran SolubleMPNN and RF3 in the same GPU job.")
        n_designs = _param(pipeline, "n_designs") if pipeline else None
        n_batches = (_param(pipeline, "n_batches") if pipeline else None) or (
            _param(rfd3, "n_batches") if rfd3 else None
        )
        if n_designs and n_batches:
            approx = int(n_batches) * 8
            bits.append(
                f"It asked for {n_designs} backbones, which RFD3 ran as {n_batches} "
                f"batches of 8 (~{approx} CIFs)."
            )
        elif n_batches:
            bits.append(
                f"RFD3 ran {n_batches} batches (8 designs each)."
            )
        ckpt = ckpt_basename(
            (_param(rfd3, "ckpt") if rfd3 else None)
            or (_param(pipeline, "rfd3_ckpt") if pipeline else None)
        )
        if ckpt:
            bits.append(f"RFD3 weights: {ckpt}.")
        paragraphs.append(" ".join(bits))

    if mpnn:
        n_seqs = _param(mpnn, "n_seqs") or (_param(pipeline, "n_seqs") if pipeline else None)
        chains = _param(mpnn, "designed_chains") or "A"
        n_structures = _param(mpnn, "n_structures")
        iface = flag_is_on(_param(mpnn, "exclude_interface"))
        mpnn_bits = [
            "ProteinMPNN / SolubleMPNN then assigned sequences on "
            f"chain {chains} (the binder)."
        ]
        if n_seqs:
            mpnn_bits.append(f"It wrote {n_seqs} sequences per input structure.")
        if n_structures:
            mpnn_bits.append(f"The MPNN stage reported {n_structures} input structures.")
        if iface is False:
            mpnn_bits.append(
                "exclude_interface is off, so the whole designed chain was sequenced, "
                "including residues at the binder–target interface "
                "(the pipeline default is to skip those)."
            )
        elif iface is True:
            cutoff = _param(mpnn, "interface_cutoff") or "5.0"
            mpnn_bits.append(
                f"Interface residues (any atom within {cutoff} Å of another chain) "
                "were left as the RFD3 sequence; only non-interface residues were redesigned."
            )
        ckpt = ckpt_basename(_param(mpnn, "ckpt"))
        if ckpt:
            mpnn_bits.append(f"Sequence weights: {ckpt}.")
        paragraphs.append(" ".join(mpnn_bits))

    if rf3:
        n_inputs = _param(rf3, "n_inputs")
        steps = _param(rf3, "num_steps")
        batch = _param(rf3, "diffusion_batch_size")
        rf3_bits = ["RF3 refolded those sequences."]
        if n_inputs:
            rf3_bits.append(f"It folded {n_inputs} inputs.")
            n_seqs = _param(mpnn, "n_seqs") if mpnn else None
            if n_seqs and n_inputs.isdigit() and n_seqs.isdigit() and int(n_seqs):
                n_bb = int(n_inputs) // int(n_seqs)
                if n_bb:
                    rf3_bits.append(
                        f"That is about {n_bb} backbones × {n_seqs} sequences."
                    )
        speed = []
        if steps:
            speed.append(f"{steps} diffusion steps")
        if batch:
            speed.append(f"{batch} sample per input" if batch == "1" else f"{batch} samples per input")
        if speed:
            rf3_bits.append(
                "Sampling used "
                + " and ".join(speed)
                + " (faster than RF3’s defaults of 200 steps and 5 samples)."
            )
        ckpt = ckpt_basename(_param(rf3, "ckpt") or _param(rf3, "ckpt_path"))
        if ckpt:
            rf3_bits.append(f"Fold weights: {ckpt}.")
        paragraphs.append(" ".join(rf3_bits))

    return "\n\n".join(paragraphs).strip() + "\n"


def wrap_plain_language_for_png(text: str, width: int = 58) -> str:
    blocks = []
    for para in text.strip().split("\n\n"):
        blocks.append(textwrap.fill(para.strip(), width=width))
    return "\n\n".join(blocks)


def yaml_annotation_lines(job: JobSpec) -> list[str]:
    loopy = job.is_non_loopy
    if isinstance(loopy, bool):
        loopy_text = "true" if loopy else "false"
    elif loopy is None:
        loopy_text = "not set"
    else:
        loopy_text = str(loopy)

    lines = [
        f"job: {job.run_name} / {job.target}",
        f"is_non_loopy: {loopy_text}",
        f"contig: {job.contig if job.contig is not None else 'not set'}",
    ]
    if job.fixed_atoms:
        lines.append("flexible residues:")
        lines.extend(collapse_flexible_residues(job.fixed_atoms))
    hotspots = hotspot_residue_labels(job)
    if hotspots:
        lines.append("hotspots:")
        lines.extend(
            f"  {line}" for line in textwrap.wrap(", ".join(hotspots), width=56)
        )
    return lines


def _param(step: SlurmStepParams | None, key: str) -> str | None:
    if step is None:
        return None
    value = step.params.get(key)
    return value if value not in (None, "") else None


def compact_slurm_lines(slurm: RunSlurmParams) -> list[str]:
    lines: list[str] = []
    pipeline = slurm.get("pipeline")
    if pipeline:
        lines.append("pipeline")
        n_designs = _param(pipeline, "n_designs")
        n_batches = _param(pipeline, "n_batches")
        if n_designs and n_batches:
            lines.append(f"  n_designs: {n_designs} ({n_batches} batches)")
        elif n_designs:
            lines.append(f"  n_designs: {n_designs}")
        if _param(pipeline, "n_seqs"):
            lines.append(f"  n_seqs: {pipeline.params['n_seqs']}")
        if _param(pipeline, "exclude_interface"):
            lines.append(f"  exclude_interface: {pipeline.params['exclude_interface']}")
        if _param(pipeline, "start_stage"):
            lines.append(f"  start_stage: {pipeline.params['start_stage']}")

    rfd3 = slurm.get("rfd3")
    if rfd3:
        lines.append("rfd3")
        if _param(rfd3, "n_batches"):
            lines.append(f"  n_batches: {rfd3.params['n_batches']}")
        ckpt = ckpt_basename(_param(rfd3, "ckpt"))
        if ckpt:
            lines.append(f"  ckpt: {ckpt}")

    mpnn = slurm.get("mpnn")
    if mpnn:
        lines.append("mpnn")
        bits = []
        if _param(mpnn, "n_structures"):
            bits.append(f"n_structures: {mpnn.params['n_structures']}")
        if _param(mpnn, "designed_chains"):
            bits.append(f"chains: {mpnn.params['designed_chains']}")
        if _param(mpnn, "n_seqs"):
            bits.append(f"n_seqs: {mpnn.params['n_seqs']}")
        if bits:
            lines.append("  " + "  ".join(bits))
        extra = []
        if _param(mpnn, "exclude_interface"):
            extra.append(f"exclude_interface: {mpnn.params['exclude_interface']}")
        if _param(mpnn, "interface_cutoff"):
            extra.append(f"cutoff: {mpnn.params['interface_cutoff']}")
        if extra:
            lines.append("  " + "  ".join(extra))
        ckpt = ckpt_basename(_param(mpnn, "ckpt"))
        if ckpt:
            lines.append(f"  ckpt: {ckpt}")

    rf3 = slurm.get("rf3")
    if rf3:
        lines.append("rf3")
        bits = []
        if _param(rf3, "n_inputs"):
            bits.append(f"n_inputs: {rf3.params['n_inputs']}")
        if _param(rf3, "num_steps"):
            bits.append(f"steps: {rf3.params['num_steps']}")
        if bits:
            lines.append("  " + "  ".join(bits))
        ckpt = ckpt_basename(_param(rf3, "ckpt") or _param(rf3, "ckpt_path"))
        if ckpt:
            lines.append(f"  ckpt: {ckpt}")
    return lines


def full_slurm_lines(slurm: RunSlurmParams) -> list[str]:
    order = ("pipeline", "rfd3", "mpnn", "rf3")
    lines: list[str] = []
    for kind in order:
        step = slurm.get(kind)
        if step is None:
            continue
        lines.append(f"=== {kind} ({step.path}) ===")
        for key, value in step.params.items():
            lines.append(f"  {key}: {value}")
        lines.append("")
    if lines and lines[-1] == "":
        lines.pop()
    return lines


def annotation_text(
    job: JobSpec,
    slurm: RunSlurmParams | None = None,
    *,
    compact: bool = True,
) -> str:
    lines = yaml_annotation_lines(job)
    if slurm and slurm.steps:
        lines.append("")
        lines.extend(compact_slurm_lines(slurm) if compact else full_slurm_lines(slurm))
    return "\n".join(lines).rstrip() + "\n"


def view_json_path_for_pdb(pdb_path: Path) -> Path:
    return pdb_path.with_name(f"{pdb_path.stem}_pymol_view.json")


def discover_yaml_files(root: Path) -> list[Path]:
    root = root.resolve()
    if not root.is_dir():
        raise FileNotFoundError(f"YAML directory does not exist: {root}")

    top = sorted(
        path
        for path in root.iterdir()
        if path.is_file() and path.suffix.lower() in YAML_SUFFIXES
    )
    if top:
        return top

    found: list[Path] = []
    for path in sorted(root.rglob("*")):
        if not path.is_file() or path.suffix.lower() not in YAML_SUFFIXES:
            continue
        if path.parent.name == "rfd3":
            parent_run = path.parent.parent
            if any(
                child.suffix.lower() in YAML_SUFFIXES
                for child in parent_run.iterdir()
                if child.is_file()
            ):
                continue
        found.append(path)
    return found


def run_dir_for_yaml(yaml_path: Path) -> Path:
    if yaml_path.parent.name == "rfd3":
        return yaml_path.parent.parent
    return yaml_path.parent


def load_jobs(yaml_path: Path) -> list[JobSpec]:
    with yaml_path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"Expected a mapping of job names in {yaml_path}")

    run_dir = run_dir_for_yaml(yaml_path)
    jobs: list[JobSpec] = []
    for target, config in data.items():
        if not isinstance(config, dict):
            continue
        hotspots = config.get("select_hotspots") or {}
        if not isinstance(hotspots, dict):
            raise ValueError(f"{yaml_path} target {target}: select_hotspots must be a mapping")
        fixed = config.get("select_fixed_atoms")
        if fixed is not None and not isinstance(fixed, dict):
            raise ValueError(f"{yaml_path} target {target}: select_fixed_atoms must be a mapping")
        jobs.append(
            JobSpec(
                yaml_path=yaml_path,
                run_dir=run_dir,
                run_name=run_dir.name,
                target=str(target),
                contig=config.get("contig"),
                is_non_loopy=config.get("is_non_loopy"),
                hotspots=hotspots,
                fixed_atoms=fixed,
            )
        )
    if not jobs:
        raise ValueError(f"No job targets found in {yaml_path}")
    return jobs


def write_pml(job: JobSpec, pdb_path: Path) -> Path:
    chains = receptor_chains(job)
    chain_sel = " or ".join(f"(chain {chain})" for chain in chains)
    ns1a_color = ns1a_hotspot_color(job)
    receptor_color = "gray" if ns1a_color else "green"
    lines = [
        "reinitialize",
        f'load "{pdb_path.resolve()}", receptor',
        "hide everything, receptor",
        "show cartoon, receptor",
        f"color {receptor_color}, {chain_sel}",
        "bg_color white",
        "set ray_opaque_background, on",
        "set cartoon_fancy_helices, 1",
    ]
    bkbn_keys, full_flex_keys = flexible_keys(job.fixed_atoms)
    if bkbn_keys:
        lines.extend(
            [
                f"select flex_bkbn, {residue_selection(bkbn_keys)}",
                "color blue, flex_bkbn and name n+ca+c+o",
            ]
        )
    if full_flex_keys:
        lines.extend(
            [
                f"select flex_full, {residue_selection(full_flex_keys)}",
                "show sticks, flex_full",
                "color blue, flex_full",
            ]
        )
    residue_keys = hotspot_residue_keys(job)
    selection = hotspot_selection(job) if residue_keys else ""
    if selection and ns1a_color:
        lines.extend(
            [
                f"select hotspot_residues, {residue_selection(residue_keys)}",
                "show sticks, hotspot_residues",
                f"color {ns1a_color}, hotspot_residues",
            ]
        )
    elif selection:
        lines.extend(
            [
                f"select hotspot_residues, {residue_selection(residue_keys)}",
                "show sticks, hotspot_residues",
                f"select hotspots, {selection}",
                "color magenta, hotspots",
            ]
        )
        if bkbn_keys:
            lines.append("color blue, flex_bkbn and name n+ca+c+o")
    else:
        lines.append("select hotspots, none")
    lines.append("orient receptor")
    job.pml_path.parent.mkdir(parents=True, exist_ok=True)
    job.pml_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return job.pml_path


def save_pse_with_view(pml_path: Path, pse_path: Path, view: list[float]) -> Path:
    pse_path.parent.mkdir(parents=True, exist_ok=True)
    helper_script = f'''
from pymol import cmd
cmd.do({view_command(view)!r})
cmd.save({str(pse_path.resolve())!r})
cmd.quit()
'''
    with tempfile.NamedTemporaryFile(mode="w", suffix="_save_pse.py", delete=False) as tmp:
        tmp.write(helper_script)
        helper_path = tmp.name
    try:
        run_pymol_command(
            ["pymol", "-c", "-q", str(pml_path), "-r", helper_path],
            context=f"save_pse:{pse_path.name}",
        )
    finally:
        if os.path.exists(helper_path):
            os.remove(helper_path)
    if not pse_path.exists():
        raise RuntimeError(f"Failed to save PyMOL session: {pse_path}")
    return pse_path


def _annotation_font(size: int = 18):
    from PIL import ImageFont

    candidates = (
        "/System/Library/Fonts/Menlo.ttc",
        "/System/Library/Fonts/Supplemental/Courier New.ttf",
        "/Library/Fonts/Menlo.ttc",
        "Menlo.ttc",
        "DejaVuSansMono.ttf",
    )
    for path in candidates:
        try:
            return ImageFont.truetype(path, size)
        except OSError:
            continue
    return ImageFont.load_default()


# Match PyMOL named colors: green, magenta, blue.
LEGEND_RGB = {
    "receptor": (0, 255, 0),
    "hotspot": (255, 0, 255),
    "flexible": (0, 0, 255),
}


def legend_entries(job: JobSpec) -> list[tuple[str, tuple[int, int, int]]]:
    ns1a_color = ns1a_hotspot_color(job)
    if ns1a_color:
        entries = [("receptor", PYMOL_RGB["gray"])]
        if job.hotspots:
            entries.append(("hotspot", PYMOL_RGB[ns1a_color]))
    else:
        entries = [("receptor", LEGEND_RGB["receptor"])]
        if job.hotspots:
            entries.append(("hotspot", LEGEND_RGB["hotspot"]))
    bkbn_keys, full_flex_keys = flexible_keys(job.fixed_atoms)
    if bkbn_keys or full_flex_keys:
        entries.append(("flexible", LEGEND_RGB["flexible"]))
    return entries


def draw_pymol_legend(png_path: Path, entries: list[tuple[str, tuple[int, int, int]]]) -> Path:
    from PIL import Image, ImageDraw

    img = Image.open(png_path).convert("RGB")
    draw = ImageDraw.Draw(img)
    font = _annotation_font(45)
    swatch = 50
    pad = 34
    gap = 22
    row_h = 64
    labels = [name for name, _ in entries]
    text_w = 0
    for label in labels:
        bbox = draw.textbbox((0, 0), label, font=font)
        text_w = max(text_w, bbox[2] - bbox[0])
    box_w = pad * 2 + swatch + gap + text_w
    box_h = pad * 2 + row_h * len(entries)
    margin = 16
    x0 = img.size[0] - margin - box_w
    y0 = margin
    draw.rectangle([x0, y0, x0 + box_w, y0 + box_h], fill="white", outline="black", width=1)
    for i, (name, color) in enumerate(entries):
        y = y0 + pad + i * row_h
        sx = x0 + pad
        sy = y + (row_h - swatch) // 2
        draw.rectangle([sx, sy, sx + swatch, sy + swatch], fill=color, outline="black")
        draw.text((sx + swatch + gap, y + 2), name, fill="black", font=font)
    img.save(png_path)
    return png_path


def compose_annotated_png(pymol_png: Path, text: str, output_path: Path) -> Path:
    """Paste the PyMOL PNG at native size; grow the canvas if text is taller."""
    try:
        from PIL import Image, ImageDraw
    except ImportError:
        return _compose_annotated_png_matplotlib(pymol_png, text, output_path)

    structure = Image.open(pymol_png).convert("RGB")
    pymol_w, pymol_h = structure.size
    panel_w = 2016
    padding = 28
    spacing = 6
    font = _annotation_font(50)

    scratch = Image.new("RGB", (1, 1), "white")
    bbox = ImageDraw.Draw(scratch).multiline_textbbox(
        (padding, padding), text, font=font, spacing=spacing
    )
    text_h = bbox[3] + padding
    canvas_h = max(pymol_h, text_h)
    canvas = Image.new("RGB", (panel_w + pymol_w, canvas_h), "white")
    ImageDraw.Draw(canvas).multiline_text(
        (padding, padding), text, fill="black", font=font, spacing=spacing
    )
    canvas.paste(structure, (panel_w, 0))
    canvas.save(output_path)
    return output_path


def _compose_annotated_png_matplotlib(pymol_png: Path, text: str, output_path: Path) -> Path:
    """Fallback paste that still does not scale the PyMOL pixels."""
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg

    img = mpimg.imread(pymol_png)
    pymol_h, pymol_w = img.shape[0], img.shape[1]
    panel_w = 2016
    padding = 28
    fontsize = 31
    line_px = int(fontsize * 1.35) + 4
    text_h = padding * 2 + (text.count("\n") + 1) * line_px
    canvas_h = max(pymol_h, text_h)
    canvas_w = panel_w + pymol_w
    dpi = 100
    fig = plt.figure(figsize=(canvas_w / dpi, canvas_h / dpi), dpi=dpi)
    fig.patch.set_facecolor("white")
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set_xlim(0, canvas_w)
    ax.set_ylim(canvas_h, 0)
    ax.axis("off")
    ax.imshow(
        img,
        extent=(panel_w, panel_w + pymol_w, pymol_h, 0),
        origin="upper",
        interpolation="nearest",
        resample=False,
    )
    ax.text(
        padding,
        padding,
        text,
        va="top",
        ha="left",
        family="monospace",
        fontsize=fontsize,
        color="black",
        linespacing=1.25,
    )
    fig.savefig(output_path, dpi=dpi, facecolor="white")
    plt.close(fig)
    return output_path


def collect_jobs(yaml_dir: Path) -> list[JobSpec]:
    jobs: list[JobSpec] = []
    yaml_files = discover_yaml_files(yaml_dir)
    if not yaml_files:
        raise FileNotFoundError(f"No YAML files found under {yaml_dir}")
    for yaml_path in yaml_files:
        jobs.extend(load_jobs(yaml_path))
    return jobs


def render_jobs(
    jobs: Iterable[JobSpec],
    pdb_path: Path,
    *,
    recapture_view: bool = False,
    plain_language: bool = False,
    force: bool = False,
) -> None:
    pdb_path = pdb_path.resolve()
    if not pdb_path.is_file():
        raise FileNotFoundError(f"PDB not found: {pdb_path}")

    job_list = list(jobs)
    if not job_list:
        raise ValueError("No jobs to render")

    for job in job_list:
        write_pml(job, pdb_path)

    view_path = view_json_path_for_pdb(pdb_path)
    first = job_list[0]
    view = resolve_view(
        view_json=view_path,
        csp_mask_pml=first.pml_path,
        view_save_path=view_path,
        interactive=recapture_view or not view_path.exists(),
        label=first.stem,
    )

    for job in job_list:
        print(f"[JOB] {job.stem}")
        slurm = collect_run_slurm_params(job.run_dir)
        png_text = annotation_text(job, slurm, compact=True)
        txt_text = annotation_text(job, slurm, compact=False)
        job.annotation_txt_path.write_text(txt_text, encoding="utf-8")
        if plain_language:
            if job.plain_language_txt_path.is_file() and not force:
                plain = job.plain_language_txt_path.read_text(encoding="utf-8")
            else:
                plain = plain_language_text(job, slurm)
                job.plain_language_txt_path.write_text(plain, encoding="utf-8")
            png_text = (
                png_text.rstrip()
                + "\n\nplain language\n"
                + wrap_plain_language_for_png(plain)
                + "\n"
            )
        save_pse_with_view(job.pml_path, job.pse_path, view)
        render_pymol_panel_with_view(job.pml_path, view, job.pymol_png_path)
        draw_pymol_legend(job.pymol_png_path, legend_entries(job))
        compose_annotated_png(job.pymol_png_path, png_text, job.composite_png_path)
        print(f"  pse  {job.pse_path}")
        print(f"  png  {job.composite_png_path}")
        print(f"  txt  {job.annotation_txt_path}")
        if plain_language:
            print(f"  plain {job.plain_language_txt_path}")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create a unique PyMOL session and annotated rendering for each run YAML. "
            "Residues with a hotspot atom are shown as sticks; only the YAML-listed atoms are magenta; the receptor chain is green cartoon; "
            "BKBN flexible backbone is blue cartoon; fully flexible ([]) residues are blue sticks."
        )
    )
    parser.add_argument(
        "yaml_dir",
        type=Path,
        help="Directory containing a run YAML, or a parent directory of run folders",
    )
    parser.add_argument(
        "--pdb",
        type=Path,
        required=True,
        help="Receptor PDB loaded into every session",
    )
    parser.add_argument(
        "--recapture-view",
        action="store_true",
        help="Open PyMOL and press F5 to replace the saved camera JSON",
    )
    parser.add_argument(
        "--plain-language",
        action="store_true",
        help=(
            "Write a docs-based explanation of YAML and pipeline knobs next to the "
            "PNG and draw it under the left-panel annotation. Reuses an existing "
            "sidecar unless --force is set."
        ),
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Regenerate the plain-language sidecar even if it already exists",
    )
    parser.epilog = (
        "Camera JSON is written next to the PDB as <pdb_stem>_pymol_view.json. "
        "If it is missing, PyMOL opens so you can set the view and press F5. "
        "Annotated PNGs need matplotlib (or Pillow)."
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    jobs = collect_jobs(args.yaml_dir)
    render_jobs(
        jobs,
        args.pdb,
        recapture_view=args.recapture_view,
        plain_language=args.plain_language,
        force=args.force,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
