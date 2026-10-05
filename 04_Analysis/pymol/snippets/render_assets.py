"""
Orchestrate CSV -> .pml -> case_study_assets PNG rendering.

Source: scripts/case_study.py generate_case_study_figure render loop (336-368),
without matplotlib composition.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Literal, Optional

from .csv_adapter import load_master_alignment
from .pymol_render import render_pymol_panel_with_view, resolve_view
from .pymol_script_writers import write_all_case_study_scripts

Variant = Literal["v1", "v2", "both"]


@dataclass
class RenderAssetsConfig:
    """Configuration for one target's case-study PyMOL asset render."""

    master_csv: Path
    structure_pdb: Path
    occlusion_pdb: Path
    output_dir: Path
    pml_dir: Optional[Path] = None
    view_json: Optional[Path] = None
    view_save_path: Optional[Path] = None
    receptor_chain: Optional[str] = None
    ligand_chain: Optional[str] = None
    interactive_view: bool = False
    variant: Variant = "both"
    ray_width: int = 1600
    ray_height: int = 1200
    ray_dpi: int = 300


def _asset_names(variant: Variant) -> Dict[str, str]:
    """Map panel key -> output PNG basename for v1 / v2 naming."""
    if variant == "v1":
        prefix = "case_study"
    elif variant == "v2":
        prefix = "case_study2"
    else:
        raise ValueError("Use variant='both' via render_case_study_assets, not _asset_names")

    return {
        "csp_mask": f"{prefix}_color_csp_mask.png",
        "occlusion": f"{prefix}_color_occlusion.png",
        "classification": f"{prefix}_csp_classification_original.png",
    }


def render_case_study_assets(config: RenderAssetsConfig) -> Dict[str, Path]:
    """
    Full pipeline for one variant (v1 or v2 when config.variant is v1/v2).

    Steps:
      1. Load master_alignment.csv
      2. Write three .pml scripts
      3. Resolve camera view (JSON reuse or interactive F5)
      4. Headless ray-trace three PNGs into output_dir
    """
    if config.variant not in ("v1", "v2"):
        raise ValueError("config.variant must be 'v1' or 'v2' for a single render call")

    data = load_master_alignment(
        config.master_csv,
        receptor_chain=config.receptor_chain,
        ligand_chain=config.ligand_chain,
    )

    pml_dir = config.pml_dir or config.output_dir.parent
    scripts = write_all_case_study_scripts(
        data,
        structure_pdb=config.structure_pdb,
        occlusion_pdb=config.occlusion_pdb,
        output_dir=pml_dir,
    )

    view = resolve_view(
        view_json=config.view_json,
        csp_mask_pml=scripts["csp_mask"],
        view_save_path=config.view_save_path,
        interactive=config.interactive_view,
        label=data.holo_pdb or config.master_csv.stem,
    )

    config.output_dir.mkdir(parents=True, exist_ok=True)
    names = _asset_names(config.variant)
    rendered: Dict[str, Path] = {}

    panel_scripts = {
        "csp_mask": scripts["csp_mask"],
        "occlusion": scripts["occlusion"],
        "classification": scripts["classification"],
    }

    for key, pml_path in panel_scripts.items():
        png_path = config.output_dir / names[key]
        render_pymol_panel_with_view(
            pml_path,
            view,
            png_path,
            width=config.ray_width,
            height=config.ray_height,
            dpi=config.ray_dpi,
        )
        rendered[key] = png_path
        print(f"[RENDER] Saved: {png_path}")

    return rendered


def render_all_variants(config: RenderAssetsConfig) -> Dict[str, Path]:
    """Render v1 and/or v2 asset sets according to config.variant."""
    all_paths: Dict[str, Path] = {}
    variants: List[Literal["v1", "v2"]] = (
        ["v1", "v2"] if config.variant == "both" else [config.variant]  # type: ignore[list-item]
    )

    for variant in variants:
        sub_config = RenderAssetsConfig(
            master_csv=config.master_csv,
            structure_pdb=config.structure_pdb,
            occlusion_pdb=config.occlusion_pdb,
            output_dir=config.output_dir,
            pml_dir=config.pml_dir,
            view_json=config.view_json,
            view_save_path=config.view_save_path,
            receptor_chain=config.receptor_chain,
            ligand_chain=config.ligand_chain,
            interactive_view=config.interactive_view and variant == variants[0],
            variant=variant,
            ray_width=config.ray_width,
            ray_height=config.ray_height,
            ray_dpi=config.ray_dpi,
        )
        paths = render_case_study_assets(sub_config)
        for key, path in paths.items():
            all_paths[f"{variant}_{key}"] = path

    return all_paths
