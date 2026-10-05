"""
Headless PyMOL rendering and optional interactive camera capture.

Source: scripts/case_study.py (lines 83-184).
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path
from typing import List, Optional


def run_pymol_command(cmd: List[str], *, context: str = "") -> None:
    """Run PyMOL subprocess; raise on non-zero exit."""
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.stdout:
        print(f"[PYMOL] {context} stdout:\n{result.stdout}")
    if result.stderr:
        print(f"[PYMOL] {context} stderr:\n{result.stderr}")
    if result.returncode != 0:
        raise RuntimeError(
            f"PyMOL failed for {context or cmd} (exit code {result.returncode})"
        )


def load_view(view_json_path: str | Path) -> List[float]:
    """Load 18-float PyMOL camera view from JSON."""
    with Path(view_json_path).open("r", encoding="utf-8") as handle:
        view = json.load(handle)
    if not isinstance(view, list) or len(view) != 18:
        raise ValueError(
            f"Invalid PyMOL view in {view_json_path}; expected list of 18 floats."
        )
    return [float(x) for x in view]


def view_command(view: List[float]) -> str:
    """Format PyMOL ``set_view`` command."""
    return "set_view (" + ", ".join(f"{value:.9f}" for value in view) + ")"


def capture_user_view_interactive(
    pml_path: str | Path,
    view_output_path: str | Path,
    *,
    label: str = "structure",
) -> None:
    """
    Open PyMOL GUI with ``pml_path``; user orients model and presses F5 to save view.

    Requires a display (not headless). Saves 18-float view JSON to view_output_path.
    """
    view_path = Path(view_output_path)
    view_path.parent.mkdir(parents=True, exist_ok=True)

    helper_script = f'''
from pymol import cmd
import json

VIEW_PATH = r"""{view_path.resolve()}"""

def _save_case_study_view():
    view = list(cmd.get_view())
    with open(VIEW_PATH, "w", encoding="utf-8") as handle:
        json.dump(view, handle)
    print(f"[CASE_STUDY] Saved view to {{VIEW_PATH}}")
    cmd.quit()

cmd.set_key("F5", _save_case_study_view)
print("[CASE_STUDY] ------------------------------------------------------------")
print("[CASE_STUDY] Case-study view capture for {label}")
print("[CASE_STUDY] 1) Set your desired model orientation.")
print("[CASE_STUDY] 2) Press F5 to save the view and close PyMOL.")
print("[CASE_STUDY] ------------------------------------------------------------")
'''

    with tempfile.NamedTemporaryFile(
        mode="w", suffix="_case_study_capture.py", delete=False
    ) as tmp:
        tmp.write(helper_script)
        helper_path = tmp.name

    try:
        run_pymol_command(
            ["pymol", str(pml_path), "-r", helper_path],
            context=f"capture_user_view_interactive:{label}",
        )
    finally:
        if os.path.exists(helper_path):
            os.remove(helper_path)

    if not view_path.exists():
        raise RuntimeError(
            "No view was captured. Re-run and press F5 in the PyMOL window."
        )


def render_pymol_panel_with_view(
    pml_path: str | Path,
    view: List[float],
    output_png_path: str | Path,
    *,
    width: int = 1600,
    height: int = 1200,
    dpi: int = 300,
) -> Path:
    """
    Ray-trace a PyMOL script headlessly using a fixed camera view.

    Equivalent to case_study.render_pymol_panel_with_view.
    """
    out = Path(output_png_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    run_pymol_command(
        [
            "pymol",
            "-c",
            "-q",
            str(pml_path),
            "-d",
            view_command(view),
            "-d",
            f"ray {width}, {height}",
            "-d",
            f"png {out.resolve()}, dpi={dpi}",
            "-d",
            "quit",
        ],
        context=f"render:{out.name}",
    )

    if not out.exists():
        raise RuntimeError(f"Failed to render PyMOL panel: {out}")
    return out


def resolve_view(
    *,
    view_json: Optional[str | Path] = None,
    csp_mask_pml: Optional[str | Path] = None,
    view_save_path: Optional[str | Path] = None,
    interactive: bool = False,
    label: str = "structure",
) -> List[float]:
    """
    Load cached view JSON, or capture interactively when missing / requested.

    Args:
        view_json: Existing view file to load
        csp_mask_pml: .pml opened for interactive capture (typically color_csp_mask.pml)
        view_save_path: Where to write view after interactive capture
        interactive: Force interactive capture even if view_json exists
        label: Label printed during interactive capture
    """
    if view_json and Path(view_json).exists() and not interactive:
        print(f"[PYMOL] Reusing saved view: {view_json}")
        return load_view(view_json)

    if csp_mask_pml is None or view_save_path is None:
        raise FileNotFoundError(
            "No view JSON found and interactive capture requires csp_mask_pml + view_save_path."
        )

    capture_user_view_interactive(csp_mask_pml, view_save_path, label=label)
    return load_view(view_save_path)
