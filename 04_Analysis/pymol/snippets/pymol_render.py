"""
Headless PyMOL rendering and optional interactive camera capture.

Source: scripts/case_study.py (lines 83-184).
"""

from __future__ import annotations

import json
import math
import os
import subprocess
import tempfile
from pathlib import Path
from typing import List, Optional, Sequence

_REPO = Path(__file__).resolve().parents[3]
ET_REFERENCE = (
    _REPO
    / "design"
    / "rfd3"
    / "foundry"
    / "runs"
    / "Brd4ET_hot_7"
    / "Brd4ET_closed_state_fix.pdb"
)
ET_VIEW = (
    _REPO
    / "design"
    / "rfd3"
    / "foundry"
    / "runs"
    / "Brd4ET_hot_7"
    / "Brd4ET_closed_state_fix_pymol_view.json"
)


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


def _view_coordinates(point: Sequence[float], view: Sequence[float]) -> tuple[float, float, float]:
    """Map a model-space point into the saved camera frame."""
    px = point[0] - view[12]
    py = point[1] - view[13]
    pz = point[2] - view[14]
    vx = view[0] * px + view[1] * py + view[2] * pz
    vy = view[3] * px + view[4] * py + view[5] * pz
    vz = view[6] * px + view[7] * py + view[8] * pz
    return vx, vy, vz


def _pullback_for_point(
    point: Sequence[float],
    view: Sequence[float],
    *,
    width: float,
    height: float,
) -> tuple[float, float]:
    """Clip scale and lateral scale that keep one atom inside the saved frustum."""
    vx, vy, vz = _view_coordinates(point, view)
    camera_z = view[11]
    front = view[15]
    back = view[16]
    near = front + camera_z
    far = back + camera_z
    if vz >= 0:
        clip = vz / far if far > 0 else 1.0
    else:
        clip = vz / near if near < 0 else 1.0
    fov = abs(view[17]) if view[17] < 0 else abs(view[17] or 20.0)
    tangent = math.tan(math.radians(fov / 2.0))
    depth_scale = -camera_z * tangent
    aspect = width / height if height else 1.0
    if depth_scale <= 1e-6:
        return max(clip, 1.0), 1.0
    lateral_y = (abs(vy - view[10]) - vz * tangent) / depth_scale
    lateral_x = (abs(vx - view[9]) / aspect - vz * tangent) / depth_scale
    return max(clip, 1.0), max(lateral_x, lateral_y, 1.0)


def scale_view_to_fit(
    view: Sequence[float],
    points: Sequence[Sequence[float]],
    *,
    width: float,
    height: float,
    reference_points: Sequence[Sequence[float]] = (),
    margin: float = 1.04,
) -> List[float]:
    """Pull the camera back until every point fits, without moving ET on screen.

    Only view elements 11, 15, and 16 change. Rotation and the look-at point stay
    as captured. The scale is at least 1, so the saved ET face is never cropped
    tighter than the reference camera. Lateral slack already present on the
    closed-state ET is not counted again, so that face keeps its captured size
    until a binder actually sticks out past it.
    """
    fitted = [float(value) for value in view]
    if len(fitted) != 18:
        raise ValueError("PyMOL view must contain 18 floats.")

    def limits(samples: Sequence[Sequence[float]]) -> tuple[float, float]:
        clip = 1.0
        lateral = 1.0
        for point in samples:
            point_clip, point_lateral = _pullback_for_point(
                point, fitted, width=width, height=height
            )
            clip = max(clip, point_clip)
            lateral = max(lateral, point_lateral)
        return clip, lateral

    clip, lateral = limits(points)
    _, reference_lateral = limits(reference_points)
    if reference_lateral > 1.0:
        lateral = max(1.0, lateral / reference_lateral)
    factor = max(clip, lateral, 1.0) * margin
    fitted[11] *= factor
    fitted[15] *= factor
    fitted[16] *= factor
    return fitted


def anchor_et_view(mobile: str, et_chain: str = "A", reference: str = "receptor") -> List[float]:
    """Superpose one complex onto the closed-state ET and return that camera.

    Call this from inside PyMOL. The closed-state receptor is loaded once,
    hidden, and left in place so later complexes share the same frame.
    The returned view is the saved camera with no ``orient`` or ``zoom``.
    """
    from pymol import cmd

    if not ET_REFERENCE.is_file() or not ET_VIEW.is_file():
        raise FileNotFoundError("Closed-state ET structure or camera view is missing")
    if reference not in cmd.get_object_list():
        cmd.load(str(ET_REFERENCE.resolve()), reference)
    cmd.disable(reference)
    cmd.align(
        f"{mobile} and chain {et_chain} and name CA",
        f"{reference} and chain A and name CA",
    )
    view = load_view(ET_VIEW)
    cmd.set_view(view)
    return view


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
