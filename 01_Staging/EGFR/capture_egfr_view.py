"""Capture the EGFR hotspot camera with F5, then bake it into the PML and PSE.

Uses capture_user_view_interactive from pymol_visualization_snippets:
PyMOL opens the hotspot PML, the helper binds F5 after that script loads,
F5 writes an 18-float view JSON and quits.
"""

from __future__ import annotations

import sys
from pathlib import Path

WEEK1 = Path(__file__).resolve().parent
SNIPPETS = WEEK1.parent / "pymol_visualization_snippets" / "snippets"
PML_PATH = WEEK1 / "EGFR_seg_medoid_hotspots.pml"
PSE_PATH = WEEK1 / "EGFR_seg_medoid_hotspots.pse"
VIEW_PATH = WEEK1 / "EGFR_seg_medoid_hotspots_view.json"

sys.path.insert(0, str(SNIPPETS))
from pymol_render import (  # noqa: E402
    capture_user_view_interactive,
    load_view,
    run_pymol_command,
    view_command,
)


def replace_set_view(pml_path: Path, command: str) -> None:
    text = pml_path.read_text(encoding="utf-8")
    start = text.find("set_view (")
    end = text.find("\n", start)
    if start < 0 or end < 0:
        raise RuntimeError(f"No set_view line in {pml_path}")
    pml_path.write_text(text[:start] + command + text[end:], encoding="utf-8")


def save_pse(pml_path: Path, pse_path: Path) -> None:
    run_pymol_command(
        [
            "pymol",
            "-c",
            "-q",
            str(pml_path),
            "-d",
            f"save {pse_path}",
            "-d",
            "quit",
        ],
        context=f"save:{pse_path.name}",
    )
    if not pse_path.exists():
        raise RuntimeError(f"PyMOL did not write {pse_path}")


def main() -> None:
    if VIEW_PATH.exists():
        VIEW_PATH.unlink()
    print("Rotate the EGFR session, then press F5. PyMOL will close and save that camera.")
    capture_user_view_interactive(
        PML_PATH,
        VIEW_PATH,
        label="EGFR medoid hotspots",
    )
    view = load_view(VIEW_PATH)
    command = view_command(view)
    replace_set_view(PML_PATH, command)
    save_pse(PML_PATH, PSE_PATH)
    print(command)
    print(f"Updated {PML_PATH.name} and {PSE_PATH.name}")


if __name__ == "__main__":
    main()
