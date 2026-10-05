#!/usr/bin/env python3
"""Render EGFR targeting slides from the stored hotspot camera.

Each scene loads EGFR_seg_medoid_hotspots.pse, shows a different set of
objects, then applies the 18-float view in EGFR_seg_medoid_hotspots_view.json.
The session file is not modified. Edit SCENES to change the next deck.
"""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from datetime import datetime
from pathlib import Path

EGFR_DIR = Path(__file__).resolve().parent
SNIPPETS = EGFR_DIR.parent / "pymol_visualization_snippets" / "snippets"
PSE_PATH = EGFR_DIR / "EGFR_seg_medoid_hotspots.pse"
VIEW_PATH = EGFR_DIR / "EGFR_seg_medoid_hotspots_view.json"
PPTX_DIR = EGFR_DIR.parent / "pptx"
PYMOL = Path("/opt/homebrew/bin/pymol")
CHALLENGE_URL = "https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/egfr"

RAY_WIDTH = 1600
RAY_HEIGHT = 1200

# Visibility only. The stored set_view is applied after these commands.
SCENES: list[dict[str, object]] = [
    {
        "slug": "01_human_cartoon",
        "title": "Human Domain III",
        "caption": (
            "Lightorange cartoon of the human EGFR segment medoid "
            "(file residues 1–192, EGFR 310–501). Domain III is the "
            "recommended epitope, also targeted by cetuximab and panitumumab.\n\n"
            "Shown: human cartoon."
        ),
        "commands": [
            "enable human",
            "show cartoon, human and polymer.protein",
        ],
    },
    {
        "slug": "02_mouse_cartoon",
        "title": "Mouse ortholog",
        "caption": (
            "Skyblue cartoon of the mouse EGFR segment medoid in the same camera. "
            "Mouse is already aligned onto human.\n\n"
            "Shown: mouse cartoon."
        ),
        "commands": [
            "enable mouse",
            "show cartoon, mouse and polymer.protein",
        ],
    },
    {
        "slug": "03_human_mouse_overlay",
        "title": "Human and mouse overlay",
        "caption": (
            "Both medoids in the stored camera. The shared Domain III fold is "
            "the cross-reactive surface.\n\n"
            "Shown: human (lightorange) and mouse (skyblue) cartoons."
        ),
        "commands": [
            "enable human",
            "enable mouse",
            "show cartoon, (human or mouse) and polymer.protein",
        ],
    },
    {
        "slug": "04_mismatches",
        "title": "Sequence mismatches",
        "caption": (
            "Green sticks mark positions that differ between human and mouse "
            "in this segment. File residue N is EGFR residue N+309.\n\n"
            "Shown: both cartoons, mismatch sticks."
        ),
        "commands": [
            "enable human",
            "enable mouse",
            "show cartoon, (human or mouse) and polymer.protein",
            "show sticks, mismatches",
        ],
    },
    {
        "slug": "05_histidines",
        "title": "Histidines",
        "caption": (
            "Purple sticks are histidines on both species. These are the main "
            "candidates for binding at pH 6.5 and not at pH 7.4.\n\n"
            "Shown: both cartoons, His sticks."
        ),
        "commands": [
            "enable human",
            "enable mouse",
            "show cartoon, (human or mouse) and polymer.protein",
            "show sticks, Human_HIS or Mouse_HIS",
        ],
    },
    {
        "slug": "06_arg_lys_his",
        "title": "Arg, Lys, and His",
        "caption": (
            "Red sticks are arginine and lysine. Purple sticks are histidine. "
            "Basic residues sit next to the pH-sensitive set.\n\n"
            "Shown: both cartoons, Arg/Lys and His sticks."
        ),
        "commands": [
            "enable human",
            "enable mouse",
            "show cartoon, (human or mouse) and polymer.protein",
            "show sticks, Human_KRL or Mouse_KRL or Human_HIS or Mouse_HIS",
        ],
    },
    {
        "slug": "07_glycans",
        "title": "N-linked glycans",
        "caption": (
            "Yellow sticks are NAG glycans on human and mouse. Glycans can "
            "cover parts of the Domain III surface.\n\n"
            "Shown: both cartoons, NAG sticks."
        ),
        "commands": [
            "enable human",
            "enable mouse",
            "show cartoon, (human or mouse) and polymer.protein",
            "show sticks, Human_NAG or Mouse_NAG",
        ],
    },
    {
        "slug": "08_6aru_egfr",
        "title": "6ARU EGFR chain A",
        "caption": (
            "PDB 6ARU chain A, rainbow from N to C, aligned onto human "
            "residues 310–501. Cetuximab is hidden.\n\n"
            "Shown: 6ARU EGFR cartoon."
        ),
        "commands": [
            "enable 6aru",
            "show cartoon, 6aru_EGFR",
        ],
    },
    {
        "slug": "09_6aru_cetuximab",
        "title": "6ARU EGFR plus cetuximab",
        "caption": (
            "6ARU EGFR chain A with the cetuximab Fab in the same camera. "
            "Light chain is B. Heavy chain is C.\n\n"
            "Shown: 6ARU EGFR and cetuximab cartoons."
        ),
        "commands": [
            "enable 6aru",
            "show cartoon, 6aru and polymer.protein",
        ],
    },
    {
        "slug": "10_human_cetuximab",
        "title": "Human medoid plus cetuximab",
        "caption": (
            "Human Domain III medoid with the cetuximab chains. 6ARU EGFR is "
            "hidden so the clinical epitope sits on the design target.\n\n"
            "Shown: human cartoon, cetuximab light and heavy chains."
        ),
        "commands": [
            "enable human",
            "enable 6aru",
            "show cartoon, human and polymer.protein",
            "show cartoon, 6aru_cetux_light or 6aru_cetux_heavy",
        ],
    },
]


def ensure_pptx() -> None:
    vendor = EGFR_DIR / ".vendor"
    if vendor.is_dir():
        sys.path.insert(0, str(vendor))
    try:
        import pptx  # noqa: F401
    except ImportError:
        vendor.mkdir(exist_ok=True)
        subprocess.check_call(
            [sys.executable, "-m", "pip", "install", "--target", str(vendor), "python-pptx"]
        )
        if str(vendor) not in sys.path:
            sys.path.insert(0, str(vendor))
        import pptx  # noqa: F401


def pymol_executable() -> str:
    if PYMOL.is_file():
        return str(PYMOL)
    return "pymol"


def render_pngs(png_dir: Path, view_command_text: str) -> None:
    png_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "pse": str(PSE_PATH),
        "view_command": view_command_text,
        "width": RAY_WIDTH,
        "height": RAY_HEIGHT,
        "scenes": [
            {
                "slug": scene["slug"],
                "png": str(png_dir / f"{scene['slug']}.png"),
                "commands": scene["commands"],
            }
            for scene in SCENES
        ],
    }
    script = f"""
import json
import os
from pymol import cmd

JOB = json.loads({json.dumps(json.dumps(payload))})

cmd.load(JOB["pse"])
cmd.set("auto_zoom", 0)
cmd.set("orthoscopic", 1)
cmd.set("ray_opaque_background", 1)
cmd.bg_color("white")

for scene in JOB["scenes"]:
    print("SCENE", scene["slug"], flush=True)
    cmd.disable("all")
    cmd.hide("everything", "all")
    cmd.set("ray_opaque_background", 1)
    cmd.set("orthoscopic", 1)
    cmd.bg_color("white")
    cmd.viewport(JOB["width"], JOB["height"])
    for command in scene["commands"]:
        cmd.do(command)
    cmd.do(JOB["view_command"])
    cmd.ray(JOB["width"], JOB["height"])
    cmd.png(scene["png"], dpi=300)
    if not os.path.exists(scene["png"]) or os.path.getsize(scene["png"]) < 1000:
        raise RuntimeError("PyMOL did not write " + scene["png"])
    print("WROTE", scene["png"], flush=True)

cmd.quit()
"""
    with tempfile.NamedTemporaryFile(
        mode="w", suffix="_egfr_targeting.py", delete=False
    ) as handle:
        handle.write(script)
        script_path = handle.name
    try:
        result = subprocess.run(
            [pymol_executable(), "-c", "-q", "-r", script_path],
            check=False,
        )
    finally:
        Path(script_path).unlink(missing_ok=True)
    if result.returncode != 0:
        raise RuntimeError(f"PyMOL failed while rendering EGFR scenes (exit {result.returncode})")
    missing = [png_dir / f"{scene['slug']}.png" for scene in SCENES]
    missing = [path for path in missing if not path.exists()]
    if missing:
        raise RuntimeError(f"Missing renders: {', '.join(path.name for path in missing)}")


def _add_title_slide(prs) -> None:
    from pptx.util import Inches, Pt

    slide = prs.slides.add_slide(prs.slide_layouts[6])
    title = slide.shapes.add_textbox(Inches(0.7), Inches(0.55), Inches(12.0), Inches(1.0))
    title_run = title.text_frame.paragraphs[0].add_run()
    title_run.text = "EGFR targeting"
    title_run.font.size = Pt(36)
    title_run.font.bold = True

    body = slide.shapes.add_textbox(Inches(0.7), Inches(1.8), Inches(11.8), Inches(4.8))
    frame = body.text_frame
    frame.word_wrap = True
    paragraphs = [
        (
            "This week’s goal is a pH-sensitive binder for Domain III of human and mouse EGFR, "
            "for Challenge 01: EGFR on ProteinBase (Adaptyv and Anthropic)."
        ),
        "Every structure slide uses the camera stored in EGFR_seg_medoid_hotspots.pse. "
        "File residue N on the medoids is EGFR residue N+309 (segment 310–501).",
        "Human cartoon is lightorange. Mouse cartoon is skyblue. Mismatches are green sticks. "
        "Arg and Lys are red. His is purple. NAG is yellow. 6ARU chain A is rainbow. "
        "Cetuximab light chain is B and heavy chain is C.",
    ]
    frame.paragraphs[0].text = paragraphs[0]
    for line in paragraphs[1:]:
        paragraph = frame.add_paragraph()
        paragraph.text = line
    for paragraph in frame.paragraphs:
        paragraph.font.size = Pt(18)
        paragraph.space_after = Pt(12)

    link_box = slide.shapes.add_textbox(Inches(0.7), Inches(6.5), Inches(11.8), Inches(0.5))
    link_run = link_box.text_frame.paragraphs[0].add_run()
    link_run.text = CHALLENGE_URL
    link_run.font.size = Pt(14)
    link_run.hyperlink.address = CHALLENGE_URL


def _add_scene_slide(prs, scene: dict[str, object], png_path: Path) -> None:
    from pptx.util import Inches, Pt

    slide = prs.slides.add_slide(prs.slide_layouts[6])
    picture_width = 8.2
    picture_height = picture_width * RAY_HEIGHT / RAY_WIDTH
    slide.shapes.add_picture(
        str(png_path),
        Inches(0.25),
        Inches((7.5 - picture_height) / 2),
        Inches(picture_width),
        Inches(picture_height),
    )

    title = slide.shapes.add_textbox(Inches(8.65), Inches(0.45), Inches(4.4), Inches(1.1))
    title_frame = title.text_frame
    title_frame.word_wrap = True
    title_run = title_frame.paragraphs[0].add_run()
    title_run.text = str(scene["title"])
    title_run.font.size = Pt(22)
    title_run.font.bold = True

    caption = slide.shapes.add_textbox(Inches(8.65), Inches(1.7), Inches(4.4), Inches(5.2))
    caption_frame = caption.text_frame
    caption_frame.word_wrap = True
    blocks = str(scene["caption"]).split("\n\n")
    caption_frame.paragraphs[0].text = blocks[0]
    for block in blocks[1:]:
        paragraph = caption_frame.add_paragraph()
        paragraph.text = block
    for paragraph in caption_frame.paragraphs:
        paragraph.font.size = Pt(15)
        paragraph.space_after = Pt(10)


def build_pptx(png_dir: Path, pptx_path: Path) -> None:
    from pptx import Presentation
    from pptx.util import Inches

    prs = Presentation()
    prs.slide_width = Inches(13.333)
    prs.slide_height = Inches(7.5)
    _add_title_slide(prs)
    for scene in SCENES:
        _add_scene_slide(prs, scene, png_dir / f"{scene['slug']}.png")
    pptx_path.parent.mkdir(parents=True, exist_ok=True)
    prs.save(pptx_path)


def main() -> None:
    sys.path.insert(0, str(SNIPPETS))
    from pymol_render import load_view, view_command

    if not PSE_PATH.exists():
        raise FileNotFoundError(PSE_PATH)
    if not VIEW_PATH.exists():
        raise FileNotFoundError(VIEW_PATH)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    png_dir = PPTX_DIR / "renders" / timestamp
    pptx_path = PPTX_DIR / f"EGFR_targeting_{timestamp}.pptx"
    render_pngs(png_dir, view_command(load_view(VIEW_PATH)))
    ensure_pptx()
    build_pptx(png_dir, pptx_path)
    print(f"Wrote {pptx_path}")


if __name__ == "__main__":
    main()
