#!/usr/bin/env python3
"""Scatter plot of DeltaForge Kd on Brd4ET vs Brd3ET for matched binder sequences."""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

RFD3_OUT = Path(__file__).resolve().parents[3] / "analysis" / "outputs" / "rfd3"
SCORES = RFD3_OUT
BRD4ET_FILENAME = "deltaforge_scores.json"
BRD3ET_FILENAME = "Brd3ET_deltaforge_scores.json"
OUTPUT_PATH = RFD3_OUT / "brd4et_vs_brd3et_kd.png"

FOLDER_COLORS = {
    "Brd4ET_design": "#2563eb",
    "Brd4ET_design_2": "#7c3aed",
    "Brd4ET_hot_4": "#16a34a",
    "Brd4ET_hot_4_flexible_30_50": "#ea580c",
    "Brd4ET_hot_4_loopy": "#db2777",
    "Brd4ET_hot_4_loopy_flexible_30_50": "#0891b2",
}
FALLBACK_COLORS = (
    "#64748b",
    "#ca8a04",
    "#0f766e",
    "#be123c",
    "#4338ca",
)


def load_rows(path: Path) -> list[dict[str, Any]]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, list):
        raise ValueError(f"Expected a JSON array in {path}")
    return [row for row in data if isinstance(row, dict)]


def _as_positive_float(value: Any) -> float | None:
    if value is None:
        return None
    try:
        kd = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(kd) or kd <= 0:
        return None
    return kd


def kd_from_result(row: dict[str, Any]) -> float | None:
    result = row.get("result")
    if not isinstance(result, dict):
        return None
    kd = _as_positive_float(result.get("kd_nm"))
    if kd is not None:
        return kd
    deltaforge = result.get("deltaforge")
    if isinstance(deltaforge, dict):
        return _as_positive_float(deltaforge.get("kd_nm"))
    return None


def index_by_sequence(rows: list[dict[str, Any]]) -> dict[str, float | None]:
    out: dict[str, float | None] = {}
    for row in rows:
        sequence = row.get("binder_sequence")
        if not isinstance(sequence, str) or not sequence:
            continue
        out[sequence] = kd_from_result(row)
    return out


def discover_folders(root: Path) -> list[Path]:
    folders = {path.parent for path in root.rglob(BRD4ET_FILENAME)}
    folders.update(path.parent for path in root.rglob(BRD3ET_FILENAME))
    return sorted(folders, key=lambda p: p.relative_to(root).as_posix())


def color_for(folder_name: str, used_fallback: list[int]) -> str:
    if folder_name in FOLDER_COLORS:
        return FOLDER_COLORS[folder_name]
    idx = used_fallback[0] % len(FALLBACK_COLORS)
    used_fallback[0] += 1
    return FALLBACK_COLORS[idx]


def collect_pairs(root: Path) -> tuple[dict[str, list[tuple[float, float]]], list[str]]:
    pairs_by_folder: dict[str, list[tuple[float, float]]] = defaultdict(list)
    summary: list[str] = []

    folders = discover_folders(root)
    summary.append(f"Folders scanned: {len(folders)}")

    plotted = 0
    skipped = 0
    skip_reasons: dict[str, int] = defaultdict(int)

    for folder in folders:
        brd4_path = folder / BRD4ET_FILENAME
        brd3_path = folder / BRD3ET_FILENAME
        label = folder.relative_to(root).as_posix()

        if not brd4_path.is_file() and not brd3_path.is_file():
            continue
        if not brd4_path.is_file():
            skip_reasons["missing Brd4ET file"] += 1
            skipped += 1
            summary.append(f"  {label}: skipped (no {BRD4ET_FILENAME})")
            continue
        if not brd3_path.is_file():
            n_brd4 = len(index_by_sequence(load_rows(brd4_path)))
            skip_reasons["missing Brd3ET file"] += n_brd4
            skipped += n_brd4
            summary.append(f"  {label}: skipped (no {BRD3ET_FILENAME}; {n_brd4} Brd4ET sequences)")
            continue

        brd4 = index_by_sequence(load_rows(brd4_path))
        brd3 = index_by_sequence(load_rows(brd3_path))
        sequences = sorted(set(brd4) | set(brd3))
        folder_plotted = 0
        folder_skipped = 0

        for sequence in sequences:
            if sequence not in brd4:
                skip_reasons["sequence only in Brd3ET"] += 1
                folder_skipped += 1
                continue
            if sequence not in brd3:
                skip_reasons["sequence only in Brd4ET"] += 1
                folder_skipped += 1
                continue
            x = brd4[sequence]
            y = brd3[sequence]
            if x is None and y is None:
                skip_reasons["missing kd_nm on both"] += 1
                folder_skipped += 1
                continue
            if x is None:
                skip_reasons["missing Brd4ET kd_nm"] += 1
                folder_skipped += 1
                continue
            if y is None:
                skip_reasons["missing Brd3ET kd_nm"] += 1
                folder_skipped += 1
                continue
            pairs_by_folder[label].append((x, y))
            folder_plotted += 1

        plotted += folder_plotted
        skipped += folder_skipped
        summary.append(
            f"  {label}: plotted {folder_plotted}, skipped {folder_skipped}"
        )

    summary.append(f"Pairs plotted: {plotted}")
    summary.append(f"Pairs skipped: {skipped}")
    if skip_reasons:
        summary.append("Skip reasons:")
        for reason, count in sorted(skip_reasons.items()):
            summary.append(f"  {reason}: {count}")
    return pairs_by_folder, summary


def main() -> None:
    pairs_by_folder, summary = collect_pairs(SCORES)
    for line in summary:
        print(line)

    if not pairs_by_folder:
        raise SystemExit("No complete Brd4ET/Brd3ET kd_nm pairs found.")

    fig, ax = plt.subplots(figsize=(8, 8))
    xs: list[float] = []
    ys: list[float] = []
    used_fallback = [0]

    for folder_name in sorted(pairs_by_folder):
        points = pairs_by_folder[folder_name]
        folder_x = np.array([p[0] for p in points], dtype=float)
        folder_y = np.array([p[1] for p in points], dtype=float)
        xs.extend(folder_x.tolist())
        ys.extend(folder_y.tolist())
        ax.scatter(
            folder_x,
            folder_y,
            s=80,
            color=color_for(folder_name, used_fallback),
            edgecolors="white",
            linewidths=0.8,
            label=folder_name,
            zorder=3,
        )

    all_x = np.array(xs, dtype=float)
    all_y = np.array(ys, dtype=float)
    lim_min = min(all_x.min(), all_y.min()) * 0.5
    lim_max = max(all_x.max(), all_y.max()) * 2.0
    ax.plot(
        [lim_min, lim_max],
        [lim_min, lim_max],
        "--",
        color="#94a3b8",
        linewidth=1.2,
        label="y = x",
        zorder=1,
    )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lim_min, lim_max)
    ax.set_ylim(lim_min, lim_max)
    ax.set_xlabel("Brd4ET Kd (nM)")
    ax.set_ylabel("Brd3ET Kd (nM)")
    ax.set_title("DeltaForge Kd: Brd4ET vs Brd3ET")
    ax.legend(loc="lower right", fontsize=8)
    ax.grid(True, which="both", alpha=0.25)
    ax.set_aspect("equal", adjustable="box")

    fig.tight_layout()
    fig.savefig(OUTPUT_PATH, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved plot to {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
