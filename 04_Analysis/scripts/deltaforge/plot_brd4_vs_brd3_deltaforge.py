#!/usr/bin/env python3
"""Scatter Brd4ET vs Brd3ET DeltaForge predictions for the top 0.5% sequences."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parents[3]
RFD3_OUT = REPO / "analysis" / "outputs" / "rfd3"
INPUT_PATH = RFD3_OUT / "top_0.5pct_sequences_deltaforge.csv"
SPECIFICITY_PATH = RFD3_OUT / "top_0.5pct_sequences_deltaforge_brd4_over_brd3.csv"
OUTPUT_PATH = RFD3_OUT / "brd4et_vs_brd3et_predicted.png"
FINGERPRINT_OUTPUT_PATH = RFD3_OUT / "brd4et_vs_brd3et_predicted_isoform_fingerprint.png"
SELECTED_PATH = (
    REPO / "design" / "rfd3" / "best_Brd4ET_binders" / "best_Brd4ET_binders.csv"
)
SELECTED_OUTPUT_PATH = RFD3_OUT / "brd4et_vs_brd3et_predicted_selected.png"

PANELS = (
    ("kd_nm", "Predicted Kd (nM)", "log"),
    ("delta_g", "Predicted ΔG", "symlog"),
    ("iptm", "ipTM", "linear"),
    ("ptm", "pTM", "linear"),
    ("peptide_ipsae", "Peptide ipSAE", "linear"),
    ("mean_plddt", "Mean pLDDT", "linear"),
)


def load_unique_rows(path: Path) -> list[dict[str, str]]:
    seen: set[str] = set()
    rows: list[dict[str, str]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            sequence = row["sequence"]
            if sequence in seen:
                continue
            seen.add(sequence)
            rows.append(row)
    return rows


def load_selected(path: Path) -> list[tuple[str, str]]:
    """Return (rank, sequence) in CSV order."""
    with path.open(newline="", encoding="utf-8") as handle:
        return [(row["rank"], row["sequence"]) for row in csv.DictReader(handle)]


def load_top_specific(path: Path, n: int = 10) -> list[tuple[str, str]]:
    """Return (specificity rank, sequence) for the first n ranked rows."""
    chosen: list[tuple[int, str]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            rank = int(row["specificity_rank"])
            if rank <= n:
                chosen.append((rank, row["sequence"]))
    chosen.sort(key=lambda item: item[0])
    return [(str(rank), sequence) for rank, sequence in chosen]


def pairs_for(rows: list[dict[str, str]], field: str) -> tuple[np.ndarray, np.ndarray]:
    x = np.array([float(row[f"brd4_{field}"]) for row in rows], dtype=float)
    y = np.array([float(row[f"brd3_{field}"]) for row in rows], dtype=float)
    return x, y


def _padded_limits(x: np.ndarray, y: np.ndarray, scale: str) -> tuple[float, float]:
    lo = float(min(x.min(), y.min()))
    hi = float(max(x.max(), y.max()))
    if scale == "log":
        log_lo = np.log10(lo)
        log_hi = np.log10(hi)
        pad = 0.08 * (log_hi - log_lo)
        return 10 ** (log_lo - pad), 10 ** (log_hi + pad)
    if scale == "symlog":
        # Values are negative. Pad the log of the magnitude so both ends stay signed.
        mag_small = abs(hi)
        mag_large = abs(lo)
        log_small = np.log10(mag_small)
        log_large = np.log10(mag_large)
        pad = 0.08 * (log_large - log_small)
        return -(10 ** (log_large + pad)), -(10 ** (log_small - pad))
    span = hi - lo
    pad = 0.05 * span if span else 0.05
    return lo - pad, hi + pad


def _plot_panel(
    ax,
    x: np.ndarray,
    y: np.ndarray,
    title: str,
    scale: str,
    *,
    highlight: tuple[np.ndarray, np.ndarray, list[str]] | None = None,
    show_legend: bool = False,
    highlight_label: str = "Selected Brd4ET binders",
    marked: tuple[np.ndarray, np.ndarray] | None = None,
    marked_label: str = "isoform_fingerprint",
) -> None:
    ax.scatter(
        x,
        y,
        s=28,
        color="#2563eb",
        edgecolors="white",
        linewidths=0.4,
        alpha=0.45 if highlight is not None or marked is not None else 0.8,
        zorder=3,
        label="Other binders" if highlight is not None or marked is not None else None,
    )
    if marked is not None:
        mx, my = marked
        ax.scatter(
            mx,
            my,
            s=42,
            marker="D",
            color="#d97706",
            edgecolors="#78350f",
            linewidths=0.5,
            alpha=0.9,
            zorder=4,
            label=marked_label,
        )
    if highlight is not None:
        hx, hy, labels = highlight
        ax.scatter(
            hx,
            hy,
            s=90,
            color="#dc2626",
            edgecolors="#1e293b",
            linewidths=0.8,
            zorder=4,
            label=highlight_label,
        )
        for label, xi, yi in zip(labels, hx, hy):
            ax.annotate(
                label,
                (xi, yi),
                textcoords="offset points",
                xytext=(5, 5),
                fontsize=9,
                fontweight="bold",
                color="#991b1b",
                zorder=5,
            )
    lim_min, lim_max = _padded_limits(x, y, scale)
    ax.plot(
        [lim_min, lim_max],
        [lim_min, lim_max],
        "--",
        color="#94a3b8",
        linewidth=1.2,
        zorder=1,
    )
    if scale == "log":
        ax.set_xscale("log")
        ax.set_yscale("log")
    elif scale == "symlog":
        ax.set_xscale("symlog", linthresh=1)
        ax.set_yscale("symlog", linthresh=1)
    ax.set_xlim(lim_min, lim_max)
    ax.set_ylim(lim_min, lim_max)
    ax.set_box_aspect(1)
    ax.set_xlabel(f"Brd4ET {title}")
    ax.set_ylabel(f"Brd3ET {title}")
    ax.set_title(title)
    ax.grid(True, which="both", alpha=0.25)
    if show_legend:
        ax.legend(loc="upper left", fontsize=8, framealpha=0.9)


def _save_figure(
    rows: list[dict[str, str]],
    output_path: Path,
    selected: list[tuple[str, str]] | None = None,
    highlight_label: str = "Selected Brd4ET binders",
    marked: list[dict[str, str]] | None = None,
    marked_label: str = "isoform_fingerprint",
) -> None:
    by_sequence = {row["sequence"]: row for row in rows}
    marked_sequences = {row["sequence"] for row in marked} if marked else set()
    selected_sequences = {sequence for _, sequence in selected} if selected else set()
    background = [
        row
        for row in rows
        if row["sequence"] not in marked_sequences and row["sequence"] not in selected_sequences
    ]
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 9))
    for ax, (field, title, scale) in zip(axes.ravel(), PANELS):
        x, y = pairs_for(background, field)
        highlight = None
        if selected is not None:
            chosen = [by_sequence[sequence] for _, sequence in selected]
            labels = [rank for rank, _ in selected]
            highlight = (*pairs_for(chosen, field), labels)
        marked_pairs = pairs_for(marked, field) if marked else None
        _plot_panel(
            ax,
            x,
            y,
            title,
            scale,
            highlight=highlight,
            show_legend=(selected is not None or marked is not None) and field == "kd_nm",
            highlight_label=highlight_label,
            marked=marked_pairs,
            marked_label=marked_label,
        )
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved plot to {output_path}")


def main() -> None:
    rows = load_unique_rows(INPUT_PATH)
    top = load_top_specific(SPECIFICITY_PATH) if SPECIFICITY_PATH.is_file() else None
    if top:
        missing = [sequence for _rank, sequence in top if sequence not in {row["sequence"] for row in rows}]
        if missing:
            raise SystemExit(f"{len(missing)} highlighted sequences are missing from {INPUT_PATH}")
    _save_figure(
        rows,
        OUTPUT_PATH,
        top,
        highlight_label="Top 10 Brd4ET-specific",
    )
    if SELECTED_PATH.exists():
        _save_figure(rows, SELECTED_OUTPUT_PATH, load_selected(SELECTED_PATH))
    fingerprint = [row for row in rows if "isoform_fingerprint" in row["run"]]
    _save_figure(
        rows,
        FINGERPRINT_OUTPUT_PATH,
        top,
        highlight_label="Top 10 Brd4ET-specific",
        marked=fingerprint,
        marked_label="isoform_fingerprint",
    )


if __name__ == "__main__":
    main()
