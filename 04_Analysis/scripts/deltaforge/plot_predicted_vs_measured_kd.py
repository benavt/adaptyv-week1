#!/usr/bin/env python3
"""Scatter plot of DeltaForge predicted Kd vs experimentally measured Kd."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# (target, predicted_kd_nm, measured_kd_nm)
# Predicted values from deltaforge_json/*_selected_model.json top-level kd/kd_nm.
# Not included: ET_50, ET_879, ET_522 (no measured Kd provided).
DATA = [
    ("ET_1147", 2.28, 2600.0),
    ("ET_2183", 10.11, 120_000.0),
    ("ET_2871", 71.81, 500_000.0),
    ("ET_AP2_3", 15.20, 16.26),
    ("ET_AP2_5", 105.59, 22.42),
    ("ET_AP2_6", 1.39, 953.5),
    ("ET_AP2_9", 2.36, 7.622),
    ("ET_AP2_11", 1095.30, 527.2),
    ("ET_TP", 12.55865104966462, 10),
    ("ET_951", 1.8926694054949593, 6.7),
    ("ET_1001", 4.6125985470176705, 32820),
]

# ipSAE min from fold ensemble report (first score column after the exit-code 0).
# Panel B only includes binders that also have measured Kd in DATA.
IPSAE_MIN = {
    "ET_AP2_3": 0.870016,
    "ET_AP2_9": 0.869769,
    "ET_AP2_11": 0.860319,
    "ET_AP2_5": 0.858519,
    "ET_2871": 0.822794,
    "ET_1147": 0.764927,
    "ET_2183": 0.558551,
    "ET_AP2_6": 0.086411,
}

LAI_OUT = Path(__file__).resolve().parents[3] / "analysis" / "outputs" / "ligand_ai"
OUTPUT_PATH = LAI_OUT / "predicted_vs_measured_kd.png"
OUTPUT_PATH_AB = LAI_OUT / "predicted_vs_measured_kd_ab.png"

CATEGORIES = {
    "AlphaProteo": {"color": "#eab308", "match": lambda t: "AP2" in t},
    "Perez Binder": {"color": "#2563eb", "match": lambda t: t in {"ET_2183", "ET_2871", "ET_1147"}},
    "Native Binder": {"color": "#16a34a", "match": lambda t: t in {"ET_TP", "ET_951", "ET_1001"}},
}


def category_for(target: str) -> str:
    for name, spec in CATEGORIES.items():
        if spec["match"](target):
            return name
    raise ValueError(f"No category defined for target {target}")


def _scatter_by_category(ax, targets, x, y, *, label_points: bool = True) -> None:
    for name, spec in CATEGORIES.items():
        mask = np.array([category_for(t) == name for t in targets], dtype=bool)
        if not mask.any():
            continue
        ax.scatter(
            x[mask],
            y[mask],
            s=80,
            color=spec["color"],
            edgecolors="white",
            linewidths=0.8,
            label=name,
            zorder=3,
        )
    if label_points:
        for target, xi, yi in zip(targets, x, y):
            ax.annotate(
                target,
                (xi, yi),
                textcoords="offset points",
                xytext=(6, 6),
                fontsize=9,
                color="#1e293b",
            )


def _plot_panel_a(ax) -> None:
    targets = [row[0] for row in DATA]
    predicted = np.array([row[1] for row in DATA], dtype=float)
    measured = np.array([row[2] for row in DATA], dtype=float)

    _scatter_by_category(ax, targets, measured, predicted)

    lim_min = min(measured.min(), predicted.min()) * 0.5
    lim_max = max(measured.max(), predicted.max()) * 2.0
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
    ax.set_xlabel("Measured Kd (nM)")
    ax.set_ylabel("Predicted Kd (nM)")
    ax.legend(loc="upper left")
    ax.grid(True, which="both", alpha=0.25)
    ax.set_aspect("equal", adjustable="box")


def _panel_b_arrays():
    rows = [row for row in DATA if row[0] in IPSAE_MIN]
    targets = [row[0] for row in rows]
    predicted = np.array([row[1] for row in rows], dtype=float)
    measured = np.array([row[2] for row in rows], dtype=float)
    ipsae = np.array([IPSAE_MIN[t] for t in targets], dtype=float)
    y = predicted * ipsae
    return targets, measured, y


def _plot_panel_b(ax) -> float:
    targets, measured, y = _panel_b_arrays()
    _scatter_by_category(ax, targets, measured, y)

    log_x = np.log10(measured)
    log_y = np.log10(y)
    slope, intercept = np.polyfit(log_x, log_y, 1)
    r = float(np.corrcoef(log_x, log_y)[0, 1])

    x_line = np.array([measured.min() * 0.5, measured.max() * 2.0])
    y_line = 10 ** (slope * np.log10(x_line) + intercept)
    ax.plot(
        x_line,
        y_line,
        "--",
        color="#dc2626",
        linewidth=1.4,
        label=f"trendline (r = {r:.3f})",
        zorder=2,
    )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(measured.min() * 0.5, measured.max() * 2.0)
    ax.set_ylim(y.min() * 0.5, y.max() * 2.0)
    ax.set_xlabel("Measured Kd (nM)")
    ax.set_ylabel("Predicted Kd × ipSAE min (nM)")
    ax.legend(loc="upper left")
    ax.grid(True, which="both", alpha=0.25)
    return r


def _save_single_panel() -> None:
    fig, ax = plt.subplots(figsize=(8, 8))
    _plot_panel_a(ax)
    ax.set_title("DeltaForge Predicted vs Measured Kd")
    fig.tight_layout()
    fig.savefig(OUTPUT_PATH, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved plot to {OUTPUT_PATH}")


def _save_ab_panel() -> None:
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(14, 7))
    _plot_panel_a(ax_a)
    ax_a.set_title("A. DeltaForge Predicted vs Measured Kd")
    r = _plot_panel_b(ax_b)
    ax_b.set_title("B. Predicted Kd × ipSAE min vs Measured Kd")
    fig.tight_layout()
    fig.savefig(OUTPUT_PATH_AB, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved plot to {OUTPUT_PATH_AB} (Panel B r = {r:.3f})")


def main() -> None:
    _save_single_panel()
    _save_ab_panel()


if __name__ == "__main__":
    main()
