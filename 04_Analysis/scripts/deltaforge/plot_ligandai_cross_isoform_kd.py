#!/usr/bin/env python3
"""4x4 log-log scatter of on-target vs off-target DeltaForge Kd.

Each row is one design target. The first three columns are the other isoforms.
The fourth column is that row's legend. X is the fold Kd on the isoform the
binder was designed against. Y is the fold Kd on the off-target isoform.
The generation predicted_kd column is not used.
"""

from __future__ import annotations

import csv
import html
import json
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch
import numpy as np

REPO = Path(__file__).resolve().parents[3]
OUT = REPO / "analysis" / "outputs" / "ligand_ai" / "isoform_fingerprint"
INPUT_PATH = OUT / "top_10pct_cross_isoform_kd.csv"
RFD3_LINK = OUT / "rfd3_brd4_over_brd3.csv"
RFD3_OFF_LINK = OUT / "rfd3_top5_brd2_brdt_kd.csv"
RFD3_OFF_COLUMNS = {"Brd2ET": "Brd2ET_kd_nm", "BrdTET": "BrdTET_kd_nm"}
OUTPUT_PATH = OUT / "cross_isoform_kd_pairwise.png"
HTML_PATH = OUT / "cross_isoform_kd_pairwise.html"
SPECIFICITY_PNG = OUT / "top10_specificity_kd_pairwise.png"
SPECIFICITY_HTML = OUT / "top10_specificity_kd_pairwise.html"
SPECIFICITY_TOP_N = 10
SPECIFICITY_HIGHLIGHT_N = 5
SPECIFICITY_RED = "#dc2626"
SPECIFICITY_GRAY = "#94a3b8"
ANDRE_ROOT = REPO / "design" / "ligand_ai" / "from_andre"
ANDRE_BATCHES = (
    "seq_batch_1",
    "seq_batch_2",
    "seq_batch_3",
    "seq_batch_4",
    "seq_batch_5",
)

# gene value in the CSV, fold-column prefix, panel label.
ISOFORMS = (
    ("BRD2ET", "Brd2ET", "Brd2ET"),
    ("BRD3ET", "BRD3ET", "Brd3ET"),
    ("BRD4ET", "Brd4ET", "Brd4ET"),
    ("BRDTET", "BrdTET", "BrdTET"),
)
GOLDEN = 0.61803398875


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def positive_kd(row: dict[str, str], column: str) -> float | None:
    raw = (row.get(column) or "").strip()
    if not raw:
        return None
    try:
        kd = float(raw)
    except ValueError:
        return None
    if not np.isfinite(kd) or kd <= 0:
        return None
    return kd


def load_rfd3_top5(path: Path) -> list[dict[str, str]]:
    chosen = [
        row
        for row in load_rows(path)
        if int(row["specificity_rank"]) <= 5
    ]
    chosen.sort(key=lambda row: int(row["specificity_rank"]))
    return chosen


def binders_for(rows: list[dict[str, str]], gene: str) -> list[dict[str, str]]:
    chosen = [row for row in rows if (row.get("gene") or "").strip() == gene]
    chosen.sort(key=lambda row: int(row["rank"]))
    return chosen


def color_for(rank: int) -> tuple[float, float, float]:
    hue = (rank * GOLDEN) % 1.0
    return tuple(matplotlib.colors.hsv_to_rgb((hue, 0.72, 0.9)))


def above_diagonal_all(
    binder: dict[str, str], on_prefix: str, off_prefixes: list[str]
) -> bool:
    """True when off-target Kd is higher than on-target Kd in every panel."""
    on_kd = positive_kd(binder, f"{on_prefix}_kd_nm")
    if on_kd is None or len(off_prefixes) != 3:
        return False
    for off_prefix in off_prefixes:
        off_kd = positive_kd(binder, f"{off_prefix}_kd_nm")
        if off_kd is None or off_kd <= on_kd:
            return False
    return True


def outline_specific_legend_entries(
    fig: plt.Figure, legend_rows: list[tuple[plt.Axes, list[bool]]]
) -> None:
    renderer = fig.canvas.get_renderer()
    for ax, flagged in legend_rows:
        legend = ax.get_legend()
        if legend is None:
            continue
        for handle, text, flag in zip(legend.legend_handles, legend.get_texts(), flagged):
            if not flag:
                continue
            handle_box = handle.get_window_extent(renderer).expanded(1.02, 1.05)
            text_box = text.get_window_extent(renderer).expanded(1.01, 1.05)
            extent = matplotlib.transforms.Bbox.union([handle_box, text_box])
            box = extent.transformed(ax.transAxes.inverted())
            ax.add_patch(
                FancyBboxPatch(
                    (box.x0, box.y0),
                    box.width,
                    box.height,
                    boxstyle="square,pad=0",
                    mutation_aspect=0.4,
                    transform=ax.transAxes,
                    fill=False,
                    edgecolor="#dc2626",
                    linewidth=1.1,
                    clip_on=False,
                    zorder=5,
                )
            )


def load_collaborators() -> list[dict[str, object]]:
    """One record per sequence. The first batch that contains it supplies structures."""
    chosen: dict[str, dict[str, object]] = {}
    ordered: list[dict[str, object]] = []
    labels = [label for _gene, _prefix, label in ISOFORMS]
    for batch in ANDRE_BATCHES:
        path = ANDRE_ROOT / batch / "batch.json"
        payload = json.loads(path.read_text(encoding="utf-8"))
        for item in payload["sequences"]:
            sequence = str(item["sequence"]).strip()
            if sequence in chosen:
                continue
            isoforms = item["isoforms"]
            kds: dict[str, float] = {}
            pdbs: dict[str, str] = {}
            for label in labels:
                entry = isoforms[label]
                raw = entry.get("kd_nm")
                try:
                    kd = float(raw)
                except (TypeError, ValueError):
                    kd = None
                if kd is not None and np.isfinite(kd) and kd > 0:
                    kds[label] = kd
                pdb_name = entry.get("pdb") or ""
                pdb_path = ANDRE_ROOT / batch / pdb_name
                if pdb_name and pdb_path.is_file():
                    pdbs[label] = Path(os.path.relpath(pdb_path, OUT)).as_posix()
            record = {
                "sequence": sequence,
                "kds": kds,
                "pdbs": pdbs,
            }
            chosen[sequence] = record
            ordered.append(record)
    for index, record in enumerate(ordered, start=1):
        record["binder_id"] = "C-%02d" % index
        record["color"] = hex_color(index)
    return ordered


def collaborator_is_specific(kds: object, on_label: str) -> bool:
    """True when on_label Kd is strictly lower than the other three isoforms."""
    if not isinstance(kds, dict):
        return False
    on_kd = kds.get(on_label)
    if not isinstance(on_kd, float) or on_kd <= 0:
        return False
    for _gene, _prefix, label in ISOFORMS:
        if label == on_label:
            continue
        off_kd = kds.get(label)
        if not isinstance(off_kd, float) or off_kd <= on_kd:
            return False
    return True


SPECIFICITY_COLOR = {
    "Brd2ET": "#7c3aed",
    "Brd3ET": "#16a34a",
    "Brd4ET": "#06b6d4",
    "BrdTET": "#dc2626",
}
MONO_GRAY = "#94a3b8"


def tightest_isoform(kds: dict[str, float]) -> str | None:
    for _gene, _prefix, label in ISOFORMS:
        if collaborator_is_specific(kds, label):
            return label
    return None


def mono_catalog(
    by_gene: dict[str, list[dict[str, str]]],
    rfd3_top: list[dict[str, str]],
    rfd3_off: dict[int, dict[str, str]],
    collaborators: list[dict[str, object]],
) -> list[dict[str, object]]:
    """Every unique sequence, with a Kd on all four isoforms."""
    catalog: list[dict[str, object]] = []
    seen: set[str] = set()

    def add(sequence: str, binder_id: str, symbol: str, kds: dict[str, float]) -> None:
        sequence = sequence.strip()
        if not sequence or sequence in seen:
            return
        if any(not isinstance(kds.get(label), float) for _gene, _prefix, label in ISOFORMS):
            return
        seen.add(sequence)
        preferred = tightest_isoform(kds)
        catalog.append(
            {
                "sequence": sequence,
                "binder_id": binder_id,
                "symbol": symbol,
                "kds": kds,
                "preferred": preferred,
                "color": SPECIFICITY_COLOR[preferred] if preferred else MONO_GRAY,
            }
        )

    for gene, _prefix, _label in ISOFORMS:
        for binder in by_gene[gene]:
            kds: dict[str, float] = {}
            for _gene, prefix, label in ISOFORMS:
                kd = positive_kd(binder, "%s_kd_nm" % prefix)
                if kd is not None:
                    kds[label] = kd
            add(
                binder["sequence"],
                "L-%s-%s" % (gene, int(binder["rank"])),
                "circle",
                kds,
            )
    for item in rfd3_top:
        rank = int(item["specificity_rank"])
        off = rfd3_off.get(rank, {})
        kds = {}
        for label, raw in (
            ("Brd4ET", positive_kd(item, "brd4_kd_nm")),
            ("Brd3ET", positive_kd(item, "brd3_kd_nm")),
            ("Brd2ET", positive_kd(off, "Brd2ET_kd_nm")),
            ("BrdTET", positive_kd(off, "BrdTET_kd_nm")),
        ):
            if raw is not None:
                kds[label] = raw
        add(item["binder_sequence"], "R-%s" % rank, "square", kds)
    for collab in collaborators:
        raw_kds = collab.get("kds")
        kds = {}
        if isinstance(raw_kds, dict):
            for label, kd in raw_kds.items():
                if isinstance(kd, float):
                    kds[str(label)] = kd
        add(str(collab["sequence"]), str(collab["binder_id"]), "diamond", kds)
    return catalog


def hex_color(rank: int) -> str:
    red, green, blue = color_for(rank)
    return "#%02x%02x%02x" % (round(red * 255), round(green * 255), round(blue * 255))


def write_plotly_html(
    by_gene: dict[str, list[dict[str, str]]],
    rfd3_top: list[dict[str, str]],
    rfd3_off: dict[int, dict[str, str]],
    lim_min: float,
    lim_max: float,
    collaborators: list[dict[str, object]],
) -> None:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    titles: list[str] = []
    for gene, _on_prefix, on_label in ISOFORMS:
        for off_gene, _off_prefix, off_label in ISOFORMS:
            if off_gene == gene:
                continue
            titles.append(f"{on_label} vs {off_label}")
    fig = make_subplots(
        rows=4,
        cols=3,
        subplot_titles=titles,
        horizontal_spacing=0.12,
        vertical_spacing=0.08,
    )
    legend_sections: list[str] = []
    diag_x = [lim_min, lim_max]
    diag_y = [lim_min, lim_max]

    for row_index, (gene, on_prefix, on_label) in enumerate(ISOFORMS, start=1):
        binders = by_gene[gene]
        off_targets = [item for item in ISOFORMS if item[0] != gene]
        off_prefixes = [prefix for _gene, prefix, _label in off_targets]
        flagged = {
            int(binder["rank"]): above_diagonal_all(binder, on_prefix, off_prefixes)
            for binder in binders
        }
        items: list[str] = [
            f'<h3>{html.escape(on_label)}</h3>',
        ]
        if not binders:
            items.append('<p class="empty">no %s designs</p>' % html.escape(on_label))
        for binder in binders:
            rank = int(binder["rank"])
            items.append(
                '<div class="legend-item%s" data-binder="%s">'
                '<span class="swatch" style="background:%s"></span>'
                '<span>%s  %s</span></div>'
                % (
                    " specific" if flagged[rank] else "",
                    "L-%s-%s" % (gene, rank),
                    hex_color(rank),
                    rank,
                    html.escape(binder["sequence"].strip()),
                )
            )
        if on_label == "Brd4ET":
            for item in rfd3_top:
                rank = int(item["specificity_rank"])
                items.append(
                    '<div class="legend-item" data-binder="R-%s">'
                    '<span class="swatch square" style="background:%s"></span>'
                    '<span>R%s  %s</span></div>'
                    % (
                        rank,
                        hex_color(rank + 40),
                        rank,
                        html.escape(item["binder_sequence"].strip()),
                    )
                )
        for index, collab in enumerate(collaborators, start=1):
            specific = collaborator_is_specific(collab["kds"], on_label)
            items.append(
                '<div class="legend-item%s" data-binder="%s">'
                '<span class="swatch diamond" style="background:%s"></span>'
                '<span>%s  %s</span></div>'
                % (
                    " specific" if specific else "",
                    html.escape(str(collab["binder_id"])),
                    html.escape(str(collab["color"])),
                    index,
                    html.escape(str(collab["sequence"])),
                )
            )
        legend_sections.append('<section class="legend-row">%s</section>' % "".join(items))

        for col_index, (_off_gene, off_prefix, off_label) in enumerate(off_targets, start=1):
            fig.add_trace(
                go.Scatter(
                    x=diag_x,
                    y=diag_y,
                    mode="lines",
                    line=dict(color="#94a3b8", dash="dash", width=1),
                    hoverinfo="skip",
                    showlegend=False,
                ),
                row=row_index,
                col=col_index,
            )
            xs: list[float] = []
            ys: list[float] = []
            texts: list[str] = []
            colors: list[str] = []
            symbols: list[str] = []
            line_colors: list[str] = []
            line_widths: list[float] = []
            custom: list[list[str]] = []
            for binder in binders:
                x = positive_kd(binder, f"{on_prefix}_kd_nm")
                y = positive_kd(binder, f"{off_prefix}_kd_nm")
                if x is None or y is None:
                    continue
                rank = int(binder["rank"])
                specific = flagged[rank]
                xs.append(x)
                ys.append(y)
                texts.append(str(rank))
                colors.append(hex_color(rank))
                symbols.append("circle")
                line_colors.append("#dc2626" if specific else "#ffffff")
                line_widths.append(2.0 if specific else 1.0)
                custom.append(
                    [
                        "L-%s-%s" % (gene, rank),
                        binder["sequence"].strip(),
                        binder.get("model_id") or gene,
                    ]
                )
            if on_label == "Brd4ET":
                for item in rfd3_top:
                    rank = int(item["specificity_rank"])
                    x = positive_kd(item, "brd4_kd_nm")
                    if off_label == "Brd3ET":
                        y = positive_kd(item, "brd3_kd_nm")
                    else:
                        y = positive_kd(rfd3_off.get(rank, {}), RFD3_OFF_COLUMNS[off_label])
                    if x is None or y is None:
                        continue
                    xs.append(x)
                    ys.append(y)
                    texts.append("R%s" % rank)
                    colors.append(hex_color(rank + 40))
                    symbols.append("square")
                    line_colors.append("#ffffff")
                    line_widths.append(1.0)
                    custom.append(
                        [
                            "R-%s" % rank,
                            item["binder_sequence"].strip(),
                            "RFD3 %s" % (item.get("model_id") or rank),
                        ]
                    )
            if xs:
                fig.add_trace(
                    go.Scatter(
                        x=xs,
                        y=ys,
                        mode="markers+text",
                        text=texts,
                        textposition="top right",
                        textfont=dict(size=9, color=colors),
                        marker=dict(
                            size=9,
                            color=colors,
                            symbol=symbols,
                            line=dict(color=line_colors, width=line_widths),
                        ),
                        customdata=custom,
                        hovertemplate="Sequence: %{customdata[1]}<br>Source: %{customdata[2]}<extra></extra>",
                        showlegend=False,
                    ),
                    row=row_index,
                    col=col_index,
                )
            else:
                fig.add_annotation(
                    text="no %s designs" % on_label,
                    xref="x%s" % ((row_index - 1) * 3 + col_index),
                    yref="y%s" % ((row_index - 1) * 3 + col_index),
                    x=float(np.sqrt(lim_min * lim_max)),
                    y=float(np.sqrt(lim_min * lim_max)),
                    showarrow=False,
                    font=dict(color="#64748b", size=12),
                )
            collab_x: list[float] = []
            collab_y: list[float] = []
            collab_colors: list[str] = []
            collab_line_colors: list[str] = []
            collab_line_widths: list[float] = []
            collab_custom: list[list[str]] = []
            for collab in collaborators:
                kds = collab["kds"]
                if not isinstance(kds, dict):
                    continue
                x = kds.get(on_label)
                y = kds.get(off_label)
                if not isinstance(x, float) or not isinstance(y, float):
                    continue
                specific = collaborator_is_specific(kds, on_label)
                collab_x.append(x)
                collab_y.append(y)
                collab_colors.append(str(collab["color"]))
                collab_line_colors.append("#dc2626" if specific else "#ffffff")
                collab_line_widths.append(2.0 if specific else 1.0)
                collab_custom.append(
                    [
                        str(collab["binder_id"]),
                        str(collab["sequence"]),
                        "collaborator",
                    ]
                )
            if collab_x:
                fig.add_trace(
                    go.Scatter(
                        x=collab_x,
                        y=collab_y,
                        mode="markers",
                        marker=dict(
                            size=9,
                            color=collab_colors,
                            symbol="diamond",
                            line=dict(
                                color=collab_line_colors,
                                width=collab_line_widths,
                            ),
                        ),
                        customdata=collab_custom,
                        hovertemplate="Sequence: %{customdata[1]}<br>Source: %{customdata[2]}<extra></extra>",
                        showlegend=False,
                    ),
                    row=row_index,
                    col=col_index,
                )
            fig.update_xaxes(
                type="log",
                range=[np.log10(lim_min), np.log10(lim_max)],
                title_text="%s Kd (nM)" % on_label,
                row=row_index,
                col=col_index,
            )
            fig.update_yaxes(
                type="log",
                range=[np.log10(lim_min), np.log10(lim_max)],
                title_text="%s Kd (nM)" % off_label,
                row=row_index,
                col=col_index,
            )

    fig.update_layout(
        height=1400,
        width=1100,
        margin=dict(l=70, r=20, t=40, b=40),
        hovermode="closest",
        plot_bgcolor="white",
    )
    catalog = mono_catalog(by_gene, rfd3_top, rfd3_off, collaborators)
    mono = make_subplots(
        rows=4,
        cols=3,
        subplot_titles=titles,
        horizontal_spacing=0.12,
        vertical_spacing=0.08,
    )
    source_for = {"circle": "LigandAI", "square": "RFD3", "diamond": "collaborator"}
    for row_index, (_gene, _on_prefix, on_label) in enumerate(ISOFORMS, start=1):
        off_targets = [item for item in ISOFORMS if item[0] != _gene]
        for col_index, (_off_gene, _off_prefix, off_label) in enumerate(off_targets, start=1):
            mono.add_trace(
                go.Scatter(
                    x=diag_x,
                    y=diag_y,
                    mode="lines",
                    line=dict(color="#94a3b8", dash="dash", width=1),
                    hoverinfo="skip",
                    showlegend=False,
                ),
                row=row_index,
                col=col_index,
            )
            for symbol in ("circle", "square", "diamond"):
                xs: list[float] = []
                ys: list[float] = []
                colors: list[str] = []
                line_colors: list[str] = []
                line_widths: list[float] = []
                custom: list[list[str]] = []
                for rec in catalog:
                    if rec["symbol"] != symbol:
                        continue
                    kds = rec["kds"]
                    if not isinstance(kds, dict):
                        continue
                    x = kds.get(on_label)
                    y = kds.get(off_label)
                    if not isinstance(x, float) or not isinstance(y, float):
                        continue
                    xs.append(x)
                    ys.append(y)
                    colors.append(str(rec["color"]))
                    line_colors.append("#ffffff")
                    line_widths.append(1.0)
                    custom.append(
                        [
                            str(rec["binder_id"]),
                            str(rec["sequence"]),
                            source_for[symbol],
                        ]
                    )
                if not xs:
                    continue
                mono.add_trace(
                    go.Scatter(
                        x=xs,
                        y=ys,
                        mode="markers",
                        marker=dict(
                            size=9,
                            color=colors,
                            symbol=symbol,
                            line=dict(color=line_colors, width=line_widths),
                        ),
                        customdata=custom,
                        hovertemplate="Sequence: %{customdata[1]}<br>Source: %{customdata[2]}<extra></extra>",
                        showlegend=False,
                    ),
                    row=row_index,
                    col=col_index,
                )
            mono.update_xaxes(
                type="log",
                range=[np.log10(lim_min), np.log10(lim_max)],
                title_text="%s Kd (nM)" % on_label,
                row=row_index,
                col=col_index,
            )
            mono.update_yaxes(
                type="log",
                range=[np.log10(lim_min), np.log10(lim_max)],
                title_text="%s Kd (nM)" % off_label,
                row=row_index,
                col=col_index,
            )
    mono.update_layout(
        height=1400,
        width=1100,
        margin=dict(l=70, r=20, t=40, b=40),
        hovermode="closest",
        plot_bgcolor="white",
    )
    swatch_class = {"circle": "", "square": " square", "diamond": " diamond"}
    mono_items = [
        '<div class="mono-key">%s</div>'
        % "".join(
            '<span style="color:%s">%s</span>' % (color, html.escape(label))
            for label, color in SPECIFICITY_COLOR.items()
        )
    ]
    for rec in catalog:
        color = str(rec["color"])
        mono_items.append(
            '<div class="legend-item" data-binder="%s">'
            '<span class="swatch%s" style="background:%s"></span>'
            '<span style="color:%s">%s</span></div>'
            % (
                html.escape(str(rec["binder_id"])),
                swatch_class[str(rec["symbol"])],
                color,
                color,
                html.escape(str(rec["sequence"])),
            )
        )
    mono_legend = '<section class="legend-row">%s</section>' % "".join(mono_items)
    structure_map: dict[str, dict[str, str]] = {}
    pdb_dir = OUT / "structures"
    for item in rfd3_top:
        rank = int(item["specificity_rank"])
        files: dict[str, str] = {}
        for isoform in ("Brd2ET", "Brd3ET", "Brd4ET", "BrdTET"):
            link = pdb_dir / ("R%s_%s.pdb" % (rank, isoform))
            files[isoform] = ("structures/%s" % link.name) if link.is_file() else ""
        structure_map["R-%s" % rank] = files
    for collab in collaborators:
        pdbs = collab["pdbs"]
        files = {}
        for isoform in ("Brd2ET", "Brd3ET", "Brd4ET", "BrdTET"):
            files[isoform] = str(pdbs.get(isoform) or "") if isinstance(pdbs, dict) else ""
        structure_map[str(collab["binder_id"])] = files
    rainbow_div = fig.to_html(full_html=False, include_plotlyjs=True, div_id="kd-plot-rainbow")
    mono_div = mono.to_html(full_html=False, include_plotlyjs=False, div_id="kd-plot-mono")
    page = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Cross-isoform Kd</title>
<style>
body { margin: 0; font-family: sans-serif; background: #fff; color: #0f172a; }
.toolbar { display: flex; gap: 8px; padding: 8px 8px 0; }
.toolbar button { font: inherit; font-size: 14px; padding: 4px 12px; border: 1px solid #cbd5e1; background: #fff; border-radius: 4px; cursor: pointer; }
.toolbar button.active { border-color: #0f172a; font-weight: 700; }
.wrap { display: flex; align-items: stretch; }
.wrap.hidden { display: none; }
.legends { width: 420px; display: flex; flex-direction: column; }
.legend-row { flex: 1; overflow: auto; padding: 4px 8px; }
.legend-row h3 { margin: 4px 0; font-size: 16px; }
.mono-key { display: flex; flex-wrap: wrap; gap: 10px; font-size: 12px; font-weight: 700; margin: 4px 0 8px; }
.legend-item { display: flex; gap: 6px; align-items: center; font-size: 11px; line-height: 1.2; padding: 1px 2px; }
.legend-item.specific { box-shadow: inset 0 0 0 1px #dc2626; }
.legend-item.hover { box-shadow: inset 0 0 0 2px #facc15; }
.swatch { width: 9px; height: 9px; border-radius: 50%%; flex: none; border: 1px solid #fff; }
.swatch.square { border-radius: 1px; }
.swatch.diamond { width: 11px; height: 11px; border-radius: 0; clip-path: polygon(50%% 0, 100%% 50%%, 50%% 100%%, 0 50%%); }
.empty { color: #64748b; }
.structures { display: flex; gap: 8px; padding: 8px; }
.structure-col { flex: 1; min-width: 0; }
.structure-col h3 { margin: 0 0 4px; font-size: 16px; }
.viewer { height: 280px; border: 1px solid #e2e8f0; position: relative; }
.note { color: #64748b; font-size: 12px; min-height: 16px; }
</style>
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
</head>
<body>
<div class="toolbar">
  <button type="button" id="btn-rainbow" class="active">Rainbow</button>
  <button type="button" id="btn-mono">Monochrome</button>
</div>
<div class="wrap" id="view-rainbow">
%s
<div class="legends" id="legends-rainbow">
%s
</div>
</div>
<div class="wrap hidden" id="view-mono">
%s
<div class="legends" id="legends-mono">
%s
</div>
</div>
<div class="structures">
  <div class="structure-col"><h3>Brd2ET</h3><div id="view-Brd2ET" class="viewer"></div><div id="note-Brd2ET" class="note"></div></div>
  <div class="structure-col"><h3>Brd3ET</h3><div id="view-Brd3ET" class="viewer"></div><div id="note-Brd3ET" class="note"></div></div>
  <div class="structure-col"><h3>Brd4ET</h3><div id="view-Brd4ET" class="viewer"></div><div id="note-Brd4ET" class="note"></div></div>
  <div class="structure-col"><h3>BrdTET</h3><div id="view-BrdTET" class="viewer"></div><div id="note-BrdTET" class="note"></div></div>
</div>
<script>
var STRUCTURES = %s;
var ISOFORMS = ["Brd2ET", "Brd3ET", "Brd4ET", "BrdTET"];
var viewers = {};
var selected = null;
ISOFORMS.forEach(function (name) {
  viewers[name] = $3Dmol.createViewer("view-" + name, {backgroundColor: "white"});
});
function clearViewers() {
  ISOFORMS.forEach(function (name) {
    viewers[name].clear();
    viewers[name].render();
    document.getElementById("note-" + name).textContent = "";
  });
}
function showBinder(binderId) {
  if (selected === binderId) {
    clearViewers();
    selected = null;
    return;
  }
  selected = binderId;
  var files = STRUCTURES[binderId] || {};
  ISOFORMS.forEach(function (name) {
    var url = files[name] || "";
    var note = document.getElementById("note-" + name);
    viewers[name].clear();
    if (!url) {
      note.textContent = "not downloaded";
      viewers[name].render();
      return;
    }
    note.textContent = "";
    fetch(url).then(function (response) { return response.text(); }).then(function (pdb) {
      if (selected !== binderId) return;
      viewers[name].addModel(pdb, "pdb");
      viewers[name].setStyle({}, {cartoon: {color: "spectrum"}});
      viewers[name].zoomTo();
      viewers[name].render();
    });
  });
}
var plots = {
  rainbow: document.getElementById('kd-plot-rainbow'),
  mono: document.getElementById('kd-plot-mono')
};
var bases = { rainbow: {}, mono: {} };
function snapshotOne(name) {
  var gd = plots[name];
  var base = {};
  gd.data.forEach(function (trace, i) {
    if (!trace.customdata || !trace.marker || !trace.marker.line) return;
    var color = trace.marker.line.color;
    var width = trace.marker.line.width;
    base[i] = {
      color: Array.isArray(color) ? color.slice() : color,
      width: Array.isArray(width) ? width.slice() : width
    };
  });
  bases[name] = base;
}
function snapshot() {
  snapshotOne('rainbow');
  snapshotOne('mono');
}
function restore() {
  Object.keys(bases).forEach(function (name) {
    var gd = plots[name];
    var base = bases[name];
    Object.keys(base).forEach(function (key) {
      var i = Number(key);
      Plotly.restyle(gd, {
        'marker.line.color': [base[i].color.slice ? base[i].color.slice() : base[i].color],
        'marker.line.width': [base[i].width.slice ? base[i].width.slice() : base[i].width]
      }, [i]);
    });
  });
  document.querySelectorAll('.legend-item.hover').forEach(function (el) {
    el.classList.remove('hover');
  });
}
function highlight(binderId) {
  restore();
  if (!binderId) return;
  Object.keys(bases).forEach(function (name) {
    var gd = plots[name];
    var base = bases[name];
    gd.data.forEach(function (trace, i) {
      if (!trace.customdata || !base[i]) return;
      var colors = base[i].color.slice();
      var widths = base[i].width.slice();
      var hit = false;
      trace.customdata.forEach(function (cd, j) {
        if (cd && cd[0] === binderId) {
          colors[j] = '#facc15';
          widths[j] = 3.5;
          hit = true;
        }
      });
      if (hit) {
        Plotly.restyle(gd, {'marker.line.color': [colors], 'marker.line.width': [widths]}, [i]);
      }
    });
  });
  document.querySelectorAll('.legend-item[data-binder="' + binderId + '"]').forEach(function (item) {
    item.classList.add('hover');
  });
}
function bindPlot(name) {
  var gd = plots[name];
  gd.on('plotly_hover', function (ev) {
    var cd = ev.points && ev.points[0] && ev.points[0].customdata;
    if (cd) highlight(cd[0]);
  });
  gd.on('plotly_unhover', function () { restore(); });
  gd.on('plotly_click', function (ev) {
    var cd = ev.points && ev.points[0] && ev.points[0].customdata;
    if (cd) showBinder(cd[0]);
  });
}
function showView(name) {
  restore();
  document.getElementById('view-rainbow').classList.toggle('hidden', name !== 'rainbow');
  document.getElementById('view-mono').classList.toggle('hidden', name !== 'mono');
  document.getElementById('btn-rainbow').classList.toggle('active', name === 'rainbow');
  document.getElementById('btn-mono').classList.toggle('active', name === 'mono');
  Plotly.Plots.resize(plots[name]);
}
document.getElementById('btn-rainbow').addEventListener('click', function () { showView('rainbow'); });
document.getElementById('btn-mono').addEventListener('click', function () { showView('mono'); });
bindPlot('rainbow');
bindPlot('mono');
document.querySelectorAll('.legend-item').forEach(function (el) {
  el.addEventListener('click', function () {
    var binderId = el.getAttribute('data-binder');
    if (binderId) showBinder(binderId);
  });
});
snapshot();
</script>
</body>
</html>
""" % (
        rainbow_div,
        "".join(legend_sections),
        mono_div,
        mono_legend,
        json.dumps(structure_map),
    )
    HTML_PATH.write_text(page, encoding="utf-8")
    print("Saved plot to %s (%s monochrome sequences)" % (HTML_PATH, len(catalog)))


def isoform_kds(row: dict[str, str]) -> dict[str, float] | None:
    """Positive fold Kd on every isoform, or None when any isoform is missing."""
    kds: dict[str, float] = {}
    for _gene, prefix, label in ISOFORMS:
        kd = positive_kd(row, "%s_kd_nm" % prefix)
        if kd is None:
            return None
        kds[label] = kd
    return kds


def specificity_score(kds: dict[str, float], label: str) -> float:
    """Mean log10 distance above y=x against the other three isoforms."""
    gaps = [
        np.log10(kds[other]) - np.log10(kds[label])
        for _gene, _prefix, other in ISOFORMS
        if other != label
    ]
    return float(np.mean(gaps))


def top_specificity(rows: list[dict[str, str]]) -> dict[str, list[dict[str, object]]]:
    """Top sequences per isoform by mean log10 Kd gap versus the other three."""
    complete = []
    for row in rows:
        kds = isoform_kds(row)
        if kds is None:
            continue
        complete.append((row, kds))
    ranked: dict[str, list[dict[str, object]]] = {}
    for _gene, _prefix, label in ISOFORMS:
        ordered = sorted(
            complete,
            key=lambda item: specificity_score(item[1], label),
            reverse=True,
        )[:SPECIFICITY_TOP_N]
        ranked[label] = [
            {
                "rank": index,
                "sequence": row["sequence"].strip(),
                "kds": kds,
                "highlight": index <= SPECIFICITY_HIGHLIGHT_N,
            }
            for index, (row, kds) in enumerate(ordered, start=1)
        ]
    return ranked


def specificity_color(rank: int) -> str:
    if rank <= SPECIFICITY_HIGHLIGHT_N:
        return SPECIFICITY_RED
    return SPECIFICITY_GRAY


def specificity_legend_line(label: str, item: dict[str, object]) -> str:
    return "%s : %s : %s" % (label, item["rank"], item["sequence"])


def print_specificity_legend(ranked: dict[str, list[dict[str, object]]]) -> None:
    for _gene, _prefix, label in ISOFORMS:
        for item in ranked[label]:
            print(specificity_legend_line(label, item))


def write_specificity_png(
    ranked: dict[str, list[dict[str, object]]],
    lim_min: float,
    lim_max: float,
) -> None:
    fig, axes = plt.subplots(
        4,
        4,
        figsize=(18, 16),
        sharex="col",
        sharey="col",
        gridspec_kw={"width_ratios": [1, 1, 1, 1.7]},
    )
    for row_index, (_gene, _prefix, on_label) in enumerate(ISOFORMS):
        binders = ranked[on_label]
        off_targets = [item for item in ISOFORMS if item[2] != on_label]
        for col_index, (_off_gene, _off_prefix, off_label) in enumerate(off_targets):
            ax = axes[row_index, col_index]
            ax.plot(
                [lim_min, lim_max],
                [lim_min, lim_max],
                "--",
                color="#94a3b8",
                linewidth=1.0,
                zorder=1,
            )
            for binder in binders:
                kds = binder["kds"]
                if not isinstance(kds, dict):
                    continue
                x = kds[on_label]
                y = kds[off_label]
                rank = int(binder["rank"])
                color = specificity_color(rank)
                ax.scatter(
                    [x],
                    [y],
                    s=46,
                    color=color,
                    edgecolors="white",
                    linewidths=0.6,
                    zorder=3,
                )
                ax.annotate(
                    str(rank),
                    (x, y),
                    textcoords="offset points",
                    xytext=(4, 3),
                    fontsize=8,
                    color=color,
                    zorder=4,
                )
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlim(lim_min, lim_max)
            ax.set_ylim(lim_min, lim_max)
            ax.set_title("%s vs %s" % (on_label, off_label), fontsize=10)
            ax.set_xlabel("%s Kd (nM)" % on_label)
            ax.set_ylabel("%s Kd (nM)" % off_label)
            ax.grid(True, which="both", alpha=0.25)
            ax.set_aspect("equal", adjustable="box")

        legend = axes[row_index, 3]
        legend.axis("off")
        handles = [
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="none",
                markerfacecolor=specificity_color(int(binder["rank"])),
                markeredgecolor="white",
                markeredgewidth=0.6,
                markersize=6,
                label="%s  %s" % (binder["rank"], binder["sequence"]),
            )
            for binder in binders
        ]
        legend.legend(
            handles=handles,
            loc="center left",
            frameon=False,
            fontsize=7,
            handletextpad=0.4,
            borderaxespad=0.2,
            labelspacing=0.35,
            title=on_label,
            title_fontproperties={"size": 12, "weight": "bold"},
        )
    fig.tight_layout(w_pad=2.4)
    fig.savefig(SPECIFICITY_PNG, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("Saved plot to %s" % SPECIFICITY_PNG)


def write_specificity_html(
    ranked: dict[str, list[dict[str, object]]],
    lim_min: float,
    lim_max: float,
) -> None:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    titles: list[str] = []
    for _gene, _prefix, on_label in ISOFORMS:
        for _off_gene, _off_prefix, off_label in ISOFORMS:
            if off_label == on_label:
                continue
            titles.append("%s vs %s" % (on_label, off_label))
    fig = make_subplots(
        rows=4,
        cols=3,
        subplot_titles=titles,
        horizontal_spacing=0.08,
        vertical_spacing=0.08,
    )
    diag = [lim_min, lim_max]
    legend_sections: list[str] = []
    for row_index, (_gene, _prefix, on_label) in enumerate(ISOFORMS, start=1):
        binders = ranked[on_label]
        items = ["<h3>%s</h3>" % html.escape(on_label)]
        for binder in binders:
            rank = int(binder["rank"])
            color = specificity_color(rank)
            items.append(
                '<div class="legend-item" style="color:%s">%s</div>'
                % (color, html.escape(specificity_legend_line(on_label, binder)))
            )
        legend_sections.append('<section class="legend-row">%s</section>' % "".join(items))
        off_targets = [item for item in ISOFORMS if item[2] != on_label]
        for col_index, (_off_gene, _off_prefix, off_label) in enumerate(off_targets, start=1):
            fig.add_trace(
                go.Scatter(
                    x=diag,
                    y=diag,
                    mode="lines",
                    line=dict(color="#94a3b8", dash="dash", width=1),
                    hoverinfo="skip",
                    showlegend=False,
                ),
                row=row_index,
                col=col_index,
            )
            xs: list[float] = []
            ys: list[float] = []
            texts: list[str] = []
            colors: list[str] = []
            custom: list[list[str]] = []
            for binder in binders:
                kds = binder["kds"]
                if not isinstance(kds, dict):
                    continue
                xs.append(float(kds[on_label]))
                ys.append(float(kds[off_label]))
                texts.append(str(binder["rank"]))
                colors.append(specificity_color(int(binder["rank"])))
                custom.append([str(binder["sequence"]), str(binder["rank"])])
            fig.add_trace(
                go.Scatter(
                    x=xs,
                    y=ys,
                    mode="markers+text",
                    text=texts,
                    textposition="top right",
                    textfont=dict(size=11, color=colors),
                    marker=dict(size=11, color=colors, line=dict(color="white", width=1)),
                    customdata=custom,
                    hovertemplate="Rank %{customdata[1]}<br>%{customdata[0]}<extra></extra>",
                    showlegend=False,
                ),
                row=row_index,
                col=col_index,
            )
            fig.update_xaxes(
                type="log",
                range=[np.log10(lim_min), np.log10(lim_max)],
                title_text="%s Kd (nM)" % on_label,
                row=row_index,
                col=col_index,
            )
            fig.update_yaxes(
                type="log",
                range=[np.log10(lim_min), np.log10(lim_max)],
                title_text="%s Kd (nM)" % off_label,
                row=row_index,
                col=col_index,
            )
    fig.update_layout(
        height=1400,
        width=980,
        margin=dict(l=70, r=20, t=40, b=40),
        hovermode="closest",
        plot_bgcolor="white",
    )
    plot_div = fig.to_html(full_html=False, include_plotlyjs="cdn")
    page = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Top 10 isoform-specific Kd</title>
<style>
body { font-family: sans-serif; margin: 0; color: #0f172a; }
.layout { display: flex; align-items: flex-start; }
.plot { flex: 1 1 auto; min-width: 0; }
.legend { width: 420px; padding: 16px 16px 16px 0; }
.legend-row { margin-bottom: 28px; }
.legend-row h3 { margin: 0 0 6px; font-size: 14px; }
.legend-item { font-size: 12px; line-height: 1.45; font-family: ui-monospace, monospace; }
</style>
</head>
<body>
<div class="layout">
<div class="plot">%s</div>
<aside class="legend">%s</aside>
</div>
</body>
</html>
""" % (
        plot_div,
        "".join(legend_sections),
    )
    SPECIFICITY_HTML.write_text(page, encoding="utf-8")
    print("Saved plot to %s" % SPECIFICITY_HTML)


def write_specificity_figure(rows: list[dict[str, str]]) -> None:
    ranked = top_specificity(rows)
    values: list[float] = []
    for binders in ranked.values():
        for binder in binders:
            kds = binder["kds"]
            if isinstance(kds, dict):
                values.extend(float(kd) for kd in kds.values())
    lim_min, lim_max = shared_limits(values)
    write_specificity_png(ranked, lim_min, lim_max)
    write_specificity_html(ranked, lim_min, lim_max)
    print_specificity_legend(ranked)


def shared_limits(values: list[float]) -> tuple[float, float]:
    if not values:
        raise SystemExit(f"No positive Kd pairs in {INPUT_PATH}")
    lo = min(values)
    hi = max(values)
    log_lo = np.log10(lo)
    log_hi = np.log10(hi)
    pad = 0.08 * (log_hi - log_lo) if log_hi > log_lo else 0.5
    return 10 ** (log_lo - pad), 10 ** (log_hi + pad)


def main() -> None:
    rows = load_rows(INPUT_PATH)
    rfd3_top = load_rfd3_top5(RFD3_LINK)
    rfd3_off = {int(row["specificity_rank"]): row for row in load_rows(RFD3_OFF_LINK)}
    kd_values: list[float] = []
    for item in rfd3_top:
        for column in ("brd4_kd_nm", "brd3_kd_nm"):
            kd = positive_kd(item, column)
            if kd is not None:
                kd_values.append(kd)
    for row in rfd3_off.values():
        for column in RFD3_OFF_COLUMNS.values():
            kd = positive_kd(row, column)
            if kd is not None:
                kd_values.append(kd)
    by_gene: dict[str, list[dict[str, str]]] = {}
    for gene, _prefix, _label in ISOFORMS:
        binders = binders_for(rows, gene)
        by_gene[gene] = binders
        for row in binders:
            for _off_gene, off_prefix, _off_label in ISOFORMS:
                kd = positive_kd(row, f"{off_prefix}_kd_nm")
                if kd is not None:
                    kd_values.append(kd)

    lim_min, lim_max = shared_limits(kd_values)
    collaborators = load_collaborators()
    plotly_values = list(kd_values)
    for collab in collaborators:
        kds = collab["kds"]
        if isinstance(kds, dict):
            plotly_values.extend(kd for kd in kds.values() if isinstance(kd, float))
    plotly_min, plotly_max = shared_limits(plotly_values)
    fig, axes = plt.subplots(
        4,
        4,
        figsize=(18, 16),
        sharex="col",
        sharey="col",
        gridspec_kw={"width_ratios": [1, 1, 1, 1.25]},
    )
    legend_rows: list[tuple[plt.Axes, list[bool]]] = []

    for row_index, (gene, on_prefix, on_label) in enumerate(ISOFORMS):
        binders = by_gene[gene]
        off_targets = [item for item in ISOFORMS if item[0] != gene]
        off_prefixes = [prefix for _gene, prefix, _label in off_targets]
        flagged = [above_diagonal_all(binder, on_prefix, off_prefixes) for binder in binders]
        flagged_by_rank = {
            int(binder["rank"]): flag for binder, flag in zip(binders, flagged)
        }
        for col_index, (_off_gene, off_prefix, off_label) in enumerate(off_targets):
            ax = axes[row_index, col_index]
            ax.plot(
                [lim_min, lim_max],
                [lim_min, lim_max],
                "--",
                color="#94a3b8",
                linewidth=1.0,
                zorder=1,
            )
            plotted = 0
            for binder in binders:
                x = positive_kd(binder, f"{on_prefix}_kd_nm")
                y = positive_kd(binder, f"{off_prefix}_kd_nm")
                if x is None or y is None:
                    continue
                rank = int(binder["rank"])
                color = color_for(rank)
                specific = flagged_by_rank.get(rank, False)
                ax.scatter(
                    [x],
                    [y],
                    s=42,
                    color=color,
                    edgecolors="#dc2626" if specific else "white",
                    linewidths=1.6 if specific else 0.6,
                    zorder=4 if specific else 3,
                )
                ax.annotate(
                    str(rank),
                    (x, y),
                    textcoords="offset points",
                    xytext=(4, 3),
                    fontsize=7,
                    color=color,
                    zorder=4,
                )
                plotted += 1
            if on_label == "Brd4ET" and off_label in ("Brd2ET", "Brd3ET", "BrdTET"):
                for item in rfd3_top:
                    rank = int(item["specificity_rank"])
                    x = positive_kd(item, "brd4_kd_nm")
                    if off_label == "Brd3ET":
                        y = positive_kd(item, "brd3_kd_nm")
                    else:
                        y = positive_kd(rfd3_off.get(rank, {}), RFD3_OFF_COLUMNS[off_label])
                    if x is None or y is None:
                        continue
                    color = color_for(rank + 40)
                    ax.scatter(
                        [x],
                        [y],
                        s=48,
                        marker="s",
                        color=color,
                        edgecolors="white",
                        linewidths=0.6,
                        zorder=5,
                    )
                    ax.annotate(
                        f"R{rank}",
                        (x, y),
                        textcoords="offset points",
                        xytext=(4, 3),
                        fontsize=7,
                        color=color,
                        zorder=6,
                    )
                    plotted += 1
            if plotted == 0:
                ax.text(
                    0.5,
                    0.5,
                    f"no {on_label} designs",
                    transform=ax.transAxes,
                    ha="center",
                    va="center",
                    color="#64748b",
                    fontsize=9,
                )
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlim(lim_min, lim_max)
            ax.set_ylim(lim_min, lim_max)
            ax.set_title(f"{on_label} vs {off_label}", fontsize=10)
            ax.set_xlabel(f"{on_label} Kd (nM)")
            ax.set_ylabel(f"{off_label} Kd (nM)")
            ax.grid(True, which="both", alpha=0.25)
            ax.set_aspect("equal", adjustable="box")
            print(f"{on_label} vs {off_label}: {plotted} points")

        legend = axes[row_index, 3]
        legend.axis("off")
        if not binders:
            legend.text(
                0.0,
                0.5,
                f"no {on_label} designs",
                transform=legend.transAxes,
                ha="left",
                va="center",
                color="#64748b",
                fontsize=9,
            )
            continue
        handles = [
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="none",
                markerfacecolor=color_for(int(binder["rank"])),
                markeredgecolor="white",
                markeredgewidth=0.6,
                markersize=6,
                label=f"{int(binder['rank'])}  {binder['sequence'].strip()}",
            )
            for binder in binders
        ]
        legend_flags = list(flagged)
        if on_label == "Brd4ET":
            for item in rfd3_top:
                rank = int(item["specificity_rank"])
                handles.append(
                    Line2D(
                        [0],
                        [0],
                        marker="s",
                        linestyle="none",
                        markerfacecolor=color_for(rank + 40),
                        markeredgecolor="white",
                        markeredgewidth=0.6,
                        markersize=6,
                        label=f"R{rank}  {item['binder_sequence'].strip()}",
                    )
                )
                legend_flags.append(False)
        legend.legend(
            handles=handles,
            loc="center left",
            frameon=False,
            fontsize=6,
            handletextpad=0.4,
            borderaxespad=0.2,
            labelspacing=0.25,
            title=on_label,
            title_fontproperties={"size": 12, "weight": "bold"},
        )
        legend_rows.append((legend, legend_flags))
        ranks = [int(binder["rank"]) for binder, flag in zip(binders, flagged) if flag]
        print(f"{on_label} legend: {len(handles)} binders, above y=x in all three: {ranks}")

    fig.tight_layout(w_pad=3.2)
    fig.subplots_adjust(wspace=0.72)
    fig.canvas.draw()
    outline_specific_legend_entries(fig, legend_rows)
    fig.savefig(OUTPUT_PATH, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved plot to {OUTPUT_PATH}")
    write_plotly_html(by_gene, rfd3_top, rfd3_off, plotly_min, plotly_max, collaborators)
    write_specificity_figure(rows)


if __name__ == "__main__":
    main()
