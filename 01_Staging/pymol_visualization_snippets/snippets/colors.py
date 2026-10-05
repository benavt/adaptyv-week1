"""
TP/FP/TN/FN color definitions for PyMOL and matplotlib.

Source: scripts/config.py ClassificationColors (lines 130-136).
"""

from __future__ import annotations

from typing import Dict, Tuple

# Hex codes used in pipeline PyMOL classification scripts and bar plots.
CLASSIFICATION_COLORS: Dict[str, str] = {
    "TP": "#2ecc71",  # Green
    "FP": "#9b59b6",  # Purple
    "TN": "#3498db",  # Blue
    "FN": "#f39c12",  # Orange
}


def hex_to_rgb01(hex_color: str) -> Tuple[float, float, float]:
    """Convert ``#RRGGBB`` to PyMOL ``set_color`` RGB components in 0-1 range."""
    hex_clean = hex_color.lstrip("#")
    if len(hex_clean) != 6:
        return (0.5, 0.5, 0.5)
    return (
        int(hex_clean[0:2], 16) / 255.0,
        int(hex_clean[2:4], 16) / 255.0,
        int(hex_clean[4:6], 16) / 255.0,
    )
