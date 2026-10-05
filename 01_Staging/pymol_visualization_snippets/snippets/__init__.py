"""
Portable PyMOL case-study visualization snippets.

Adapted from CSP_UBQ scripts/visualize.py and scripts/case_study.py.
"""

from .colors import CLASSIFICATION_COLORS, hex_to_rgb01
from .csv_adapter import MasterAlignmentData, load_master_alignment
from .render_assets import RenderAssetsConfig, render_all_variants, render_case_study_assets

__all__ = [
    "CLASSIFICATION_COLORS",
    "hex_to_rgb01",
    "MasterAlignmentData",
    "load_master_alignment",
    "RenderAssetsConfig",
    "render_case_study_assets",
    "render_all_variants",
]
