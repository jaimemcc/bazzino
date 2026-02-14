"""
Configuration settings for paper figures.

This module contains shared settings used across all figure notebooks to ensure
consistency in colors, paths, and styling.
"""

from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns

# ──────────────────────────────────────────────────────────────────────
# Matplotlib Configuration
# ──────────────────────────────────────────────────────────────────────

def configure_matplotlib():
    """Set up matplotlib for publication-quality figures."""
    rcParams['font.family'] = 'Arial'
    plt.rcParams['savefig.dpi'] = 300
    plt.rcParams['savefig.transparent'] = True


# ──────────────────────────────────────────────────────────────────────
# Color Palette
# ──────────────────────────────────────────────────────────────────────

# Color scheme for conditions and infusion types
# [replete_10NaCl, replete_45NaCl, deplete_10NaCl, deplete_45NaCl]
COLORS = ["#67AFD2", "#016895", "#F4795B", "#C74632"]

# Dictionary mapping condition/infusion combinations to colors
COLOR_MAP = {
    ("replete", "10NaCl"): COLORS[0],
    ("replete", "45NaCl"): COLORS[1],
    ("deplete", "10NaCl"): COLORS[2],
    ("deplete", "45NaCl"): COLORS[3],
}

# Custom diverging colormap for heatmaps: blue (negative) - white (neutral) - red (positive)
HEATMAP_CMAP = LinearSegmentedColormap.from_list("custom_diverging", [COLORS[1], "white", COLORS[3]])

# ──────────────────────────────────────────────────────────────────────
# Paths
# ──────────────────────────────────────────────────────────────────────

# These are relative paths for use in notebooks (notebooks/ folder as working directory)
DATAFOLDER = Path("../data")
RESULTSFOLDER = Path("../results")
FIGSFOLDER = Path("../paper/figs")

# Ensure figures folder exists
FIGSFOLDER.mkdir(parents=True, exist_ok=True)

# ──────────────────────────────────────────────────────────────────────
# Visualization Parameters
# ──────────────────────────────────────────────────────────────────────

# Heatmap settings
HEATMAP_VLIM_BEHAV = 4.5  # Movement data heatmap limits
HEATMAP_VLIM_PHOTO = 0.5  # Photometry data heatmap limits

# Snips plot settings
YLIMS_BEHAV = [-3, 3]  # Movement snips y-axis limits
YLIMS_PHOTO = [-0.5, 0.5]  # Photometry snips y-axis limits

# Time series settings
TIMEPOINTS_PER_SECOND = 10  # For converting bins to seconds
INFUSION_START_BIN = 50  # Start of infusion window
INFUSION_END_BIN = 150  # End of infusion window

# ──────────────────────────────────────────────────────────────────────
# Data Settings
# ──────────────────────────────────────────────────────────────────────

# Conditions and infusion types
CONDITIONS = ["replete", "deplete"]
INFUSION_TYPES = ["10NaCl", "45NaCl"]

# Data smoothing
BEHAV_SMOOTH_WINDOW = 5  # Moving average window for behavioral data
PHOTO_SMOOTH_WINDOW = 1  # Moving average window for photometry data (no smoothing by default)

# ──────────────────────────────────────────────────────────────────────
# Figure Control
# ──────────────────────────────────────────────────────────────────────

# Set this to True to save all figures; False to display only
SAVE_FIGS = False

# Figure file format options
SAVE_PDF = True  # Save PDF for publication
SAVE_PNG = True  # Save PNG for presentations
PNG_DPI = 300  # DPI for PNG export
