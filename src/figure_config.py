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
HEATMAP_CMAP_DIV = LinearSegmentedColormap.from_list("custom_diverging", [COLORS[1], "white", COLORS[3]])
HEATMAP_CMAP_RED = LinearSegmentedColormap.from_list(
    "white_to_darkred", 
    ["white", COLORS[3]]
)
HEATMAP_CMAP_BLUE = LinearSegmentedColormap.from_list(
    "white_to_darkblue", 
    ["white", COLORS[1]]
)

# ──────────────────────────────────────────────────────────────────────
# Paths
# ──────────────────────────────────────────────────────────────────────

# These are relative paths for use in notebooks (notebooks/ folder as working directory)
DATAFOLDER = Path("../data")
RESULTSFOLDER = Path("../results")
FIGSFOLDER = Path("../paper/figs/panels")

# Ensure figures folder exists
FIGSFOLDER.mkdir(parents=True, exist_ok=True)

# ──────────────────────────────────────────────────────────────────────
# Visualization Parameters
# ──────────────────────────────────────────────────────────────────────

# Heatmap settings (STATIC SCALING — hardcoded values)
# For dynamic scaling based on data distribution, use scale_vlim_to_data() from figure_plotting
# Example: vmin, vmax = scale_vlim_to_data(snips, percentile=99)
HEATMAP_VLIM_BEHAV = 4.5  # Movement data heatmap limits (static)
HEATMAP_VLIM_PHOTO = 0.5  # Photometry data heatmap limits (static)

# Alternative: Use dynamic scaling percentile
# HEATMAP_PERCENTILE_BEHAV = 99  # Use 99th percentile of absolute values
# HEATMAP_PERCENTILE_PHOTO = 99  # Use 99th percentile of absolute values

# Snips plot settings
YLIMS_BEHAV = [-0.02, 0.02]  # Movement snips y-axis limits
YLIMS_PHOTO = [-0.5, 0.5]  # Photometry snips y-axis limits

# Time series settings
TIMEPOINTS_PER_SECOND = 10  # For converting bins to seconds
INFUSION_START_BIN = 50  # Start of infusion window
INFUSION_END_BIN = 150  # End of infusion window

# ──────────────────────────────────────────────────────────────────────
# Data Processing Status
# ──────────────────────────────────────────────────────────────────────

# NOTE: These values are determined during data assembly in assemble_all_data.py.
# They are stored in the 'metadata' field of assembled_data.pickle.
# These are the DEFAULT expected values — actual values should be read from metadata.
#
# Behavioral data (movement/angular velocity):
#   - Smoothed during assembly (Gaussian smoothing, window=10)
#   - Z-scored to baseline (pre-infusion window)
#   → Do NOT smooth or z-score again in notebooks
#
# Photometry data (dopamine):
#   - NOT smoothed during assembly
#   - Z-scored by trompy during signal processing
#   → Do NOT smooth or z-score again in notebooks

# ──────────────────────────────────────────────────────────────────────
# Data Settings
# ──────────────────────────────────────────────────────────────────────

# Conditions and infusion types
CONDITIONS = ["replete", "deplete"]
INFUSION_TYPES = ["10NaCl", "45NaCl"]

# NOTE: Behavioral smoothing is done during assembly, not during plotting
# These parameters are kept for reference only and should NOT be applied to snips
BEHAV_SMOOTH_WINDOW = 5  # Reference only — behavior is already smoothed

# Photometry data is NOT smoothed, so keeping this as reference
PHOTO_SMOOTH_WINDOW = 1  # Reference only — photometry is NOT smoothed

# ──────────────────────────────────────────────────────────────────────
# Figure Control
# ──────────────────────────────────────────────────────────────────────

# Set this to True to save all figures; False to display only
SAVE_FIGS = True

# Figure file format options
SAVE_PDF = True  # Save PDF for publication
SAVE_PNG = False  # Save PNG for presentations
PNG_DPI = 300  # DPI for PNG export
