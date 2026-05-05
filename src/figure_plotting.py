"""
Plotting utilities for paper figures.

This module contains reusable plotting functions for creating publication-quality
figures from assembled behavioral and photometry data.
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from pathlib import Path
from scipy import stats
from model_fit_helpers import fit_curve_series, _sigmoid_model, _exp_model

from distance_matrix import compute_mds_coordinates, GROUP_COLORS, GROUP_ORDER, GROUP_LABELS


def smooth_array(arr, window_size=5):
    """
    Smooth a 2D array along one dimension using a moving average.
    
    :param arr: 2D NumPy array (e.g., trials x timepoints)
    :param window_size: Size of the smoothing window
    :return: Smoothed 2D array with same shape as input
    """
    kernel = np.ones(window_size) / window_size
    smoothed = np.apply_along_axis(lambda m: np.convolve(m, kernel, mode='same'), axis=1, arr=arr)
    return smoothed


def conditionally_smooth(snips, should_smooth, window_size=5):
    """
    Conditionally apply smoothing based on metadata.
    
    Use this function to handle both pre-smoothed and non-smoothed data.
    Check the 'metadata' dict from assembled_data.pickle to determine if smoothing is needed.
    
    :param snips: 2D array of snips (samples x timepoints)
    :param should_smooth: Boolean; if True, apply smoothing; if False, return as-is
    :param window_size: Smoothing window size (only used if should_smooth=True)
    :return: Snips array (smoothed or original)
    """
    if should_smooth:
        return smooth_array(snips, window_size=window_size)
    else:
        return snips


def scale_vlim_to_data(snips, percentile=99):
    """
    Calculate symmetric colorbar limits based on data distribution.
    
    Useful for automatically scaling heatmap colormaps to the data range
    instead of using hard-coded limits. This ensures consistent visualization
    across different datasets and conditions.
    
    :param snips: 2D array of data (samples x timepoints)
    :param percentile: Percentile of absolute values to use as limit (default 99)
    :return: Tuple of (vmin, vmax) suitable for sns.heatmap(vmin=vmin, vmax=vmax)
    """
    vlim = np.percentile(np.abs(snips), percentile)
    return -vlim, vlim


def calculate_ylims(snips_arrays, pad_percentage=5):
    """
    Calculate automatic y-axis limits based on data range with padding.
    
    Computes the min and max across all provided snips arrays, then adds
    padding as a percentage of the range. Useful for setting consistent
    y-limits across multiple time series plots.
    
    :param snips_arrays: List of snips arrays (each is 2D: samples x timepoints)
    :param pad_percentage: Percentage of data range to use as padding (default 5%)
    :return: Tuple of (ymin, ymax) suitable for ax.set_ylim()
    """
    # Flatten all data and find global min/max
    all_data = np.concatenate([arr.flatten() for arr in snips_arrays])
    data_min = np.nanmin(all_data)
    data_max = np.nanmax(all_data)
    
    # Calculate range and add padding
    data_range = data_max - data_min
    padding = (data_range * pad_percentage) / 100
    
    ymin = data_min - padding
    ymax = data_max + padding
    
    return ymin, ymax


# ──────────────────────────────────────────────────────────────────────
# Data Extraction Functions
# ──────────────────────────────────────────────────────────────────────

def get_heatmap_data(snips, x_array, condition, infusiontype):
    """
    Extract snips for a specific condition and infusion type, averaged by trial.
    
    Creates a 2D array where each row is a trial with the trial's average across all samples.
    
    :param snips: 2D array of snips (samples x timepoints)
    :param x_array: DataFrame with trial metadata (must have columns: condition, infusiontype, trial)
    :param condition: "replete" or "deplete"
    :param infusiontype: "10NaCl" or "45NaCl"
    :return: 2D array of trial-averaged snips (trials x timepoints)
    """
    query_string = "condition == @condition & infusiontype == @infusiontype"
    
    heatmap_data = []
    df = x_array.query(query_string)
    for trial in df.trial.unique():
        tmp_snips = snips[x_array.query(query_string + " & trial == @trial").index]
        mean_snip = np.nanmean(tmp_snips, axis=0)
        heatmap_data.append(mean_snip)
        
    return np.array(heatmap_data)

def get_heatmap_data_by_rat(snips, x_array, condition, infusiontype):
    """
    Extract snips for a specific condition and infusion type, averaged by trial.
    
    Creates a 2D array where each row is a trial with the trial's average across all samples.
    
    :param snips: 2D array of snips (samples x timepoints)
    :param x_array: DataFrame with trial metadata (must have columns: condition, infusiontype, trial)
    :param condition: "replete" or "deplete"
    :param infusiontype: "10NaCl" or "45NaCl"
    :return: 2D array of trial-averaged snips (trials x timepoints)
    """
    query_string = "condition == @condition & infusiontype == @infusiontype"
    
    heatmap_data = []
    df = x_array.query(query_string)
    for id in df.id.unique():
        tmp_snips = snips[x_array.query(query_string + " & id == @id").index]
        mean_snip = np.nanmean(tmp_snips, axis=0)
        heatmap_data.append(mean_snip)
        
    return np.array(heatmap_data)


def get_mean_snips(snips, x_array, condition):
    """
    Get mean snips for each animal in a condition, separated by infusion type.
    
    Averages snips across trials for each unique animal ID within a condition.
    
    :param snips: 2D array of snips (samples x timepoints)
    :param x_array: DataFrame with trial metadata (must have columns: condition, infusiontype, id)
    :param condition: "replete" or "deplete"
    :return: Tuple of two arrays (snips_10, snips_45) of animal-averaged snips
    """
    query_string = "condition == @condition"
    
    snips_10, snips_45 = [], []
    for id in x_array.query(query_string + " & infusiontype == '10NaCl'").id.unique():
        snips_10.append(np.nanmean(snips[x_array.query(query_string + " & id == @id").index], axis=0))
    for id in x_array.query(query_string + " & infusiontype == '45NaCl'").id.unique():
        snips_45.append(np.nanmean(snips[x_array.query(query_string + " & id == @id").index], axis=0))
        
    return np.array(snips_10), np.array(snips_45)


def get_auc(snips, start_bin=50, end_bin=150):
    """
    Calculate area under the curve for snips within a time window using trapezoid rule.
    
    :param snips: 2D array of snips (samples x timepoints)
    :param start_bin: starting bin index (inclusive)
    :param end_bin: ending bin index (exclusive)
    :return: 1D array of AUC values, one per sample/animal
    """
    auc = []
    for snip in snips:
        auc.append(np.trapezoid(snip[start_bin:end_bin]))
    return np.array(auc)

def get_trial_data_by_rat(snips, x_array, condition, infusiontype, simba_col="simba_zscore_mean"):

    query_string = "condition == @condition & infusiontype == @infusiontype"

    trial_data = []
    for id in x_array.query(query_string).id.unique():
        trial_data.append(x_array.query(query_string + " & id == @id")[simba_col].values)  # Get index of first sample for this animal
    return np.array(trial_data)

# ──────────────────────────────────────────────────────────────────────
# Figure Initialization Functions
# ──────────────────────────────────────────────────────────────────────

def init_heatmap_figure():
    """
    Initialize figure with 2 heatmaps side-by-side and a shared colorbar.
    
    :return: Tuple of (fig, ax1, ax2, cbar_ax)
    """
    f = plt.figure(figsize=(2, 3.5))
    gs = f.add_gridspec(2, 2, hspace=0.1, wspace=0.05, width_ratios=[10, 1], right=0.85)
    
    ax1 = f.add_subplot(gs[0, 0])
    ax2 = f.add_subplot(gs[1, 0])
    cbar_ax = f.add_subplot(gs[0, 1])
    
    return f, ax1, ax2, cbar_ax


def init_snips_figure():
    """
    Initialize figure for time series snips plot.
    
    :return: Tuple of (fig, ax)
    """
    f = plt.figure(figsize=(2, 2))
    gs = f.add_gridspec(1, 2, hspace=0.1, wspace=0.05, width_ratios=[10, 1], bottom=0.2, right=0.85)
    ax = f.add_subplot(gs[0, 0])
    
    return f, ax


# ──────────────────────────────────────────────────────────────────────
# Plotting Functions
# ──────────────────────────────────────────────────────────────────────

def make_heatmap(data, ax, vlim, cbar_ax=None, inf_bar=False, cmap=None):
    """
    Create a heatmap on given axis with optional infusion window indicator.
    
    :param data: 2D array for heatmap (rows=trials, cols=timepoints)
    :param ax: matplotlib axis to plot on
    :param vlim: tuple of (vmin, vmax) for colorbar limits (e.g., from scale_vlim_to_data)
    :param cbar_ax: axis for colorbar (if None, no colorbar shown)
    :param inf_bar: if True, draw a small line at bottom indicating infusion window (50-150 bins)
    :param cmap: colormap to use (if None, uses default diverging colormap)
    """
    if cmap is None:
        from figure_config import HEATMAP_CMAP
        cmap = HEATMAP_CMAP
        
    has_cbar = cbar_ax is not None
    vmin, vmax = vlim  # Unpack tuple
        
    sns.heatmap(np.array(data),
                cmap=cmap,
                ax=ax,
                cbar=has_cbar,
                cbar_ax=cbar_ax,
                vmin=vmin, vmax=vmax
                )
    
    if cbar_ax is not None:
        cbar_ax.set_yticks([])
    
    if inf_bar:
        # Draw infusion window indicator bar at bottom of heatmap
        ax.plot([50, 150], [-3, -3], color="black", lw=2, alpha=0.5, clip_on=False)
        
        
    ax.set_xticks([])
    ax.set_yticks([])


def plot_snips(snips_10, snips_45, ax, colors_10, colors_45, ylims, yscalebar=None, xscalebar=False, fill_epoch=None):
    """
    Plot time series snips with error envelopes for two infusion types.
    
    Creates line plots (mean ± SEM) for each infusion type.
    
    :param snips_10: 2D array of 10NaCl snips (animals x timepoints)
    :param snips_45: 2D array of 45NaCl snips (animals x timepoints)
    :param ax: matplotlib axis to plot on
    :param colors_10: color for 10NaCl line
    :param colors_45: color for 45NaCl line
    :param ylims: [ymin, ymax] limits for y-axis
    :param scalebar: if not None, use the provided values to draw a scale bar on the left
    """
    for snips, col in zip([snips_10, snips_45], [colors_10, colors_45]):
        x = np.arange(snips.shape[1]) / 10  # Convert bins to seconds
        mean = np.mean(snips, axis=0)
        sd = np.std(snips, axis=0)
        sem = sd / np.sqrt(snips.shape[0])
        
        ax.plot(x, mean, color=col, lw=1.5)
        ax.fill_between(x, mean-sem, mean+sem, alpha=0.3, color=col)
        
    sns.despine(ax=ax, top=True, right=True, left=True, bottom=True)
    ax.set_xticks([])
    ax.set_yticks([])
    
    ax.set_xlim([0, 20])
    ax.set_ylim(ylims)
    
    # Time scale bar (5 seconds at bottom-right)
    if xscalebar:
        bar_y = ylims[0] + (ylims[1] - ylims[0]) * 0.05  # 5% above bottom of plot
        ax.plot([15, 20], [bar_y, bar_y], color="black", lw=2, alpha=0.2, clip_on=False)
    # ax.text(17.5, bar_y - (ylims[1] - ylims[0]) * 0.08, "5 s", ha="center", va="top", fontsize=9)
    
    # Value scale bar (only for deplete, typically)
    if yscalebar is not None:
        ax.plot([0, 0], [yscalebar[0], yscalebar[1]], color="black", lw=2, alpha=0.2, clip_on=False)

    if fill_epoch is not None:
        ax.axvspan(fill_epoch[0], fill_epoch[1], color="red", alpha=0.1, zorder=-10)
        

def plot_lag_peak_sharpness(
    ax,
    lag_diff_df,
    lag_diff_summary,
    color,
    mean_label="Mean ± SEM",
    individual_color="gray",
    individual_alpha=0.2,
    individual_linewidth=1,
    individual_markersize=4,
):
    """Draw lag-peak-sharpness traces on a provided axis."""
    ax.errorbar(
        lag_diff_summary["offset"].values,
        lag_diff_summary["mean_diff"].values,
        yerr=lag_diff_summary["sem_diff"].values,
        marker="o",
        capsize=5,
        capthick=2,
        linewidth=2,
        markersize=5,
        color=color,
        label=mean_label,
    )

    for animal in lag_diff_df["animal"].unique():
        animal_data = lag_diff_df.query("animal == @animal").sort_values("offset")
        ax.plot(
            animal_data["offset"],
            animal_data["r_difference"],
            "o-",
            color=individual_color,
            alpha=individual_alpha,
            linewidth=individual_linewidth,
            markersize=individual_markersize,
        )


def plot_auc_summary(
    aucs,
    colors,
    figsize=(2.2, 2.2),
    ylabel="AUC",
    sex_groups=None,
    sex_markers=None,
    add_sex_legend=False,
):
    """
    Create a bar plot with AUC summary and individual data points overlaid.
    
    :param aucs: List of two lists: [replete_aucs, deplete_aucs], where each contains
                 [aucs_10NaCl, aucs_45NaCl]
    :param colors: List of 4 colors for [rep_10, rep_45, dep_10, dep_45]
    :param figsize: Figure size (width, height)
    :param ylabel: Label for y-axis
    :param sex_groups: Optional nested list matching aucs shape with per-point sex labels
                      (e.g., [[sex_rep_10, sex_rep_45], [sex_dep_10, sex_dep_45]])
    :param sex_markers: Optional dict mapping sex label to marker (default: {"M": "^", "F": "o"})
    :param add_sex_legend: If True and sex_groups provided, add a marker legend
    :return: (fig, ax) tuple
    """
    f, ax = plt.subplots(figsize=figsize,
                         gridspec_kw={"left": 0.35, "right": 0.95, "top": 0.95, "bottom": 0.2})

    barx = [1, 2]
    barwidth = 0.35
    spacer = 0.2
    jitter_k = 0.03  # Jitter factor for individual points

    # Plot bars (means)
    ax.bar(barx[0] - spacer, np.mean(aucs[0][0]), color=colors[0], width=barwidth)
    ax.bar(barx[0] + spacer, np.mean(aucs[0][1]), color=colors[1], width=barwidth)
    ax.bar(barx[1] - spacer, np.mean(aucs[1][0]), color=colors[2], width=barwidth)
    ax.bar(barx[1] + spacer, np.mean(aucs[1][1]), color=colors[3], width=barwidth)

    if sex_markers is None:
        sex_markers = {"M": "^", "F": "o"}

    if sex_groups is None:
        # Overlay individual points (single marker style)
        ax.scatter([barx[0] - spacer]*len(aucs[0][0]) + np.random.normal(0, jitter_k, len(aucs[0][0])), aucs[0][0], facecolors="white", edgecolors=colors[0], alpha=0.5, s=30, zorder=2)
        ax.scatter([barx[0] + spacer]*len(aucs[0][1]) + np.random.normal(0, jitter_k, len(aucs[0][1])), aucs[0][1], facecolors="white", edgecolors=colors[1], alpha=0.5, s=30, zorder=2)
        ax.scatter([barx[1] - spacer]*len(aucs[1][0]) + np.random.normal(0, jitter_k, len(aucs[1][0])), aucs[1][0], facecolors="white", edgecolors=colors[2], alpha=0.8, s=30, zorder=2)
        ax.scatter([barx[1] + spacer]*len(aucs[1][1]) + np.random.normal(0, jitter_k, len(aucs[1][1])), aucs[1][1], facecolors="white", edgecolors=colors[3], alpha=0.5, s=30, zorder=2)
    else:
        group_meta = [
            (barx[0] - spacer, aucs[0][0], colors[0], 0.5, sex_groups[0][0]),
            (barx[0] + spacer, aucs[0][1], colors[1], 0.5, sex_groups[0][1]),
            (barx[1] - spacer, aucs[1][0], colors[2], 0.8, sex_groups[1][0]),
            (barx[1] + spacer, aucs[1][1], colors[3], 0.5, sex_groups[1][1]),
        ]

        for xpos, yvals, edge_col, alpha, sex_vals in group_meta:
            yvals = np.asarray(yvals)
            sex_vals = np.asarray(sex_vals)
            if len(yvals) != len(sex_vals):
                raise ValueError("sex_groups must match auc lengths for each subgroup")

            xvals = np.full(len(yvals), xpos) + np.random.normal(0, jitter_k, len(yvals))
            for sex_label, marker in sex_markers.items():
                mask = sex_vals == sex_label
                if np.any(mask):
                    ax.scatter(
                        xvals[mask],
                        yvals[mask],
                        marker=marker,
                        facecolors="white",
                        edgecolors=edge_col,
                        alpha=alpha,
                        s=30,
                        zorder=2,
                    )

            # Fallback marker for any unexpected/missing sex labels
            known_labels = np.array(list(sex_markers.keys()), dtype=object)
            unknown_mask = ~np.isin(sex_vals, known_labels)
            if np.any(unknown_mask):
                ax.scatter(
                    xvals[unknown_mask],
                    yvals[unknown_mask],
                    marker="o",
                    facecolors="white",
                    edgecolors=edge_col,
                    alpha=alpha,
                    s=30,
                    zorder=2,
                )

    # Styling
    sns.despine(ax=ax, top=True, right=True, left=False, bottom=True)
    ax.set_xticks([])
    ax.set_ylabel(ylabel, fontsize=10)

    if sex_groups is not None and add_sex_legend:
        legend_handles = []
        legend_labels = []
        for sex_label, marker in sex_markers.items():
            handle = plt.Line2D(
                [0],
                [0],
                marker=marker,
                color="none",
                markerfacecolor="white",
                markeredgecolor="black",
                markersize=5,
                linestyle="None",
            )
            legend_handles.append(handle)
            legend_labels.append("Male" if sex_label == "M" else "Female" if sex_label == "F" else str(sex_label))
        ax.legend(legend_handles, legend_labels, frameon=False, fontsize=8, loc="upper right")
    
    return f, ax


# ──────────────────────────────────────────────────────────────────────
# Utility Functions
# ──────────────────────────────────────────────────────────────────────

def save_figure(fig, filename, folder, save_pdf=True, save_png=True, png_dpi=300):
    """
    Save figure in one or both formats (PDF and PNG).
    
    :param fig: matplotlib figure object
    :param filename: filename without extension (e.g., "fig1_heatmap_replete")
    :param folder: Path object pointing to output folder
    :param save_pdf: if True, save as PDF
    :param save_png: if True, save as PNG
    :param png_dpi: DPI for PNG export
    """
    folder.mkdir(parents=True, exist_ok=True)
    
    if save_pdf:
        pdf_path = folder / f"{filename}.pdf"
        fig.savefig(pdf_path, bbox_inches='tight')
        
    if save_png:
        png_path = folder / f"{filename}.png"
        fig.savefig(png_path, bbox_inches='tight', dpi=png_dpi)


def save_figure_atomic(
    fig,
    filename,
    folder,
    save_pdf=True,
    save_png=True,
    png_dpi=300,
    temp_folder=None,
):
    """
    Save figure to a temporary folder, then move into place.

    This avoids downstream tools reading partially written files.

    :param fig: matplotlib figure object
    :param filename: filename without extension (e.g., "fig1_heatmap_replete")
    :param folder: Path object or string pointing to output folder
    :param save_pdf: if True, save as PDF
    :param save_png: if True, save as PNG
    :param png_dpi: DPI for PNG export
    :param temp_folder: optional temp folder (defaults to folder / "_tmp")
    """
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)

    if temp_folder is None:
        temp_folder = folder / "_tmp"
    temp_folder = Path(temp_folder)
    temp_folder.mkdir(parents=True, exist_ok=True)

    if save_pdf:
        temp_pdf = temp_folder / f"{filename}.pdf"
        final_pdf = folder / f"{filename}.pdf"
        fig.savefig(temp_pdf, bbox_inches='tight')
        temp_pdf.replace(final_pdf)

    if save_png:
        temp_png = temp_folder / f"{filename}.png"
        final_png = folder / f"{filename}.png"
        fig.savefig(temp_png, bbox_inches='tight', dpi=png_dpi)
        temp_png.replace(final_png)


def print_auc_stats(aucs, labels, title="Summary Statistics"):
    """
    Print formatted summary statistics for AUC data.
    
    :param aucs: List of AUC arrays (one per group)
    :param labels: List of labels for each group (e.g., ["Replete 10NaCl", "Replete 45NaCl", ...])
    :param title: Title for output section
    """
    print(f"\n{title}")
    print("=" * 60)
    for auc, label in zip(aucs, labels):
        mean = np.mean(auc)
        sem = np.std(auc) / np.sqrt(len(auc))
        print(f"{label:30s} (n={len(auc):2d}): {mean:7.2f} ± {sem:.2f}")


def draw_regression_line(y, ax, color):
    """
    Draw a linear regression line on an axis and return statistics.
    
    :param y: Array of y values (x values are assumed to be indices)
    :param ax: Matplotlib axis to plot on
    :param color: Color for the regression line
    :return: Tuple of (r_value, p_value)
    """
    x = np.arange(len(y))
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    y_fit = slope * x + intercept
    ax.plot(x, y_fit, color=color, lw=1.5)

    print(f"r = {r_value:.2f}, p = {p_value:.3f}")

    return r_value, p_value


def _fit_corr_stats(y_true, y_fit):
    """Return Pearson r and p between observed data and model predictions."""
    y_true = np.asarray(y_true, dtype=float)
    y_fit = np.asarray(y_fit, dtype=float)
    valid = np.isfinite(y_true) & np.isfinite(y_fit)
    if valid.sum() < 3:
        return np.nan, np.nan
    return stats.pearsonr(y_true[valid], y_fit[valid])


def make_correlation_plot_behav(inf10, inf45, col10, col45, yaxis=False):
    """
    Create a correlation plot showing AUC values across trials for two infusion types.
    
    :param inf10: Array of AUC values for 0.10M infusion
    :param inf45: Array of AUC values for 0.45M infusion
    :param col10: Color for 0.10M data points and fit line
    :param col45: Color for 0.45M data points and fit line
    :param yaxis: If True, show y-axis labels; if False, show tick marks only
    :return: Figure object
    """
    f, ax = plt.subplots(figsize=(1.8, 1.8),
                         gridspec_kw={"left": 0.28, "right": 0.9, "top": 0.85, "bottom": 0.24})

    ax.scatter(np.arange(len(inf10)), inf10, color=col10, alpha=0.5)
    ax.scatter(np.arange(len(inf45)), inf45, color=col45, alpha=0.5)

    r, p = draw_regression_line(inf10, ax, col10)
    if p < 0.001:
        p = "p<0.001"
    else:
        p = f"p={p:.3f}"
    ax.text(0, 1.1, f"0.10 M: r={r:.2f}, {p}", color=col10, fontsize=8,
            va="bottom", ha="left")
    
    r, p = draw_regression_line(inf45, ax, col45)
    if p < 0.001:
        p = "p<0.001"
    else:
        p = f"p={p:.3f}"
    ax.text(0, 1, f"0.45 M: r={r:.2f}, {p}", color=col45, fontsize=8,
            va="bottom", ha="left")

    sns.despine(ax=ax, offset=2)

    ax.set_ylim([-0.1, 1])
  
    if yaxis:
        ax.set_yticks([0, 0.5, 1])
        ax.set_ylabel("Time moving")
    else:
        ax.set_yticks([0, 0.5, 1], labels=["", "", ""])

    ax.set_xticks([0, 10, 20, 30, 40, 49], labels=["0", "10", "20", "30", "40", "50"])
    ax.set_xlabel("Trial")

    # ax.axhline(0, color="k", linestyle=":", alpha=0.7, zorder=-20)
    
    return f

def make_correlation_plot_simba(inf10, inf45, col10, col45, yaxis=False):
    """
    Create a correlation plot showing AUC values across trials for two infusion types.
    
    :param inf10: Array of AUC values for 0.10M infusion
    :param inf45: Array of AUC values for 0.45M infusion
    :param col10: Color for 0.10M data points and fit line
    :param col45: Color for 0.45M data points and fit line
    :param yaxis: If True, show y-axis labels; if False, show tick marks only
    :return: Figure object
    """
    f, ax = plt.subplots(figsize=(2.2, 1.8),
                         gridspec_kw={"left": 0.33, "right": 0.85, "top": 0.85, "bottom": 0.24})

    ax.scatter(np.arange(len(inf10)), inf10, color=col10, alpha=0.5)
    ax.scatter(np.arange(len(inf45)), inf45, color=col45, alpha=0.5)

    r, p = draw_regression_line(inf10, ax, col10)
    if p < 0.001:
        p = "p<0.001"
    else:
        p = f"p={p:.3f}"
    ax.text(0, 1.4, f"0.10 M: r={r:.2f}, {p}", color=col10, fontsize=8,
            va="bottom", ha="left")
    
    r, p = draw_regression_line(inf45, ax, col45)
    if p < 0.001:
        p = "p<0.001"
    else:
        p = f"p={p:.3f}"
    ax.text(0, 1.2, f"0.45 M: r={r:.2f}, {p}", color=col45, fontsize=8,
            va="bottom", ha="left")

    sns.despine(ax=ax, offset=2)

    ax.set_ylim([-1, 1])
  
    if yaxis:
        ax.set_yticks([-1, 0, 1])
        ax.set_ylabel("Appetitive Probability")
    else:
        ax.set_yticks([-1, 0, 1], labels=["", "", ""])

    ax.set_xticks([0, 10, 20, 30, 40, 49], labels=["0", "10", "20", "30", "40", "50"])
    ax.set_xlabel("Trial")

    ax.axhline(0, color="k", linestyle=":", alpha=0.7, zorder=-20)
    
    return f

def make_correlation_plot_simba_1group(inf, color, yaxis=False, fit="linear", return_stats=False, simba_metric=None):
    """
    Create a correlation plot showing AUC values across trials for two infusion types.
    
    :param inf10: Array of AUC values for 0.10M infusion
    :param inf45: Array of AUC values for 0.45M infusion
    :param col10: Color for 0.10M data points and fit line
    :param col45: Color for 0.45M data points and fit line
    :param yaxis: If True, show y-axis labels; if False, show tick marks only
    :param fit: Fit type ("linear" or "sigmoid")
    :param return_stats: If True, return (figure, r_value, p_value); otherwise return figure only
    :return: Figure object, or tuple (figure, r_value, p_value) when return_stats=True
    """
    f, ax = plt.subplots(figsize=(1.3, 1.8),
                         gridspec_kw={"left": 0.45, "right": 0.88, "top": 0.85, "bottom": 0.24})

    x = np.arange(len(inf))
    ax.scatter(x, inf, color=color, alpha=0.3)
    r_value, p_value = stats.linregress(x, inf)[2:4]

    if fit == "linear":
        draw_regression_line(inf, ax, color)
        
    elif fit == "sigmoid":
        fit_res = fit_curve_series(
            x,
            inf,
            model_name="sigmoidal",
            p0=[1, 1, 25, -1],
            bounds=([-np.inf, -np.inf, np.min(x), -np.inf], [np.inf, np.inf, np.max(x), np.inf]),
            maxfev=10000,
        )
        if fit_res["success"]:
            popt = fit_res["params"]
            x_fit = np.linspace(0, len(inf) - 1, 100)
            y_fit = _sigmoid_model(x_fit, *popt)
            r_value, p_value = fit_res["r"], fit_res["p"]
            print(f"Sigmoid fit parameters: L={popt[0]:.2f}, k={popt[1]:.2f}, x0={popt[2]:.2f}, b={popt[3]:.2f}")
            ax.plot(x_fit, y_fit, color=color, lw=1.5)
        
    elif fit == "exponential":
        fit_res = fit_curve_series(
            x,
            inf,
            model_name="exponential",
            p0=[1, -0.1, -1],
            bounds=([-np.inf, -np.inf, -np.inf], [np.inf, np.inf, np.inf]),
            maxfev=10000,
        )
        if fit_res["success"]:
            popt = fit_res["params"]
            x_fit = np.linspace(0, len(inf) - 1, 100)
            y_fit = _exp_model(x_fit, *popt)
            r_value, p_value = fit_res["r"], fit_res["p"]
            print(f"Exponential fit parameters: a={popt[0]:.2f}, b={popt[1]:.2f}, c={popt[2]:.2f}")
            ax.plot(x_fit, y_fit, color=color, lw=1.5)

    sns.despine(ax=ax, offset=2)
    
    if yaxis:
        ax.set_ylabel("Appetitive Probability")
        
    if simba_metric == "zscore":
        ax.set_ylim([-0.12, 0.2])
        ax.set_yticks([-0.1, 0, 0.1, 0.2])
        if not yaxis:
            ax.set_yticklabels(["", "", "", ""])
        p_text_y = 0.21
            
    elif simba_metric == "median":
        ax.set_ylim([-0.7, 0.9])
        ax.set_yticks([-0.5, 0, 0.5])
        if not yaxis:
            ax.set_yticklabels(["", "", ""])
        p_text_y = 0.99
    else:
        p_text_y = ax.get_ylim()[1] * 1.05
    
    ax.set_xlim([-5, 53])
    ax.set_xticks([0, 49])
    ax.set_xlabel("Trial")

    ax.axhline(0, color="k", linestyle=":", alpha=0.7, zorder=-20)
    
    if p_value < 0.001:
        p_text = "p<0.001"
    else:
        p_text = f"p={p_value:.3f}"
    ax.text(25, p_text_y, f"r={r_value:.2f}, {p_text}", color=color, fontsize=8,
            va="bottom", ha="center")
    
    if return_stats:
        return f, r_value, p_value
    return f, ax

def make_correlation_plot_da(inf10, inf45, col10, col45, yaxis=False):
    """
    Create a correlation plot showing AUC values across trials for two infusion types.
    
    :param inf10: Array of AUC values for 0.10M infusion
    :param inf45: Array of AUC values for 0.45M infusion
    :param col10: Color for 0.10M data points and fit line
    :param col45: Color for 0.45M data points and fit line
    :param yaxis: If True, show y-axis labels; if False, show tick marks only
    :return: Figure object
    """
    f, ax = plt.subplots(figsize=(2.2, 1.8),
                         gridspec_kw={"left": 0.33, "right": 0.85, "top": 0.85, "bottom": 0.24})

    ax.scatter(np.arange(len(inf10)), inf10, color=col10, alpha=0.5)
    ax.scatter(np.arange(len(inf45)), inf45, color=col45, alpha=0.5)

    r, p = draw_regression_line(inf10, ax, col10)
    if p < 0.001:
        p = "p<0.001"
    else:
        p = f"p={p:.3f}"
    ax.text(0, 200, f"0.10 M: r={r:.2f}, {p}", color=col10, fontsize=8,
            va="bottom", ha="left")
    
    r, p = draw_regression_line(inf45, ax, col45)
    if p < 0.001:
        p = "p<0.001"
    else:
        p = f"p={p:.3f}"
    ax.text(0, 175, f"0.45 M: r={r:.2f}, {p}", color=col45, fontsize=8,
            va="bottom", ha="left")

    sns.despine(ax=ax, offset=2)

    ax.set_ylim([-65, 180])
  
    if yaxis:
        ax.set_yticks([-50, 0, 50, 100, 150])
        ax.set_ylabel("Dopamine (AUC)")
    else:
        ax.set_yticks([-50, 0, 50, 100, 150], labels=["", "0", "", "", ""])

    ax.set_xticks([0, 10, 20, 30, 40, 49], labels=["0", "10", "20", "30", "40", "50"])
    ax.set_xlabel("Trial")

    ax.axhline(0, color="k", linestyle=":", alpha=0.7, zorder=-20)
    
    return f


def make_correlation_plot_da_1group(inf, color, yaxis=False, fit="linear", return_stats=False, print_stats=True, ylim=None):
    """
    Create a one-group correlation plot for dopamine (AUC) values across trials.

    :param inf: Array of dopamine (AUC) values for a single group
    :param color: Color for data points and fit line
    :param yaxis: If True, show y-axis labels; if False, show minimal tick labels
    :param fit: Fit type ("linear" or "sigmoid")
    :param return_stats: If True, return (figure, r_value, p_value); otherwise return figure only
    :return: Figure object, or tuple (figure, r_value, p_value) when return_stats=True
    """
    f, ax = plt.subplots(figsize=(1.3, 1.8),
                         gridspec_kw={"left": 0.43, "right": 0.90, "top": 0.85, "bottom": 0.24})

    x = np.arange(len(inf))
    ax.scatter(x, inf, color=color, alpha=0.3)
    r_value, p_value = stats.linregress(x, inf)[2:4]

    if fit == "linear":
        draw_regression_line(inf, ax, color)

    elif fit == "sigmoid":
        fit_res = fit_curve_series(
            x,
            inf,
            model_name="sigmoidal",
            p0=[250, 0.1, 25, -60],
            bounds=([-np.inf, -np.inf, np.min(x), -np.inf], [np.inf, np.inf, np.max(x), np.inf]),
            maxfev=10000,
        )
        if fit_res["success"]:
            popt = fit_res["params"]
            x_fit = np.linspace(0, len(inf) - 1, 100)
            y_fit = _sigmoid_model(x_fit, *popt)
            r_value, p_value = fit_res["r"], fit_res["p"]
            ax.plot(x_fit, y_fit, color=color, lw=1.5)
            print(f"Sigmoid fit parameters: L={popt[0]:.2f}, k={popt[1]:.2f}, x0={popt[2]:.2f}, b={popt[3]:.2f}")

    if print_stats:
        if p_value < 0.001:
            p_text = "p<0.001"
        else:
            p_text = f"p={p_value:.3f}"
        ax.text(20, 200, f"r={r_value:.2f}, {p_text}", color=color, fontsize=8,
                va="bottom", ha="center")

    sns.despine(ax=ax, offset=2)

    if yaxis:
        ax.set_yticks([-50, 0, 50, 100, 150])
        ax.set_ylabel("Dopamine (AUC)")
    else:
        ax.set_yticks([-50, 0, 50, 100, 150], labels=["", "", "", "", ""])
        
    if ylim is not None:
        ax.set_ylim(ylim)
    else:
        ax.set_ylim([-65, 190])

    ax.set_xticks([0, 49])
    ax.set_xlabel("Trial")
    ax.set_xlim([-5, 50])

    ax.axhline(0, color="k", linestyle=":", alpha=0.7, zorder=-20)

    if return_stats:
        return f, r_value, p_value
    
    return f, ax


def plot_lag_histogram(
    ax,
    lag_values,
    lag_range,
    colors=None,
    colors_lines=None,
    show_avg=False
):
    """
    Plot lag histogram with option for sign-based bar colors and centered integer bins.
    
    :param ax: matplotlib axis to plot on
    :param lag_values: Array of lag values to histogram
    :param lag_range: Tuple of (min_lag, max_lag) for x-axis range
    :param title: Title for the plot
    :param colors_lines: List of [color1, color2] for mean and median lines (if None, uses defaults)
    :param legend_loc: Location for legend
    :param pos_color: Color for positive lags
    :param neg_color: Color for negative lags
    :param zero_color: Color for zero lag
    """
    if colors_lines is None:
        colors_lines = ['#e74c3c', '#3498db']  # Default red and blue

    if colors == None:
        pos_color, zero_color, neg_color = '#d34d4d',  '#9a9a9a', '#3b6fb6'
    elif len(colors) == 1:
        pos_color, zero_color, neg_color = colors[0], colors[0], colors[0]
    elif len(colors) == 3:
        pos_color, zero_color, neg_color = colors
    else:
        print("Could not unpack colors properly. Make sure a list of 1 or 3 items")
    
    tick_min, tick_max = lag_range[0], lag_range[1]
    ticks = np.arange(tick_min, tick_max + 1, 1)
    bins = np.arange(tick_min - 0.5, tick_max + 1.5, 1)

    counts, bin_edges, patches = ax.hist(
        lag_values,
        bins=bins,
        edgecolor='k',
        linewidth=1.5,
    )

    for left, right, patch in zip(bin_edges[:-1], bin_edges[1:], patches):
        center = 0.5 * (left + right)
        if center == 0:
            patch.set_facecolor(zero_color)
        elif center > 0:
            patch.set_facecolor(pos_color)
        else:
            patch.set_facecolor(neg_color)
        patch.set_alpha(0.8)
    
    if show_avg:
        mean_val = np.nanmean(lag_values)
        median_val = np.nanmedian(lag_values)

        ax.axvline(mean_val, color=colors_lines[0], linestyle='--', linewidth=2.5,
                label=f'Mean: {mean_val:.2f} trials')
        ax.axvline(median_val, color=colors_lines[1], linestyle='--', linewidth=2,
                label=f'Median: {median_val:.1f} trials')
    ax.axvline(0, color='gray', linestyle=':', linewidth=1, alpha=0.5)

    sns.despine(ax=ax)
    
def make_euclidean_distance_heatmap_all_rats(distances, group_boundaries, cmap="Oranges"):
    
    heatmap_vmin = np.nanpercentile(distances, 5)
    heatmap_vmax = np.nanpercentile(distances, 95)
    cmap = plt.get_cmap(cmap).reversed()

    # All-animal heatmap
    fig, ax = plt.subplots(figsize=(1.8, 1.8))

    sns.heatmap(
        distances,
        ax=ax,
        cmap=cmap,
        vmin=heatmap_vmin,
        vmax=heatmap_vmax,
        cbar=False,
        xticklabels=False,
        yticklabels=False,
    )

    for boundary in group_boundaries[:-1]:
        ax.axhline(boundary, color="white", linewidth=1)
        ax.axvline(boundary, color="white", linewidth=1)

    return fig, ax


def make_euclidean_distance_heatmap_averaged(group_distance_matrix, cmap="Oranges"):
    
    # Group-averaged heatmap with separate colorbar
    fig, ax = plt.subplots(figsize=(1.8, 1.8))
    fig_cbar, ax_cbar = plt.subplots(figsize=(0.18, 1.8))
    cmap = plt.get_cmap(cmap).reversed()

    sns.heatmap(
        group_distance_matrix,
        ax=ax,
        cmap=cmap,
        annot=False,
        fmt=".2f",
        cbar_ax=ax_cbar,
        linewidths=1,
        linecolor="white",
    )
    ax.set_yticks([])
    ax.set_xticks([])
    ax_cbar.set_yticks([])
    
    return fig, ax, fig_cbar, ax_cbar


def make_mds_plot(coords, group_to_row_indices):
    fig, ax = plt.subplots(figsize=(1.8, 1.8),
                       gridspec_kw={"left": 0.15, "bottom": 0.15})

    group_centroids = []
    for group_key in GROUP_ORDER:
        row_idx = group_to_row_indices[group_key]
        group_coords = coords[row_idx]
        color = GROUP_COLORS[group_key]
        centroid = group_coords.mean(axis=0)
        group_centroids.append(centroid)

        for coord in group_coords:
            ax.plot([coord[0], centroid[0]], [coord[1], centroid[1]], color=color, alpha=0.3, zorder=0)

        ax.scatter(group_coords[:, 0], group_coords[:, 1], s=40, edgecolor=color, facecolor="none")
        ax.scatter(centroid[0], centroid[1], s=100, color=color)

    ax.set_xlabel("MDS1")
    ax.set_ylabel("MDS2")
    ax.set_xticks([])
    ax.set_yticks([])
    sns.despine(ax=ax, offset=5)
    
    return fig, ax


