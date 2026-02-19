"""
Plotting utilities for paper figures.

This module contains reusable plotting functions for creating publication-quality
figures from assembled behavioral and photometry data.
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


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
        auc.append(np.trapz(snip[start_bin:end_bin]))
    return np.array(auc)


# ──────────────────────────────────────────────────────────────────────
# Figure Initialization Functions
# ──────────────────────────────────────────────────────────────────────

def init_heatmap_figure():
    """
    Initialize figure with 2 heatmaps side-by-side and a shared colorbar.
    
    :return: Tuple of (fig, ax1, ax2, cbar_ax)
    """
    f = plt.figure(figsize=(2, 3.5))
    gs = f.add_gridspec(2, 2, hspace=0.1, wspace=0.05, width_ratios=[10, 1])
    
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
    gs = f.add_gridspec(1, 2, hspace=0.1, wspace=0.05, width_ratios=[10, 1], bottom=0.2)
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


def plot_snips(snips_10, snips_45, ax, colors_10, colors_45, ylims, scalebar=False):
    """
    Plot time series snips with error envelopes for two infusion types.
    
    Creates line plots (mean ± SEM) for each infusion type.
    
    :param snips_10: 2D array of 10NaCl snips (animals x timepoints)
    :param snips_45: 2D array of 45NaCl snips (animals x timepoints)
    :param ax: matplotlib axis to plot on
    :param colors_10: color for 10NaCl line
    :param colors_45: color for 45NaCl line
    :param ylims: [ymin, ymax] limits for y-axis
    :param scalebar: if True, draw a 1-unit scale bar at top-left
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
    bar_y = ylims[0] + (ylims[1] - ylims[0]) * 0.05  # 5% above bottom of plot
    ax.plot([15, 20], [bar_y, bar_y], color="black", lw=2, alpha=0.5, clip_on=False)
    ax.text(17.5, bar_y - (ylims[1] - ylims[0]) * 0.08, "5 s", ha="center", va="top", fontsize=9)
    
    # Value scale bar (only for deplete, typically)
    if scalebar:
        ax.plot([0, 0], [1, 2], color="black", lw=2, alpha=0.5, clip_on=False)


def plot_auc_summary(aucs, colors, figsize=(2.2, 2.2), ylabel="AUC"):
    """
    Create a bar plot with AUC summary and individual data points overlaid.
    
    :param aucs: List of two lists: [replete_aucs, deplete_aucs], where each contains
                 [aucs_10NaCl, aucs_45NaCl]
    :param colors: List of 4 colors for [rep_10, rep_45, dep_10, dep_45]
    :param figsize: Figure size (width, height)
    :param ylabel: Label for y-axis
    :return: (fig, ax) tuple
    """
    f, ax = plt.subplots(figsize=figsize,
                         gridspec_kw={"left": 0.3, "right": 0.95, "top": 0.95, "bottom": 0.2})

    barx = [1, 2]
    barwidth = 0.35
    spacer = 0.2

    # Plot bars (means)
    ax.bar(barx[0] - spacer, np.mean(aucs[0][0]), color=colors[0], width=barwidth)
    ax.bar(barx[0] + spacer, np.mean(aucs[0][1]), color=colors[1], width=barwidth)
    ax.bar(barx[1] - spacer, np.mean(aucs[1][0]), color=colors[2], width=barwidth)
    ax.bar(barx[1] + spacer, np.mean(aucs[1][1]), color=colors[3], width=barwidth)

    # Overlay individual points
    ax.scatter([barx[0] - spacer]*len(aucs[0][0]), aucs[0][0], facecolors="white", edgecolors=colors[0], alpha=0.5, s=30)
    ax.scatter([barx[0] + spacer]*len(aucs[0][1]), aucs[0][1], facecolors="white", edgecolors=colors[1], alpha=0.5, s=30)
    ax.scatter([barx[1] - spacer]*len(aucs[1][0]), aucs[1][0], facecolors="white", edgecolors=colors[2], alpha=0.8, s=30)
    ax.scatter([barx[1] + spacer]*len(aucs[1][1]), aucs[1][1], facecolors="white", edgecolors=colors[3], alpha=0.5, s=30)

    # Styling
    sns.despine(ax=ax, top=True, right=True, left=False, bottom=True)
    ax.set_xticks([])
    ax.set_ylabel(ylabel, fontsize=10)
    
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
