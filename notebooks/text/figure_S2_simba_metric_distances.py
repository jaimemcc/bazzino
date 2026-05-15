# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.2
#   kernelspec:
#     display_name: default
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Figure S2 SIMBA Metric Distances
#
# This notebook uses a scalar SIMBA metric throughout, rather than trial-window means from the snips.
#
# Supported metric modes:
# - `median_balance`: uses `simba_median_balance_score`
# - `mean_zscore_auc`: uses `simba_mean_zscore_auc`
#
# The main output is the animal-by-animal distance matrix and downstream visualizations.

# %%
import sys
import warnings
from pathlib import Path

import dill
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.ndimage import gaussian_filter1d
from scipy.spatial.distance import cdist
from sklearn.manifold import MDS

from trompy import save_figure_atomic

sys.path.insert(0, str(Path("../src").resolve()))

from pickle_compat import enable_dill_pathlib_compat
from figure_config import configure_matplotlib, COLORS, DATAFOLDER, FIGSFOLDER, HEATMAP_CMAP_DIV

from figure_plotting import (
    smooth_array, get_heatmap_data_by_rat, get_mean_snips, get_auc,
    init_heatmap_figure, init_snips_figure, make_heatmap,
    plot_snips, plot_auc_summary, print_auc_stats, get_trial_data_by_rat,
    scale_vlim_to_data, calculate_ylims,
    draw_regression_line, make_correlation_plot_simba
)

from figure_config import COLORS, FIGSFOLDER, DATAFOLDER

ROOT = Path.cwd()
if ROOT.name == "notebooks":
    ROOT = ROOT.parent

SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

# Use absolute paths so notebook behavior is robust to kernel cwd.
DATAFOLDER = (ROOT / "data").resolve()
FIGSFOLDER = (ROOT / "paper" / "figs" / "panels").resolve()
FIGSFOLDER.mkdir(parents=True, exist_ok=True)

enable_dill_pathlib_compat()
configure_matplotlib()

SIMBA_SCALAR_METRIC = "da_auc"
DISTANCE_METRIC = "euclidean"
APPLY_TRIAL_SMOOTHING = True
SMOOTH_SIGMA = 1.2
SAVE_FIGS = True

GROUP_ORDER = [
    ("deplete", "10NaCl"),
    ("deplete", "45NaCl"),
    ("replete", "10NaCl"),
    ("replete", "45NaCl"),
]

GROUP_LABELS = {
    ("replete", "10NaCl"): "Replete 10NaCl",
    ("replete", "45NaCl"): "Replete 45NaCl",
    ("deplete", "10NaCl"): "Deplete 10NaCl",
    ("deplete", "45NaCl"): "Deplete 45NaCl",
}

GROUP_COLORS = {group: COLORS[idx] for idx, group in zip([2, 3, 0, 1], GROUP_ORDER)}
custom_cmap = HEATMAP_CMAP_DIV  # Use shared diverging colormap

# %%
assembled_data_path = DATAFOLDER / "assembled_data.pickle"

with open(assembled_data_path, "rb") as f:
    data = dill.load(f)

x_array = data["x_array"].copy()
params = data.get("params", {})
metadata = data.get("metadata", {})

metric_options = {
    "median_balance": {
        "column": "simba_median_balance",
        "label": "SIMBA median-balance score",
        "file_stub": "median_balance",
    },
    "mean_zscore_auc": {
        "column": "simba_zscore_mean",
        "label": "SIMBA mean z-score AUC",
        "file_stub": "mean_zscore_auc",
    },
    "da_auc": {
        "column": "auc_snips",
        "label": "Dopamine AUC",
        "file_stub": "da_auc",
    },
}

snips_simba = data["snips_simba_baseline"]

if SIMBA_SCALAR_METRIC not in metric_options:
    raise ValueError(f"Unknown SIMBA_SCALAR_METRIC: {SIMBA_SCALAR_METRIC}")

metric_config = metric_options[SIMBA_SCALAR_METRIC]
metric_col = metric_config["column"]
metric_label = metric_config["label"]
metric_file_stub = metric_config["file_stub"]

required_cols = ["id", "condition", "infusiontype", "trial", metric_col]
missing_cols = [col for col in required_cols if col not in x_array.columns]
if missing_cols:
    raise KeyError(f"x_array is missing required columns: {missing_cols}")

print(f"Loaded {assembled_data_path}")
print(f"Using metric column: {metric_col}")
print(f"Metric label: {metric_label}")
print(f"Rows in x_array: {len(x_array)}")
print(f"Metric range: {x_array[metric_col].min():.4f} to {x_array[metric_col].max():.4f}")

# %%
scalar_df = x_array.loc[:, ["id", "condition", "infusiontype", "trial", metric_col]].copy()
scalar_df = scalar_df.dropna(subset=[metric_col])
scalar_df["group_key"] = list(zip(scalar_df["condition"], scalar_df["infusiontype"]))
scalar_df = scalar_df[scalar_df["group_key"].isin(GROUP_ORDER)].copy()
scalar_df["group_key"] = pd.Categorical(scalar_df["group_key"], categories=GROUP_ORDER, ordered=True)
scalar_df = scalar_df.sort_values(["group_key", "id", "trial"])

trial_counts = scalar_df.groupby(["group_key", "id"]).trial.nunique().sort_index()
print("Trials per animal/group series:")
print(trial_counts.groupby(level=0).agg(["count", "min", "max"]))

metric_matrix_df = scalar_df.pivot_table(
    index=["group_key", "id"],
    columns="trial",
    values=metric_col,
    aggfunc="mean",
).sort_index()

if metric_matrix_df.isna().any().any():
    missing_by_series = metric_matrix_df.isna().sum(axis=1)
    missing_by_series = missing_by_series[missing_by_series > 0]
    raise ValueError(
        "Some animal/group series are missing trial values. Missing counts:\n"
        f"{missing_by_series.to_string()}"
    )

metric_matrix = metric_matrix_df.to_numpy(dtype=float)
if APPLY_TRIAL_SMOOTHING:
    metric_matrix_for_distance = gaussian_filter1d(metric_matrix, sigma=SMOOTH_SIGMA, axis=1)
else:
    metric_matrix_for_distance = metric_matrix.copy()

series_index = list(metric_matrix_df.index)
series_labels = [f"{GROUP_LABELS[group_key]} | {animal_id}" for group_key, animal_id in series_index]
group_names = [GROUP_LABELS[group_key] for group_key, _ in series_index]
group_to_row_indices = {
    group_key: [idx for idx, (series_group, _) in enumerate(series_index) if series_group == group_key]
    for group_key in GROUP_ORDER
}
group_sizes = {GROUP_LABELS[group_key]: len(indices) for group_key, indices in group_to_row_indices.items()}
group_boundaries = np.cumsum([len(group_to_row_indices[group_key]) for group_key in GROUP_ORDER])

print(f"Metric matrix shape: {metric_matrix.shape}")
print(f"Distance input shape: {metric_matrix_for_distance.shape}")
print("Series per group:")
print(group_sizes)

metric_matrix_df.head()

# %%
snips_to_plot = np.array(smooth_array(snips_simba))  # Use smoothed and z-scored data for all subsequent analyses and plotting

# Parameters for visualization - MOVEMENT
# Use dynamic scaling based on actual data ranges (asymmetric, not symmetric)
# This ensures heatmaps use the full range of your data without artificial symmetry
vmin = np.nanpercentile(snips_to_plot, 5)  # 5th percentile for lower bound
vmax = np.nanpercentile(snips_to_plot, 95)  # 95th percentile for upper bound
print(f"Calculated asymmetric vlims: vmin={vmin:.4f}, vmax={vmax:.4f}")

vlim = (-0.2, 0.2)

# Calculate dynamic y-limits based on all snips data - MOVEMENT
print("\n" + "="*60)
print("CALCULATING Y-LIMITS FOR TIME SERIES PLOTS")
print("="*60)

# Get all condition/infusion combinations to compute limits
rep_10, rep_45 = get_mean_snips(snips_simba, x_array, "replete")
dep_10, dep_45 = get_mean_snips(snips_simba, x_array, "deplete")

# Calculate y-limits based on all snips with 5% padding
calc_ylims = calculate_ylims([rep_10, rep_45, dep_10, dep_45], pad_percentage=5)

print(f"  Data range (min, max): ({np.nanmin(snips_simba):.4f}, {np.nanmax(snips_simba):.4f})")
print(f"  Calculated y-limits: {calc_ylims}")
print("="*60 + "\n")

# Use calculated limits or override with fixed values if preferred
ylims = calc_ylims
ylims = (0, 8)


# %%
### 2A-D Heatmaps and time series of SIMBA snips by rat by group

heatmap_data_rep_10 = get_heatmap_data_by_rat(snips_to_plot, x_array, "replete", "10NaCl")
heatmap_data_rep_45 = get_heatmap_data_by_rat(snips_to_plot, x_array, "replete", "45NaCl")
heatmap_data_dep_10 = get_heatmap_data_by_rat(snips_to_plot, x_array, "deplete", "10NaCl")
heatmap_data_dep_45 = get_heatmap_data_by_rat(snips_to_plot, x_array, "deplete", "45NaCl")

f, ax = plt.subplots(ncols=5, figsize=(7, 1.8),
                     gridspec_kw={"width_ratios": [1, 1, 1, 1, 0.1],
                                  "left": 0.1, "right": 0.95, "wspace": 0.1})

sns.heatmap(heatmap_data_rep_10, ax=ax[0], cmap=custom_cmap, vmin=vlim[0], vmax=vlim[1], cbar=False)
sns.heatmap(heatmap_data_rep_45, ax=ax[1], cmap=custom_cmap, vmin=vlim[0], vmax=vlim[1], cbar=False)
sns.heatmap(heatmap_data_dep_10, ax=ax[2], cmap=custom_cmap, vmin=vlim[0], vmax=vlim[1], cbar=False)
sns.heatmap(heatmap_data_dep_45, ax=ax[3], cmap=custom_cmap, vmin=vlim[0], vmax=vlim[1], cbar=True, cbar_ax=ax[4])



for axis in ax:
    axis.set_xticks([])
    axis.set_yticks([])
    axis.plot((149, 198), (10.5, 10.5), color="black", lw=2, alpha=0.5, clip_on=False)  # Add time scale bar (0.5s at 300Hz)
    axis.plot((49, 148), (-0.5, -0.5), color="black", lw=2, alpha=0.5, clip_on=False)  # Add time scale bar (0.5s at 300Hz)

ax[0].set_yticks([0.5, 9.5], labels=["1", "10"], fontsize=10, rotation=0)
ax[0].set_ylabel("Rats")

if SAVE_FIGS:
    save_figure_atomic(f, "figure_S2_simba_groups_heatmaps", FIGSFOLDER)

# %%
fig, axes = plt.subplots(1, 5, figsize=(7, 1.8), sharey=True,
                         gridspec_kw={"width_ratios": [1, 1, 1, 1, 0.1],
                                      "left": 0.1, "right": 0.95, "bottom": 0.3, "wspace": 0.1})
trial_numbers = metric_matrix_df.columns.to_numpy()

for ax, group_key in zip(axes[:4], GROUP_ORDER):
    row_indices = group_to_row_indices[group_key]
    group_series = metric_matrix[row_indices]
    color = GROUP_COLORS[group_key]

    for series in group_series:
        ax.plot(trial_numbers, gaussian_filter1d(series, sigma=SMOOTH_SIGMA), color=color, alpha=0.2, linewidth=1)

    ax.plot(trial_numbers, np.nanmean(group_series, axis=0), color=color, linewidth=2)
    ax.axhline(0, color="gray", linestyle="--", linewidth=1)
    sns.despine(ax=ax, offset=5)
    ax.set_xlabel("Trial")
    ax.set_xticks([0, 49])

axes[0].set_ylabel("Appetitive probability")
axes[4].axis("off")

if SAVE_FIGS:
    save_figure_atomic(fig, f"figure_S2_simba_{metric_file_stub}_trajectories", folder=FIGSFOLDER)


# %%
distances = cdist(metric_matrix_for_distance, metric_matrix_for_distance, metric=DISTANCE_METRIC)
distances_df = pd.DataFrame(distances, index=series_labels, columns=series_labels)

group_distance_matrix = pd.DataFrame(
    index=[GROUP_LABELS[group_key] for group_key in GROUP_ORDER],
    columns=[GROUP_LABELS[group_key] for group_key in GROUP_ORDER],
    dtype=float,
)

for row_group in GROUP_ORDER:
    row_idx = group_to_row_indices[row_group]
    for col_group in GROUP_ORDER:
        col_idx = group_to_row_indices[col_group]
        group_distance_matrix.loc[GROUP_LABELS[row_group], GROUP_LABELS[col_group]] = distances[np.ix_(row_idx, col_idx)].mean()

print(f"Distance metric: {DISTANCE_METRIC}")
print(f"Mean off-diagonal distance: {distances[~np.eye(distances.shape[0], dtype=bool)].mean():.4f}")
group_distance_matrix

# %%
heatmap_vmin = np.nanpercentile(distances, 5)
heatmap_vmax = np.nanpercentile(distances, 95)
cmap = plt.get_cmap("Oranges").reversed()

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

if SAVE_FIGS:
    save_figure_atomic(fig, f"figure_S2_simba_{metric_file_stub}_distance_heatmap_all", folder=FIGSFOLDER)

# Group-averaged heatmap with separate colorbar
fig2, ax2 = plt.subplots(figsize=(1.8, 1.8))
fig2_cbar, ax2_cbar = plt.subplots(figsize=(0.18, 1.8))

sns.heatmap(
    group_distance_matrix,
    ax=ax2,
    cmap=cmap,
    annot=False,
    fmt=".2f",
    cbar_ax=ax2_cbar,
    linewidths=1,
    linecolor="white",
)
ax2.set_yticks([])
ax2.set_xticks([])

if SAVE_FIGS:
    save_figure_atomic(fig2, f"figure_S2_simba_{metric_file_stub}_distance_heatmap", folder=FIGSFOLDER)
    save_figure_atomic(fig2_cbar, f"figure_S2_simba_{metric_file_stub}_distance_heatmap_colorbar", folder=FIGSFOLDER)

plt.show()


# %%
distance_scale = np.nanmax(distances)

if distance_scale > 0:
    distances_for_mds = distances / distance_scale
else:
    distances_for_mds = distances.copy()

with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning, module="sklearn.manifold._mds")
    mds = MDS(
        n_components=2,
        metric=True,
        dissimilarity="precomputed",
        random_state=42,
        n_init=8,
        max_iter=500,
        init="random",
    )
    coords = mds.fit_transform(distances_for_mds)

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

if SAVE_FIGS:
    save_figure_atomic(fig, f"figure_S2_simba_{metric_file_stub}_distance_mds", folder=FIGSFOLDER)

plt.show()


# %%
nearest_neighbors = []
for idx, label in enumerate(series_labels):
    row = distances[idx].copy()
    row[idx] = np.inf
    neighbor_idx = int(np.argmin(row))
    nearest_neighbors.append({
        "series": label,
        "group": group_names[idx],
        "nearest_neighbor": series_labels[neighbor_idx],
        "nearest_neighbor_group": group_names[neighbor_idx],
        "distance": row[neighbor_idx],
    })

nearest_neighbors_df = pd.DataFrame(nearest_neighbors).sort_values(["group", "distance"])
nearest_neighbors_df.head(20)

# %%
