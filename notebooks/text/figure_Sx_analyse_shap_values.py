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

# %%
# %load_ext autoreload
# %autoreload 2

from pathlib import Path
import sys
import re

# Register dill/pathlib compatibility shim BEFORE importing dill
sys.path.insert(0, str(Path("../src").resolve()))
from pickle_compat import enable_dill_pathlib_compat
enable_dill_pathlib_compat()

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import dill

from matplotlib.colors import to_hex
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

from trompy import save_figure_atomic

from figure_config import (
    configure_matplotlib, COLORS, HEATMAP_CMAP_DIV,
    DATAFOLDER, RESULTSFOLDER, FIGSFOLDER,
    SAVE_FIGS
)

# from utils import make_realigned_trials
from realignment_helpers import get_realigned_data

SAVE_FIGS = True

# Configure matplotlib
configure_matplotlib()
colors = COLORS  # Use shared color palette
custom_cmap = HEATMAP_CMAP_DIV  # Use shared colormap

# 11 evenly spaced colors from custom_cmap
sampled_hex = [to_hex(custom_cmap(x)) for x in np.linspace(0, 1, 11)]
DA_COLOR = sampled_hex[0]  # Choose the 3rd color for DA
BEHAV_COLOR = sampled_hex[-1]  # Choose the 6th color for behavior


# %%
from assemble_shap_dfs import assemble_shap_dfs, assign_group

DATAFOLDER = Path("../data/shap_values")
df_shap, df_shap_raw, feature_summary = assemble_shap_dfs(DATAFOLDER)

# Compatibility views used by the exploratory plots below.
# feature_groups = feature_summary[["feature", "importance", "group"]].copy()
# cumul_summary = feature_summary.copy()

# %%
feature_summary


# %%
# Plot all features with cumulative SHAP importance on x and feature rows on y.

bodypart_colors = {
    "nose": "#d95f02",
    "tail base": "#1b9e77",
    "head base": "#7570b3",
    "ears": "#e7298a",
    "all body parts": "#e6ab02",
    "whole mouse": "#a6761d",
    "other": "#bdbdbd",
}
bodypart_codes = {name: code for code, name in enumerate(bodypart_colors)}
bodypart_values = feature_summary["bodypart"].map(bodypart_codes).to_numpy()

window_order = ["none", "2", "5", "6", "7.5", "10"]
window_colors = ["#d9d9d9", "#440154", "#31688e", "#35b779", "#90d743", "#fde725"]
window_codes = {name: code for code, name in enumerate(window_order)}
window_values = feature_summary["timewindow"].map(window_codes).to_numpy()

f, [ax1, ax_window, ax2] = plt.subplots(
    figsize=(5, 5),
    ncols=3,
    sharey=True,
    gridspec_kw={"width_ratios": [0.05, 0.05, 1], "wspace": 0.08},
)

ax1.imshow(
    bodypart_values[:, None],
    aspect="auto",
    interpolation="none",
    cmap=ListedColormap(list(bodypart_colors.values())),
    vmin=-0.5,
    vmax=len(bodypart_colors) - 0.5,
)
ax_window.imshow(
    window_values[:, None],
    aspect="auto",
    interpolation="none",
    cmap=ListedColormap(window_colors),
    vmin=-0.5,
    vmax=len(window_order) - 0.5,
)

for axis in [ax1, ax_window]:
    axis.set_xticks([])
    axis.set_yticks([])
    axis.set_xlim(-0.5, 0.5)
    axis.set_ylim(len(feature_summary) - 0.5, -0.5)

ax1.set_title("Body\npart", fontsize=8)
ax_window.set_title("Time\nwindow", fontsize=8)

for ytick, row in enumerate(feature_summary.itertuples(index=False)):
    color = {"movement": "red", "geometry": "blue"}.get(row.group, "grey")
    ax2.scatter(
        row.cumulative_importance,
        ytick,
        edgecolors=color,
        facecolors="w",
        alpha=0.5,
        s=30,
        clip_on=False,
    )

ax2.set_xlabel("Cumulative SHAP importance")

bodypart_legend = [
    Patch(color=color, label=label)
    for label, color in bodypart_colors.items()
]
bodypart_legend_artist = ax2.legend(
    handles=bodypart_legend,
    loc="upper right",
    bbox_to_anchor=(1, 1),
    frameon=False,
    fontsize=8,
)

window_legend = [
    Patch(color=color, label=f"{window} s" if window != "none" else "none")
    for window, color in zip(window_order, window_colors)
]
ax2.legend(
    handles=window_legend,
    loc="lower left",
    bbox_to_anchor=(0, 0),
    frameon=False,
    fontsize=7,
)
ax2.add_artist(bodypart_legend_artist)

ax1.set_ylabel("Features (in order of importance)")

sns.despine(ax=ax1, left=True, bottom=True)
sns.despine(ax=ax_window, left=True, bottom=True)
sns.despine(ax=ax2, offset=5)
for axis in [ax1, ax_window, ax2]:
    axis.set_yticks([])
    
save_figure_atomic(f, "figSx_cumul_shap_importance", FIGSFOLDER)

# %%
# Bar + scatter plot comparing signed SHAP importance for movement vs geometry.
# Reverse-coded features have their sign flipped back (aligned_signed_importance)
# so they combine meaningfully with raw features on the same bar; they're
# plotted as squares, raw features as circles. Sign: negative means the
# (non-reverse-coded) feature must go DOWN to drive Appetitive up.

group_order = ["tortuosity", "geometry", "movement"]
group_colors = {"geometry": "#1b9e77", "movement": "#d95f02", "tortuosity": "#7570b3"}
group_positions = {name: i for i, name in enumerate(group_order)}

plot_df = feature_summary[feature_summary["group"].isin(group_order)]
group_means = plot_df.groupby("group")["aligned_signed_importance"].mean().reindex(group_order)

rng = np.random.default_rng(42)

f, ax = plt.subplots(figsize=(4, 3),
                     gridspec_kw={"left": 0.3, "bottom": 0.2})

for group_name in group_order:
    y = group_positions[group_name]
    ax.barh(
        y,
        group_means[group_name],
        color=group_colors[group_name],
        alpha=0.4,
        height=0.6,
        zorder=1,
    )

    group_data = plot_df[plot_df["group"] == group_name]
    jitter = rng.uniform(-0.15, 0.15, size=len(group_data))
    for marker, is_reverse in [("s", True), ("o", False)]:
        mask = group_data["reverse_coded"] == is_reverse
        ax.scatter(
            group_data.loc[mask, "aligned_signed_importance"],
            y + jitter[mask.to_numpy()],
            marker=marker,
            color=group_colors[group_name],
            edgecolors="k",
            linewidths=0.3,
            s=25,
            alpha=0.7,
            zorder=2,
        )

ax.axvline(0, color="0.3", linewidth=0.8, zorder=0)
ax.set_yticks(list(group_positions.values()))
ax.set_yticklabels(["Tortuosity", "Geometry", "Movement"])
ax.set_xlabel("Signed SHAP importance")
ax.set_ylabel("Feature group")

sns.despine(ax=ax, offset=5)

save_figure_atomic(f, "figure_Sx_signed_shap_importance_by_group", FIGSFOLDER)

# %%
# Movement features where signed_importance is positive, i.e. the raw value
# must go UP to drive the Appetitive prediction up. Grouped by reverse_coded
# first so raw vs already-flipped features are easy to tell apart.
positive_movement = (
    feature_summary[
        (feature_summary["group"] == "movement")
        & (feature_summary["signed_importance"] > 0)
    ]
    .sort_values(["reverse_coded", "signed_importance"], ascending=[False, False])
)

print(f"{len(positive_movement)} movement features with positive signed_importance:")
positive_movement[["feature", "importance", "raw_shap_corr", "signed_importance", "reverse_coded"]]


# %%
COLORS

# %%
feature_groups

# %%
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# Build a full feature matrix from the SHAP values table.
feature_matrix = df_shap.drop(
    columns=["Appetitive", "Unnamed: 0", "Prediction_probability", "Sum", "Expected_value"],
    errors="ignore",
)

# Standardise features before PCA.
X = StandardScaler().fit_transform(feature_matrix)

# Fit a full PCA to compute the explained-variance curve.
pca_full = PCA(n_components=None, random_state=42)
pca_full.fit(X)

explained_ratio = pca_full.explained_variance_ratio_
cumulative_ratio = explained_ratio.cumsum()
n_components_95 = int((cumulative_ratio >= 0.95).sum())

print("Explained variance ratios (first 10 components):")
print(explained_ratio[:10])
print("\nCumulative explained variance (first 10 components):")
print(cumulative_ratio[:10])
print(f"\nNumber of components needed for 95% variance: {n_components_95}")

# Also fit a 2-component PCA for a quick scatter plot.
pca_2 = PCA(n_components=2, random_state=42)
X_pca = pca_2.fit_transform(X)

pca_df = pd.DataFrame(
    X_pca,
    columns=["PC1", "PC2"],
    index=feature_matrix.index,
)
pca_df["Appetitive"] = df_shap["Appetitive"].values

ax = pca_df.plot.scatter(x="PC1", y="PC2", c="Appetitive", cmap="viridis", figsize=(7, 5), s=40)
ax.set_title("PCA of SHAP feature matrix")
ax.set_xlabel("PC1")
ax.set_ylabel("PC2")
plt.tight_layout()
plt.show()

pca_df.head()

# Feature loadings for PC1/PC2, coloured by semantic group so we can see
# whether movement and geometry features load onto different components.
# Marker size reflects |signed SHAP importance|, so the most influential
# features (in either direction) stand out within each group.
loading_df = pd.DataFrame(
    pca_2.components_.T,
    index=feature_matrix.columns,
    columns=["PC1", "PC2"],
)
loading_df = loading_df.join(
    feature_summary.set_index("feature")[["group", "signed_importance"]]
)

loading_group_colors = {
    "movement": "#d95f02",
    "geometry": "#1b9e77",
    "tortuosity": "#7570b3",
    "other": "#bdbdbd",
}
max_abs_importance = loading_df["signed_importance"].abs().max()

f, ax = plt.subplots(figsize=(5, 4),
                     gridspec_kw={"left": 0.3, "bottom": 0.2})
for group_name, group_data in loading_df.groupby("group", sort=False):
    sizes = 200 * group_data["signed_importance"].abs() / max_abs_importance
    ax.scatter(
        group_data["PC1"],
        group_data["PC2"],
        s=sizes.clip(lower=10),
        color=loading_group_colors.get(group_name, "#bdbdbd"),
        alpha=0.6,
        edgecolors="k",
        linewidths=0.3,
        label=group_name.title(),
    )

ax.axhline(0, color="0.7", linewidth=0.8)
ax.axvline(0, color="0.7", linewidth=0.8)
# ax.set_title("PCA feature loadings by group\n(marker size = |signed SHAP importance|)")
ax.set_xlabel("PC1 loading")
ax.set_ylabel("PC2 loading")
ax.legend(frameon=False, loc="lower left")
sns.despine(ax=ax, offset=5)

save_figure_atomic(f, "figSx_pca_feature_loadings", FIGSFOLDER)

# %%
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# Build a full feature matrix from the raw SHAP values table.
# Keep the Appetitive label as a separate column for colouring later.
feature_matrix = df_shap.drop(columns=["Appetitive"], errors="ignore")
# feature_matrix = df_shap

# Standardise features before PCA, since they may be on different scales.
X = StandardScaler().fit_transform(feature_matrix)

# Fit PCA and reduce to 2 components.
pca = PCA(n_components=2, random_state=42)
X_pca = pca.fit_transform(X)

pca_df = pd.DataFrame(
    X_pca,
    columns=["PC1", "PC2"],
    index=feature_matrix.index,
)
pca_df["Appetitive"] = df_shap["Appetitive"].values

print("Explained variance ratio:", pca.explained_variance_ratio_)
print("Cumulative explained variance:", pca.explained_variance_ratio_.sum())

# Scatter plot of the first two PCs
ax = pca_df.plot.scatter(x="PC1", y="PC2", c="Appetitive", cmap="viridis", figsize=(7, 5), s=40)
ax.set_title("PCA of SHAP feature matrix")
ax.set_xlabel("PC1")
ax.set_ylabel("PC2")
plt.tight_layout()
plt.show()

pca_df.head()

# %%
# Inspect and plot feature loadings for the first two principal components.
# Each point is a feature; colour indicates its semantic feature group.

loading_df = pd.DataFrame(
    pca_2.components_.T,
    index=feature_matrix.columns,
    columns=["PC1", "PC2"],
)

loading_df["abs_PC1"] = loading_df["PC1"].abs()
loading_df["abs_PC2"] = loading_df["PC2"].abs()
loading_df["group"] = loading_df.index.to_series().apply(assign_group)

print("Top features for PC1:")
print(loading_df.sort_values("abs_PC1", ascending=False).head(15)[["PC1"]])
print("\nTop features for PC2:")
print(loading_df.sort_values("abs_PC2", ascending=False).head(15)[["PC2"]])

feature_group_colors = {
    "movement": "#d95f02",
    "geometry": "#1b9e77",
    "other": "#7570b3",
}

fig, ax = plt.subplots(figsize=(8, 6))
for group_name, group_data in loading_df.groupby("group", sort=False):
    ax.scatter(
        group_data["PC1"],
        group_data["PC2"],
        label=group_name.title(),
        color=feature_group_colors.get(group_name, "#bdbdbd"),
        alpha=0.75,
        s=35,
    )

ax.axhline(0, color="0.7", linewidth=0.8)
ax.axvline(0, color="0.7", linewidth=0.8)
ax.set_title("PCA feature loadings coloured by feature group")
ax.set_xlabel("PC1 loading")
ax.set_ylabel("PC2 loading")
ax.legend(frameon=False)
plt.tight_layout()
plt.show()

loading_df.sort_values("abs_PC1", ascending=False).head(10)[["PC1", "PC2", "group"]]

# %%
# Use df_shap_raw to inspect how SHAP values change across the raw feature distribution.
# Instead of splitting by positive vs negative SHAP, we split the raw feature into quartiles
# and compare the average SHAP value in each quartile.
#
# This is often more informative for questions like: does a larger raw value correspond to
# a larger positive contribution, or does the relationship saturate or reverse?

# Restrict to the class of interest.
subset_mask = df_shap["Appetitive"] == 1
shap_subset = df_shap.loc[subset_mask].copy()
raw_subset = df_shap_raw.loc[subset_mask].copy()

# Keep only the feature columns from the raw table.
raw_features = [c for c in raw_subset.columns if c not in ["Appetitive", "Unnamed: 0", "Prediction_probability", "Sum", "Expected_value"]]

# Build a per-feature summary across quartiles of the raw feature values.
quartile_summary = []
for feat in raw_features:
    if feat not in shap_subset.columns:
        continue
    raw_vals = raw_subset[feat]
    shap_vals = shap_subset[feat]
    if raw_vals.isna().all() or shap_vals.isna().all():
        continue
    try:
        q = pd.qcut(raw_vals, q=4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
    except ValueError:
        q = pd.qcut(raw_vals, q=3, labels=["Q1", "Q2", "Q3"], duplicates="drop")
    if q.nunique() < 2:
        continue
    grouped = shap_vals.groupby(q).mean().rename("mean_shap")
    ordered_labels = [lab for lab in ["Q1", "Q2", "Q3", "Q4"] if lab in grouped.index]
    grouped = grouped.reindex(ordered_labels)
    quartile_summary.append({
        "feature": feat,
        "q1_mean_shap": grouped.get("Q1", float("nan")),
        "q2_mean_shap": grouped.get("Q2", float("nan")),
        "q3_mean_shap": grouped.get("Q3", float("nan")),
        "q4_mean_shap": grouped.get("Q4", float("nan")),
        "trend": grouped.iloc[-1] - grouped.iloc[0] if len(grouped.dropna()) >= 2 else float("nan"),
    })

quartile_summary = pd.DataFrame(quartile_summary)

# Rank by the overall increase from Q1 to Q4.
quartile_summary = quartile_summary.sort_values("trend", ascending=False)

print("Top features where SHAP increases across the raw-value quartiles within Appetitive == 1:")
print(quartile_summary.head(20)[["feature", "q1_mean_shap", "q2_mean_shap", "q3_mean_shap", "q4_mean_shap", "trend"]])

print("\nTop features where SHAP decreases across the raw-value quartiles within Appetitive == 1:")
print(quartile_summary.sort_values("trend", ascending=True).head(20)[["feature", "q1_mean_shap", "q2_mean_shap", "q3_mean_shap", "q4_mean_shap", "trend"]])

# Plot one feature with the strongest increasing trend and one with the strongest decreasing trend.
if not quartile_summary.empty:
    top_increasing = quartile_summary.sort_values("trend", ascending=False).iloc[0]
    top_decreasing = quartile_summary.sort_values("trend", ascending=True).iloc[0]
    plot_df = pd.DataFrame({
        "quartile": ["Q1", "Q2", "Q3", "Q4"],
        top_increasing["feature"]: [top_increasing["q1_mean_shap"], top_increasing["q2_mean_shap"], top_increasing["q3_mean_shap"], top_increasing["q4_mean_shap"]],
        top_decreasing["feature"]: [top_decreasing["q1_mean_shap"], top_decreasing["q2_mean_shap"], top_decreasing["q3_mean_shap"], top_decreasing["q4_mean_shap"]],
    })
    ax = plot_df.set_index("quartile").plot(marker="o", figsize=(7, 4))
    ax.set_title("Example SHAP-vs-raw quartile trends")
    ax.set_ylabel("Mean SHAP value")
    ax.set_xlabel("Raw feature quartile")
    plt.tight_layout()
    plt.show()

quartile_summary.head(20)

# %% [markdown]
# Interpretation of the SHAP/PCA results
#
# The main conclusion from the PCA is that the structure is dominated by two broad patterns:
#
# 1. A geometry/spread pattern (PC1), driven by hull-distance and body-shape features.
# 2. A movement pattern (PC2), driven by tail-base, head-base, nose and rolling-window movement features.
#
# This means that movement is the more prominent axis in the feature space, while body geometry contributes a separate but secondary axis. The relative-to-average or reverse-coded features do not form a clearly separate dominant PCA axis on their own; instead, they are mixed into the broader structure rather than acting as a standalone component.
#
# This interpretation is consistent with the SHAP results. The strongest SHAP effects are concentrated in movement-related features, and the direction of effect is determined by the feature construction:
#
# - raw, sustained movement features tend to be associated with lower appetitive scores,
# - reverse-coded relative-to-average features tend to be associated with higher appetitive scores when their value becomes larger.
#
# The same broad structure is visible when the PCA is repeated on the full dataset (all frames, not only Appetitive == 1). The main distinction remains geometry/spread versus movement, with movement remaining the more dominant axis.

# %%
# PCA on the full dataset (all frames, not just Appetitive == 1)
# This repeats the PCA analysis using all rows so we can check whether the same structure appears outside the appetitive subset.

feature_matrix_all = df_shap.drop(
    columns=["Appetitive", "Unnamed: 0", "Prediction_probability", "Sum", "Expected_value"],
    errors="ignore",
)

# Standardise the full feature matrix.
X_all = StandardScaler().fit_transform(feature_matrix_all)

# Fit PCA with two components for comparison.
pca_all = PCA(n_components=2, random_state=42)
X_pca_all = pca_all.fit_transform(X_all)

pca_df_all = pd.DataFrame(
    X_pca_all,
    columns=["PC1", "PC2"],
    index=feature_matrix_all.index,
)
pca_df_all["Appetitive"] = df_shap["Appetitive"].values

# Inspect the loadings for the first two components.
loading_df_all = pd.DataFrame(
    pca_all.components_.T,
    index=feature_matrix_all.columns,
    columns=["PC1", "PC2"],
)
loading_df_all["abs_PC1"] = loading_df_all["PC1"].abs()
loading_df_all["abs_PC2"] = loading_df_all["PC2"].abs()

print("Explained variance ratio for full-data PCA:")
print(pca_all.explained_variance_ratio_)
print("\nTop loadings for PC1 (full data):")
print(loading_df_all.sort_values("abs_PC1", ascending=False).head(15)[["PC1"]])
print("\nTop loadings for PC2 (full data):")
print(loading_df_all.sort_values("abs_PC2", ascending=False).head(15)[["PC2"]])

# Scatter plot of the first two PCs coloured by Appetitive label.
ax = pca_df_all.plot.scatter(x="PC1", y="PC2", c="Appetitive", cmap="viridis", figsize=(7, 5), s=40)
ax.set_title("PCA of SHAP feature matrix (all frames)")
ax.set_xlabel("PC1")
ax.set_ylabel("PC2")
plt.tight_layout()
plt.show()

pca_df_all.head()

# %%
# Compare rolling-window suffix importance across feature groups
# We parse the suffix in each feature name (e.g. 2, 5, 6, 7.5, 10) and summarize
# mean absolute SHAP importance for each feature group.

subset = df_shap.query("Appetitive == 1").copy()

# Build a long-form summary of feature importance.
summary_long = pd.DataFrame({
    "feature": summary.index,
    "importance": summary.abs().values,
})
summary_long["group"] = summary_long["feature"].apply(assign_group)

# Extract the rolling-window suffix from feature names.
def get_window_suffix(feature: str):
    m = re.search(r"_(\d+(?:\.\d+)?)$", feature)
    if m:
        return float(m.group(1))
    return None

summary_long["window"] = summary_long["feature"].apply(get_window_suffix)

# Keep only features with a rolling-window suffix.
rolling_summary = summary_long.dropna(subset=["window"]).copy()

# Show the distribution of importance by suffix for each group.
window_group_summary = (
    rolling_summary.groupby(["group", "window"], as_index=False)["importance"]
    .mean()
    .sort_values(["group", "window"])
)

print("Mean absolute SHAP importance by feature group and rolling-window suffix:")
print(window_group_summary.to_string(index=False))

# Also show the top features for each window size.
for window in sorted(rolling_summary["window"].unique()):
    top_for_window = (
        rolling_summary[rolling_summary["window"] == window]
        .sort_values("importance", ascending=False)
        .head(10)
    )
    print(f"\nWindow {window} top features:")
    print(top_for_window[["group", "feature", "importance"]].to_string(index=False))

# A simple plot: one line per group, with window size on x-axis.
plt.figure(figsize=(8, 4))
for group_name, grp in window_group_summary.groupby("group"):
    plt.plot(grp["window"], grp["importance"], marker="o", label=group_name)
plt.xlabel("Rolling-window suffix")
plt.ylabel("Mean |SHAP| importance")
plt.title("Importance of rolling-window suffixes by feature group")
plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left")
plt.tight_layout()
plt.show()
