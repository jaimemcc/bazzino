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
from pathlib import Path
import re
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
import seaborn as sns

from figure_config import (
    configure_matplotlib, COLORS, HEATMAP_CMAP_DIV,
    # HEATMAP_CMAP_RED, HEATMAP_CMAP_BLUE,
    # DATAFOLDER, RESULTSFOLDER, FIGSFOLDER,
    # HEATMAP_VLIM_BEHAV, YLIMS_BEHAV,
    # BEHAV_SMOOTH_WINDOW, SAVE_FIGS
)



# %%
DATAFOLDER = Path("../data/shap_values")

df_shap, df_shap_raw = [], []

for folder in [f for f in DATAFOLDER.iterdir() if f.is_dir()]:
    print(f"Folder: {folder.name}")
    tmp = pd.read_csv(folder / "SHAP_values_Appetitive.csv")
    df_shap.append(tmp)
    tmp = pd.read_csv(folder / "RAW_SHAP_feature_values_Appetitive.csv")
    df_shap_raw.append(tmp)


df_shap = pd.concat(df_shap, axis=0).reset_index(drop=True)
df_shap_raw = pd.concat(df_shap_raw, axis=0).reset_index(drop=True)

df_shap = (
    df_shap
    .drop(columns=["Unnamed: 0", "Prediction_probability", "Sum", "Expected_value"], errors="ignore")
)

df_shap_raw = (
    df_shap_raw
    .drop(columns=["Unnamed: 0"], errors="ignore")
)




# %%
df_shap.shape


# %%
df_shap.query("Appetitive == 1").mean().sort_values(ascending=False)

# %%

summary = (
    df_shap.query("Appetitive == 1")
    .drop(columns=["Appetitive", "Unnamed: 0", "Prediction_probability", "Sum", "Expected_value"], errors="ignore")
    .abs()
    .mean()
    .sort_values(ascending=False)
)

top_features = summary.head(15)

ax = top_features.plot.barh(figsize=(8, 6), color="#1f77b4")
ax.invert_yaxis()
ax.set_title("Top 15 features by mean absolute SHAP value")
ax.set_xlabel("Mean |SHAP value|")
ax.set_ylabel("Feature")
plt.tight_layout()
plt.show()

summary.head(15)


# %%
# Group features into movement, body/hull geometry, or other.

def assign_group(col: str) -> str:
    if re.search(
        r"^(Movement_mouse_(nose|tail_base|left_ear|right_ear|head_base)|"
        r"Total_movement_all_bodyparts_M1|Total_movement_M1_|"
        r"Tail_base_movement_M1_|Head_base_movement_M1_|Nose_movement_M1_)"
        r"|(_deviation|_percentile_rank)$",
        col,
    ):
        return "movement"
    if re.search(
        r"^(Mouse_nose_to_tail|Mouse_head_to_tail|Mouse_Ear_distance|"
        r"M1_|Mouse1_(smallest|largest|mean)_euclid_distances_)",
        col,
    ):
        return "body_hull_geometry"
    return "other"

feature_groups = pd.DataFrame({
    "feature": summary.index,
    "importance": summary.abs().values,
})
feature_groups["group"] = feature_groups["feature"].apply(assign_group)

# %%
cumul_summary


# %%
# Add bodypart groups and timewindow
 
def assign_bodypart(feature: str) -> str:
    feature_lower = feature.lower()
    if "all_bodyparts" in feature_lower:
        return "all body parts"
    for key, label in [
        ("nose", "nose"),
        ("tail_base", "tail base"),
        ("head_base", "head base"),
        ("left_ear", "left ear"),
        ("ear_left", "left ear"),
        ("right_ear", "right ear"),
        ("ear_right", "right ear"),
    ]:
        if key in feature_lower:
            return label
    if feature_lower.startswith(("mouse_", "mouse1_")) or "total" in feature_lower:
        return "whole mouse"
    print(f"Feature: {feature}, no body part assigned")
    return "other"

def get_time_window(feature: str):
    match = re.search(r"_(\d+(?:\.\d+)?)$", feature)
    return match.group(1) if match else "none"

cumul_summary = (
    feature_groups
    .assign(cumulative_importance=lambda x: x.importance.cumsum() / np.sum(x.importance))
)
cumul_summary["bodypart"] = cumul_summary["feature"].apply(assign_bodypart)
cumul_summary["time_window"] = cumul_summary["feature"].apply(get_time_window)

# %%
COLORS

# %%
# Plot all features with cumulative SHAP importance on x and feature rows on y.

bodypart_colors = {
    "nose": "#d95f02",
    "tail base": "#1b9e77",
    "head base": "#7570b3",
    "left ear": "#e7298a",
    "right ear": "#66a61e",
    "all body parts": "#e6ab02",
    "whole mouse": "#a6761d",
    "other": "#bdbdbd",
}
bodypart_codes = {name: code for code, name in enumerate(bodypart_colors)}
bodypart_values = cumul_summary["bodypart"].map(bodypart_codes).to_numpy()

window_order = ["none", "2", "5", "6", "7.5", "10"]
window_colors = ["#d9d9d9", "#440154", "#31688e", "#35b779", "#90d743", "#fde725"]
window_codes = {name: code for code, name in enumerate(window_order)}
window_values = cumul_summary["time_window"].map(window_codes).to_numpy()

f, [ax1, ax_window, ax2] = plt.subplots(
    figsize=(5.5, 6),
    ncols=3,
    sharey=True,
    gridspec_kw={"width_ratios": [0.05, 0.05, 1], "wspace": 0.08},
)

# One-column heatmaps: one colored cell for each feature row.
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
    axis.set_ylim(len(cumul_summary) - 0.5, -0.5)

ax1.set_title("Body\npart", fontsize=8)
ax_window.set_title("Time\nwindow", fontsize=8)

for ytick, row in enumerate(cumul_summary.itertuples(index=False)):
    group = str(row.group)
    if "movement" in group:
        color = "red"
    elif "body_hull_geometry" in group:
        color = "blue"
    else:
        color = "grey"
    ax2.scatter(
        row.cumulative_importance,
        ytick,
        edgecolors=color,
        facecolors="w",
        alpha=0.5,
        s=30,
        clip_on=False
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

ax1.set_ylabel("Features")

sns.despine(ax=ax1, left=True, bottom=True)
sns.despine(ax=ax_window, left=True, bottom=True)
sns.despine(ax=ax2, offset=5)
for axis in [ax1, ax_window, ax2]:
    axis.set_yticks([])

# %%
# make fig showing aggregated SHAP importance for movement

# need to get info from RAW file and place by features in df

f, ax = plt.subplots()



# %%
# maybe make figure where an erro is shown representing deviation/variability in each feature's importance


# %%
# Sum importance within each semantic family
summary_by_group = (
    feature_groups.groupby("group", as_index=False)["importance"]
    .sum()
    .sort_values("importance", ascending=False)
)

summary_by_group["group"] = summary_by_group["group"].replace({
    "movement": "Movement",
    "body_hull_geometry": "Body/hull geometry",
    "other": "Other",
})

summary_by_group

ax = summary_by_group.set_index("group").plot.barh(figsize=(8, 4), color="#2ca02c")
ax.invert_yaxis()
ax.set_title("SHAP importance by semantic feature group")
ax.set_xlabel("Total absolute SHAP importance")
ax.set_ylabel("Feature group")
plt.tight_layout()
plt.show()

# %%
feature_groups

# %%
# Create a readable table showing which features were grouped together
feature_group_table = (
    feature_groups
    .assign(group_name=feature_groups["group"].replace({
        "movement": "Movement",
        "body_hull_geometry": "Body/hull geometry",
        "other": "Other",
    }))
    .sort_values(["group_name", "feature"])
    [["group_name", "feature", "importance"]]
)

feature_group_table.head(50)

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
# Inspect the feature loadings for the first two principal components.
# Larger absolute values indicate stronger contribution to that component.

loading_df = pd.DataFrame(
    pca_2.components_.T,
    index=feature_matrix.columns,
    columns=["PC1", "PC2"],
)

loading_df["abs_PC1"] = loading_df["PC1"].abs()
loading_df["abs_PC2"] = loading_df["PC2"].abs()

print("Top features for PC1:")
print(loading_df.sort_values("abs_PC1", ascending=False).head(15)[["PC1"]])
print("\nTop features for PC2:")
print(loading_df.sort_values("abs_PC2", ascending=False).head(15)[["PC2"]])

# Optional: show a compact heatmap-like view of the top contributors
loading_df.sort_values("abs_PC1", ascending=False).head(10)[["PC1", "PC2"]]

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
