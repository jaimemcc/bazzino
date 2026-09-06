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

# %%
# %load_ext autoreload
# %autoreload 2

from pathlib import Path
import sys

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
assembled_data_path = DATAFOLDER / "assembled_data.pickle"

with open(assembled_data_path, "rb") as f:
    data = dill.load(f)

# Extract main components
x_array = data["x_array"].copy()

# reverses cluster assignments
x_array["cluster_photo"] = x_array["cluster_photo"].map({0: 1, 1: 0})

print(f"Loaded assembled data from {assembled_data_path}")
print(f"\nData structure:")
print(f"  - x_array shape: {x_array.shape}")

# %%
# additional analysis for revision

da_transition = []
behav_transition = []
for rat in x_array.id.unique():

    rat_data = (x_array
    .query("condition == 'deplete' and infusiontype == '45NaCl' and id == @rat and realigned_trials_da_velocity == 0")
    .rename(columns={"trial": "da_transition"})
    )
    if len(rat_data) > 0:
        da_transition.append(rat_data.loc[:, ["id", "sex", "da_transition"]])

    rat_data = (x_array
    .query("condition == 'deplete' and infusiontype == '45NaCl' and id == @rat and realigned_trials_behav == 0")
    .rename(columns={"trial": "behav_transition"})
    )
    if len(rat_data) > 0:
        behav_transition.append(rat_data.loc[:, ["id", "sex", "behav_transition"]])

da_transition = pd.concat(da_transition).reset_index(drop=True)
behav_transition = pd.concat(behav_transition).reset_index(drop=True)

df_transition = pd.merge(left=da_transition, right=behav_transition)
df_transition["bodyweight"] = [247, 237, 266, 470, 332, 304, 454, 540, 281, 278]

# %%
df_transition


# %%
f, ax = plt.subplots(figsize=(4, 1.5))

for _, row in df_transition.iterrows():
    if row.sex == "M":
        marker = "^"
    elif row.sex == "F":
        marker = "o"

    jitter = np.random.random(1) * 0.1
    print(jitter)

    ax.scatter(row.da_transition, 1+jitter[0],
               marker=marker,
               facecolors="w", edgecolors=DA_COLOR)
    ax.scatter(row.behav_transition, 2+jitter[0],
               marker=marker,
               facecolors="w", edgecolors=BEHAV_COLOR)

    ax.plot((row.da_transition, row.behav_transition),
            (1+jitter[0], 2+jitter[0]),
            color="k", alpha=0.3,
            zorder=-10)


ax.axvline(15, linestyle="--", alpha=0.5, zorder=-20, color="k")

ax.set_yticks([1,2], labels=["DA transition", "Behav transition"])
ax.set_xlabel("Trial")
sns.despine(offset=5)



# %%
# plot with body weight
f, [ax1, ax2] = plt.subplots(figsize=(3, 1.5), ncols=2, sharey=True)

for _, row in df_transition.iterrows():
    if row.sex == "M":
        marker = "^"
    elif row.sex == "F":
        marker = "o"

    ax1.scatter(row.da_transition, row.bodyweight,
               marker=marker,
               facecolors="w", edgecolors=DA_COLOR)

    ax2.scatter(row.behav_transition, row.bodyweight,
            marker=marker,
            facecolors="w", edgecolors=BEHAV_COLOR)

ax1.set_ylabel("Body weight (g)")
ax1.set_xlabel("DA transition")
ax2.set_xlabel("Behav transition")

for axis in [ax1, ax2]:
    sns.despine(ax=axis, offset=5)
    axis.set_xlim(0, 50)




# %%
