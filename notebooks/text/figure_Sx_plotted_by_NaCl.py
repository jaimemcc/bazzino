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
# ## Plotting behaviour and dopamine as a function of NaCl consumed

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

from figure_plotting import (make_correlation_plot_simba
)

from model_fit_helpers import (
    fit_sigmoid,
    _sigmoid_model,
)

# from utils import make_realigned_trials
from realignment_helpers import get_realigned_data

SAVE_FIGS = True
FIT_TO_RAW_DATA = True

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

print(f"Loaded assembled data from {assembled_data_path}")
print(f"\nData structure:")
print(f"  - x_array shape: {x_array.shape}")

# %%
x_array.columns


# %%
def make_auc_data_per_trial(x_array, condition, infusiontype, y_col):
    return (
        x_array
        .query("condition == @condition & infusiontype == @infusiontype")
        .groupby("trial")[y_col]
        .mean()
        .values
    )
    

deplete_10 = make_auc_data_per_trial(x_array, "deplete", "10NaCl", "simba_median_balance")
deplete_45 = make_auc_data_per_trial(x_array, "deplete", "45NaCl", "simba_median_balance")

# %%
f = make_correlation_plot_simba(deplete_10, deplete_45, colors[0], colors[1], yaxis=True)

# %%
from scipy.stats import linregress

scaling_factor_100mM = 2.3376
scaling_factor_450mM = 2.3376 * 4.5

f, ax = plt.subplots(figsize=(1.8, 1.8))

x_vals_100 = np.arange(0,49)*scaling_factor_100mM
# ax.scatter(x_vals_100, deplete_10, color=colors[0], label="10 mM NaCl", alpha=0.3)

ax.plot(x_vals_100, deplete_10, color=colors[2], alpha=0.7)

# a, b, r_value, p_value, std_err = linregress(x_vals_100, deplete_10)
# ax.plot(x_vals_100, a*x_vals_100 + b, color=colors[0], alpha=1)

x_vals_450 = np.arange(0,49)*scaling_factor_450mM
# ax.scatter(x_vals_450, deplete_45, color=colors[1], label="45 mM NaCl", alpha=0.3)
ax.plot(x_vals_450, deplete_45, color=colors[3], alpha=0.7)

# a, b, r_value, p_value, std_err = linregress(x_vals_450, deplete_45)
# ax.plot(x_vals_450, a*x_vals_450 + b, color=colors[1], alpha=1)

ax.set_xlabel("NaCl consumed (mg)")
ax.set_ylabel("Appetitive behaviour")

sns.despine(ax=ax, offset=5)

# %%
deplete_10_behav = make_auc_data_per_trial(x_array, "deplete", "10NaCl", "simba_median_balance")
deplete_45_behav = make_auc_data_per_trial(x_array, "deplete", "45NaCl", "simba_median_balance")

deplete_10_da = make_auc_data_per_trial(x_array, "deplete", "10NaCl", "auc_snips")
deplete_45_da = make_auc_data_per_trial(x_array, "deplete", "45NaCl", "auc_snips")

# %%
from scipy.stats import linregress

scaling_factor_100mM = 2.3376
scaling_factor_450mM = 2.3376 * 4.5

f, ax = plt.subplots(figsize=(1.8, 1.8))

x_vals_100 = np.arange(0,49)*scaling_factor_100mM
# ax.scatter(x_vals_100, deplete_10, color=colors[2], label="10 mM NaCl", alpha=0.3)

ax.plot(x_vals_100, deplete_10, color=colors[2], alpha=0.7)

# a, b, r_value, p_value, std_err = linregress(x_vals_100, deplete_10)
# ax.plot(x_vals_100, a*x_vals_100 + b, color=colors[0], alpha=1)

x_vals_450 = np.arange(0,49)*scaling_factor_450mM
# ax.scatter(x_vals_450, deplete_45, color=colors[3], label="45 mM NaCl", alpha=0.3)
ax.plot(x_vals_450, deplete_45, color=colors[3], alpha=0.7)

# a, b, r_value, p_value, std_err = linregress(x_vals_450, deplete_45)
# ax.plot(x_vals_450, a*x_vals_450 + b, color=colors[1], alpha=1)

ax.set_xlabel("NaCl consumed (mg)")
ax.set_ylabel("Dopamine (AUC)")

sns.despine(ax=ax, offset=5)

# %%
f, ax = plt.subplots(ncols=2, figsize=(3, 1.8), sharey=True, sharex=True,
                     gridspec_kw={"wspace": 0.5})
ax[0].scatter(x_vals_100, deplete_10, color=colors[2], label="10 mM NaCl", alpha=0.3)
ax[1].scatter(x_vals_450, deplete_45, color=colors[3], label="45 mM NaCl", alpha=0.3)

for axis in ax:
    axis.set_xlabel("NaCl consumed (mg)")
    sns.despine(ax=axis, offset=5)
    
ax[0].set_ylabel("Dopamine (AUC)")


# %%
def make_triple_plot(x_vals_100, deplete_10, x_vals_450, deplete_45, colors=colors):
    f, ax = plt.subplots(ncols=3, figsize=(5, 1.8), sharey=True,
                         gridspec_kw={"wspace": 0.5, "width_ratios": [0.2, 0.2, 0.6]})
    
    ax[0].scatter(x_vals_100, deplete_10, color=colors[2], label="10 mM NaCl", alpha=0.3, clip_on=False)
    ax[1].scatter(x_vals_450, deplete_45, color=colors[3], label="45 mM NaCl", alpha=0.3, clip_on=False)
    
    ax[1].axvspan(np.min(x_vals_100), np.max(x_vals_100), color="k", alpha=0.05, zorder=-10)
    
    # ax[0].set_xlim(ax[1].get_xlim())

    ax[2].plot(x_vals_100, deplete_10, color=colors[2], alpha=0.5, marker="o")
    
    x_vals_450_red, deplete_45_red = zip(*[(x, y) for x, y in zip(x_vals_450, deplete_45) if x <= np.max(x_vals_100)])
    ax[2].plot(x_vals_450_red, deplete_45_red, color=colors[3], alpha=0.7, marker="o")
    
    for axis in ax:
        axis.set_xlabel("NaCl (mg)")
        sns.despine(ax=axis, offset=5)

    return f, ax

deplete_10_behav = make_auc_data_per_trial(x_array, "deplete", "10NaCl", "simba_median_balance")
deplete_45_behav = make_auc_data_per_trial(x_array, "deplete", "45NaCl", "simba_median_balance")

deplete_10_da = make_auc_data_per_trial(x_array, "deplete", "10NaCl", "auc_snips")
deplete_45_da = make_auc_data_per_trial(x_array, "deplete", "45NaCl", "auc_snips")

f, ax = make_triple_plot(x_vals_100, deplete_10_behav, x_vals_450, deplete_45_behav, colors=colors)
ax[0].set_ylabel("Appetitive behaviour")

f, ax = make_triple_plot(x_vals_100, deplete_10_da, x_vals_450, deplete_45_da, colors=colors)
ax[0].set_ylabel("Dopamine (AUC)")

# %%
