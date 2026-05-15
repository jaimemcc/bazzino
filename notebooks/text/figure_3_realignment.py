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
# # Paper Figures: Figure 3 - Data Aligned to Sigmoidal Transitions
#
# This notebook generates Figure 3 showing dopamine (photometry) response aligned to the transition point where animals shift cluster membership. Uses sigmoidal transition points calculated in `src/assemble_all_data.py`.
#
# **Figure 3: Transition-Aligned Dopamine Response** — Photometry AUC heatmaps and time series for each rat, centered at their individual transition points.

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
from scipy.stats import mannwhitneyu

from matplotlib.colors import to_hex

from trompy import save_figure_atomic

from figure_config import (
    configure_matplotlib, COLORS, HEATMAP_CMAP_DIV,
    DATAFOLDER, RESULTSFOLDER, FIGSFOLDER,
    SAVE_FIGS
)

from figure_plotting import (plot_auc_and_sigmoid, plot_realigned_behaviour,
                             make_realignment_panel_behav_and_da,
                             make_realignment_panel_only_one_parameter
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


# %% [markdown]
# ## Load Assembled Data
# Load the complete dataset with transition points and trial-aligned indices.

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
x_array.columns

# %% [markdown]
# ## Realignment figure

# %%
## prepare shared data
from realignment_helpers import only_keep_complete_trials

def prepare_shared_data(df, realignment_col, da_col, simba_col, initial_da_k=-1, initial_behav_k=-1):
    """Compute shared fits/arrays used by all Figure 4B panels."""
    
    ntrials_for_slope = 5
    df_realigned = only_keep_complete_trials(df, realignment_col).copy()
    df_reali_for_slope = df_realigned.query(f"{realignment_col} >= -{ntrials_for_slope} and {realignment_col} <= {ntrials_for_slope}").copy()
    
    ntrials = df_realigned[realignment_col].nunique()
    midway_point = int(df[df[realignment_col] == 0].trial.mean())
    min_trial, max_trial = midway_point - ntrials // 2, midway_point + ntrials // 2
    
    print(midway_point, min_trial, max_trial)
    
    df_original = df.query("trial >= @min_trial and trial <= @max_trial").copy()
    
    min_trial, max_trial = midway_point - ntrials_for_slope, midway_point + ntrials_for_slope
    df_orig_for_slope = df_original.query("trial >= @min_trial and trial <= @max_trial").copy()

    popt_behav_orig = fit_sigmoid(df, "trial", simba_col, fit_to_raw_data=FIT_TO_RAW_DATA, initial_k=initial_behav_k)
    popt_behav_realigned = fit_sigmoid(df_realigned, realignment_col, simba_col, fit_to_raw_data=FIT_TO_RAW_DATA, initial_k=initial_behav_k)
    popt_da_orig = fit_sigmoid(df, "trial", da_col, fit_to_raw_data=FIT_TO_RAW_DATA, initial_k=initial_da_k)
    popt_da_realigned = fit_sigmoid(df_realigned, realignment_col, da_col, fit_to_raw_data=FIT_TO_RAW_DATA, initial_k=initial_da_k)

    k_behav_orig = popt_behav_orig[2] if np.all(np.isfinite(popt_behav_orig)) else np.nan
    k_behav_realigned = popt_behav_realigned[2] if np.all(np.isfinite(popt_behav_realigned)) else np.nan
    k_da_orig = popt_da_orig[2] if np.all(np.isfinite(popt_da_orig)) else np.nan
    k_da_realigned = popt_da_realigned[2] if np.all(np.isfinite(popt_da_realigned)) else np.nan

    return {
        "df_original": df_original,
        "df_orig_for_slope": df_orig_for_slope,
        "df_realigned": df_realigned,
        "df_reali_for_slope": df_reali_for_slope,
        "k_behav_orig": k_behav_orig,
        "k_behav_realigned": k_behav_realigned,
        "k_da_orig": k_da_orig,
        "k_da_realigned": k_da_realigned,
        "popt_behav_orig": popt_behav_orig,
        "popt_behav_realigned": popt_behav_realigned,
        "popt_da_orig": popt_da_orig,
        "popt_da_realigned": popt_da_realigned,
        "ntrials_for_slope": ntrials_for_slope,
    }



# %%
# new figure 3

DA_COLUMN = "auc_snips"
SIMBA_COLUMN = "simba_median_balance"
REALIGNMENT_COLUMN = "realigned_trials_behav"
RATS_TO_EXCLUDE = ["PB27", "PB71"]  # These rats have poor fits of behavioral data

subset_aligned = get_realigned_data(x_array, REALIGNMENT_COLUMN, rats_to_exclude=RATS_TO_EXCLUDE, verbose=False)
data_for_figure = prepare_shared_data(subset_aligned, REALIGNMENT_COLUMN, DA_COLUMN, SIMBA_COLUMN)

fig, ax = make_realignment_panel_only_one_parameter(data_for_figure, REALIGNMENT_COLUMN, SIMBA_COLUMN, "behav",
                                          color=BEHAV_COLOR, do_not_plot_sigmoid=True)
save_figure_atomic(fig, FIGSFOLDER / "fig3_realignment_behaviour")

fig, ax = make_realignment_panel_only_one_parameter(data_for_figure, REALIGNMENT_COLUMN, DA_COLUMN, "da",
                                          color=DA_COLOR, do_not_plot_sigmoid=True)
save_figure_atomic(fig, FIGSFOLDER / "fig3_realignment_da")


# %%
from figure_plotting import make_realignment_slope_panel

f, ax, slopes =make_realignment_slope_panel(data_for_figure, SIMBA_COLUMN, REALIGNMENT_COLUMN,
                             color=BEHAV_COLOR)

ax.set_yticks([-0.5, 0, 0.5, 1])
ax.set_ylabel("App. behaviour")

save_figure_atomic(f, FIGSFOLDER / "fig3_realignment_slope_behav")


print("Behaviour slopes")
print("Orig. slope = {: .3f}".format(slopes["trial"]["slope"]))
print("Reali. slope = {: .3f}".format(slopes["realigned_trials_behav"]["slope"]))

f, ax, slopes =make_realignment_slope_panel(data_for_figure, DA_COLUMN, REALIGNMENT_COLUMN,
                             color=DA_COLOR)


# ax.set_yticks([-0.5, 0, 0.5, 1])
ax.set_ylabel("Dopamine AUC")
# ax.text(0.5, 0.9, "Orig. slope = {: .3f}".format(slopes["trial"]["slope"]), transform=ax.transAxes, ha="left", va="bottom", fontsize=10)
# ax.text(0.5, 0.8, "Reali. slope = {: .3f}".format(slopes["realigned_trials_behav"]["slope"]), transform=ax.transAxes, ha="left", va="bottom", fontsize=10)
save_figure_atomic(f, FIGSFOLDER / "fig3_realignment_slope_da")

print("Dopamine slopes")
print("Orig. slope = {: .3f}".format(slopes["trial"]["slope"]))
print("Reali. slope = {: .3f}".format(slopes["realigned_trials_behav"]["slope"]))

# %%
from bootstrap_and_shuffle_helpers import get_bootstrapped_distribution, get_bootstrapped_distribution_using_slopes

DA_COLUMN = "auc_snips"
SIMBA_COLUMN = "simba_median_balance"
REALIGNMENT_COLUMN = "realigned_trials_behav"
N_BOOTSTRAPS = 100

df_orig_for_slope = data_for_figure["df_orig_for_slope"]
df_reali_for_slope = data_for_figure["df_reali_for_slope"]

bootstrapped_behav_orig = get_bootstrapped_distribution_using_slopes(df_orig_for_slope, "trial", SIMBA_COLUMN, n_bootstraps=N_BOOTSTRAPS)
bootstrapped_behav_realigned = get_bootstrapped_distribution_using_slopes(df_reali_for_slope, REALIGNMENT_COLUMN, SIMBA_COLUMN, n_bootstraps=N_BOOTSTRAPS)

bootstrapped_da_orig = get_bootstrapped_distribution_using_slopes(df_orig_for_slope, "trial", DA_COLUMN, n_bootstraps=N_BOOTSTRAPS)
bootstrapped_da_realigned = get_bootstrapped_distribution_using_slopes(df_reali_for_slope, REALIGNMENT_COLUMN, DA_COLUMN, n_bootstraps=N_BOOTSTRAPS)

# %%

box_data = [
    bootstrapped_behav_orig,
    bootstrapped_behav_realigned,
]

positions = [0, 1]
labels = ["Behavior Original", "Behavior Realigned"]
box_colors = ["grey", BEHAV_COLOR]

f, ax = plt.subplots(figsize=(1.8, 2.1),
                     gridspec_kw={"left": 0.35, "right": 0.95, "top": 0.8, "bottom": 0.25},)
bp = ax.boxplot(
    box_data,
    positions=positions,
    widths=0.55,
    showfliers=False,
    patch_artist=True,
)

for patch, color in zip(bp["boxes"], box_colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.6)
    patch.set_edgecolor("black")

for key in ["medians", "whiskers", "caps"]:
    for artist in bp[key]:
        artist.set_color("black")

ax.set_xticks([])
# ax.set_xticklabels(labels, rotation=25)
ax.set_ylabel("Slope")
ax.set_xlabel("")
ax.set_yticks([-0.1, 0, 0.1])
sns.despine(ax=ax, offset=5, bottom=True)

save_figure_atomic(f, FIGSFOLDER / "fig3_realignment_bootstrap_behaviour")

# %%
box_data = [
    bootstrapped_da_orig,
    bootstrapped_da_realigned,
]

positions = [0, 1]
labels = ["DA Original", "DA Realigned"]
box_colors = ["grey", DA_COLOR]

f, ax = plt.subplots(figsize=(1.8, 2.1),
                     gridspec_kw={"left": 0.35, "right": 0.95, "top": 0.8, "bottom": 0.25},)
bp = ax.boxplot(
    box_data,
    positions=positions,
    widths=0.55,
    showfliers=False,
    patch_artist=True,
)

for patch, color in zip(bp["boxes"], box_colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.6)
    patch.set_edgecolor("black")

for key in ["medians", "whiskers", "caps"]:
    for artist in bp[key]:
        artist.set_color("black")

ax.set_xticks([])
# ax.set_xticklabels(labels, rotation=25)
ax.set_ylabel("Slope")
ax.set_xlabel("")
ax.set_yticks([-1, 0, 1])
sns.despine(ax=ax, offset=5, bottom=True)

save_figure_atomic(f, FIGSFOLDER / "fig3_realignment_bootstrap_da")

# %%
stat, p = mannwhitneyu(bootstrapped_da_orig, bootstrapped_da_realigned, alternative="two-sided")
print(f"Dopamine - Mann-Whitney U: U = {stat:.1f}, p = {p:.4f}")

stat, p = mannwhitneyu(bootstrapped_behav_orig, bootstrapped_behav_realigned, alternative="two-sided")
print(f"Behavior - Mann-Whitney U: U = {stat:.1f}, p = {p:.4f}")

# %% [markdown]
# ## Old stuff below here

# %%
# plot data religned to behavior, full dopamine

DA_COLUMN = "auc_snips"
SIMBA_COLUMN = "simba_median_balance"
REALIGNMENT_COLUMN = "realigned_trials_behav"
RATS_TO_EXCLUDE = ["PB27", "PB71"]  # These rats have poor fits of behavioral data

subset_aligned = get_realigned_data(x_array, REALIGNMENT_COLUMN, rats_to_exclude=RATS_TO_EXCLUDE, verbose=False)
data_for_figure = prepare_shared_data(subset_aligned, REALIGNMENT_COLUMN, DA_COLUMN, SIMBA_COLUMN)


fig, ax, ax2 = make_realignment_panel_behav_and_da(data_for_figure, "trial", DA_COLUMN, SIMBA_COLUMN,
                                      da_color=DA_COLOR, behav_color=BEHAV_COLOR, behav_label=False
                                      )
if SAVE_FIGS:
    save_figure_atomic(fig, FIGSFOLDER / "figure_3_original_full_epoch")

fig, ax, ax2 = make_realignment_panel_behav_and_da(data_for_figure, REALIGNMENT_COLUMN, DA_COLUMN, SIMBA_COLUMN,
                                      orig_ks=False, da_color=DA_COLOR, behav_color=BEHAV_COLOR,
                                        da_label=False
                                      )
if SAVE_FIGS:
    save_figure_atomic(fig, FIGSFOLDER / "figure_3_realigned_full_epoch")




# %%
data_for_figure["df_original"]



# %%
# plot single parameters

DA_COLUMN = "auc_snips"
SIMBA_COLUMN = "simba_median_balance"
REALIGNMENT_COLUMN = "realigned_trials_behav"

# prepare realigned data
subset_aligned = get_realigned_data(x_array, REALIGNMENT_COLUMN,
                                    rats_to_exclude=[
                                        "PB27", # behavior does not appear sigmoidal even though it can be fitted
                                        # "PB48", # dopamine signal cannot be fitted with sigmoid
                                        "PB71", # behavior does not look sigmoidal
                                        ],
                                    # verbose=False
                                    )

data_for_figure = prepare_shared_data(subset_aligned, REALIGNMENT_COLUMN, DA_COLUMN, SIMBA_COLUMN)

fig, ax = make_realignment_panel_only_one_parameter(data_for_figure, REALIGNMENT_COLUMN, SIMBA_COLUMN, "behav",
                                          color=BEHAV_COLOR)

ax[0].set_yticks([-1, 0, 1])
save_figure_atomic(fig, FIGSFOLDER / "figure_3_realigned_behaviour_only")

fig, ax = make_realignment_panel_only_one_parameter(data_for_figure, REALIGNMENT_COLUMN, DA_COLUMN, "da",
                                          color=DA_COLOR)
ax[0].set_ylim(-10, 18)
save_figure_atomic(fig, FIGSFOLDER / "figure_3_realigned_da_fullepoch_only")

DA_COLUMN = "auc_snips_early"
data_for_figure = prepare_shared_data(subset_aligned, REALIGNMENT_COLUMN, DA_COLUMN, SIMBA_COLUMN)
fig, ax = make_realignment_panel_only_one_parameter(data_for_figure, REALIGNMENT_COLUMN, DA_COLUMN, "da",
                                          color=DA_COLOR)
ax[0].set_ylabel("Dopamine (Early) AUC")
ax[0].set_ylim(-4, 8)
save_figure_atomic(fig, FIGSFOLDER / "figure_3_realigned_da_early_epoch_only")

DA_COLUMN = "auc_snips_late"
data_for_figure = prepare_shared_data(subset_aligned, REALIGNMENT_COLUMN, DA_COLUMN, SIMBA_COLUMN)
fig, ax = make_realignment_panel_only_one_parameter(data_for_figure, REALIGNMENT_COLUMN, DA_COLUMN, "da",
                                          color=DA_COLOR)
ax[0].set_ylabel("Dopamine (Late) AUC")
ax[0].set_ylim(-5, 11)
save_figure_atomic(fig, FIGSFOLDER / "figure_3_realigned_da_late_epoch_only")

DA_COLUMN = "cluster_photo"

data_for_figure = prepare_shared_data(subset_aligned, REALIGNMENT_COLUMN, DA_COLUMN, SIMBA_COLUMN)
fig, ax = make_realignment_panel_only_one_parameter(data_for_figure, REALIGNMENT_COLUMN, DA_COLUMN, "da",
                                          color=DA_COLOR)
ax[0].set_ylabel("Cluster 1 Probability")
save_figure_atomic(fig, FIGSFOLDER / "figure_3_realigned_cluster_photo_only")


# %% [markdown]
# ## Normalized realignment

# %%
from figure_plotting import make_normalized_realignment_plot

DA_COLUMN = "cluster_photo"
SIMBA_COLUMN = "simba_median_balance"
REALIGNMENT_COLUMN = "realigned_trials_behav"

subset_aligned = get_realigned_data(x_array, REALIGNMENT_COLUMN, rats_to_exclude=["PB27", "PB71"], verbose=False)

data_for_figure = prepare_shared_data(subset_aligned, REALIGNMENT_COLUMN, DA_COLUMN, SIMBA_COLUMN)

fig, ax = make_normalized_realignment_plot(data_for_figure, da_color=DA_COLOR, behav_color=BEHAV_COLOR)
save_figure_atomic(fig, FIGSFOLDER / "figure_3_realigned_normalized_behav_and_da")

# %% [markdown]
# ## Bootstrapped histogram

# %%
from bootstrap_and_shuffle_helpers import get_bootstrapped_distribution

DA_COLUMN = "cluster_photo"
SIMBA_COLUMN = "simba_median_balance"
REALIGNMENT_COLUMN = "realigned_trials_behav"
N_BOOTSTRAPS = 100

subset_aligned = get_realigned_data(x_array, REALIGNMENT_COLUMN, rats_to_exclude=["PB27", "PB71"], verbose=False)

bootstrapped_behav_orig = get_bootstrapped_distribution(subset_aligned, "trial", SIMBA_COLUMN, n_bootstraps=N_BOOTSTRAPS)
bootstrapped_behav_realigned = get_bootstrapped_distribution(subset_aligned, REALIGNMENT_COLUMN, SIMBA_COLUMN, n_bootstraps=N_BOOTSTRAPS)

# %%
bootstrapped_da_orig = get_bootstrapped_distribution(subset_aligned, "trial", DA_COLUMN, n_bootstraps=N_BOOTSTRAPS)
bootstrapped_da_realigned = get_bootstrapped_distribution(subset_aligned, REALIGNMENT_COLUMN, DA_COLUMN, n_bootstraps=N_BOOTSTRAPS)

# %%
BANDWIDTH=0.25

fig, ax = plt.subplots(figsize=(2.3, 2.1),
                     gridspec_kw={"left": 0.25, "right": 0.9, "top": 0.8, "bottom": 0.25})

sns.kdeplot(bootstrapped_behav_orig, color="grey", ax=ax, cut=0, fill=True, bw_adjust=BANDWIDTH)
sns.kdeplot(bootstrapped_behav_realigned, color=BEHAV_COLOR, ax=ax, cut=0, fill=True, bw_adjust=BANDWIDTH)

ax.set_xlabel("Steepness (k)")
ax.set_xticks([-5, 0])
ax.set_xlim(-5.5, 0.5)
ax.set_yticks([0, 1])
sns.despine(ax=ax, offset=5)

save_figure_atomic(fig, FIGSFOLDER / "fig3_realigned_bootstrap_behaviour")

fig, ax = plt.subplots(figsize=(2.3, 2.1),
                     gridspec_kw={"left": 0.25, "right": 0.9, "top": 0.8, "bottom": 0.25})

sns.kdeplot(bootstrapped_da_orig, color="grey", ax=ax, cut=0, fill=True, bw_adjust=BANDWIDTH)
sns.kdeplot(bootstrapped_da_realigned, color=DA_COLOR, ax=ax, cut=0, fill=True, bw_adjust=BANDWIDTH)

ax.set_xlabel("Steepness (k)")
ax.set_xticks([-5, 0])
ax.set_xlim(-5.5, 0.5)
ax.set_yticks([0, 1])
sns.despine(ax=ax, offset=5)

save_figure_atomic(fig, FIGSFOLDER / "fig3_realigned_bootstrap_dopamine")

# %%
stat, p = mannwhitneyu(bootstrapped_da_orig, bootstrapped_da_realigned, alternative="two-sided")
print(f"Dopamine - Mann-Whitney U: U = {stat:.1f}, p = {p:.4f}")

stat, p = mannwhitneyu(bootstrapped_behav_orig, bootstrapped_behav_realigned, alternative="two-sided")
print(f"Behavior - Mann-Whitney U: U = {stat:.1f}, p = {p:.4f}")

# %% [markdown]
# ## Bootstrapped boxplots

# %%


box_data = [
    bootstrapped_da_orig,
    bootstrapped_da_realigned,
    bootstrapped_behav_orig,
    bootstrapped_behav_realigned,
]
positions = [0, 1, 3, 4]
labels = ["DA Original", "DA Realigned", "Behavior Original", "Behavior Realigned"]
box_colors = ["grey", DA_COLOR, "grey", BEHAV_COLOR]

f, ax = plt.subplots(figsize=(1.7, 2.1),
                     gridspec_kw={"left": 0.3, "right": 0.95, "top": 0.8, "bottom": 0.25},)
bp = ax.boxplot(
    box_data,
    positions=positions,
    widths=0.55,
    showfliers=False,
    patch_artist=True,
)

for patch, color in zip(bp["boxes"], box_colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.6)
    patch.set_edgecolor("black")

for key in ["medians", "whiskers", "caps"]:
    for artist in bp[key]:
        artist.set_color("black")

ax.set_xticks([])
# ax.set_xticklabels(labels, rotation=25)
ax.set_ylabel("Steepness of transition (k)")
ax.set_xlabel("")
sns.despine(ax=ax, offset=5, bottom=True)

# if SAVE_FIGS:
#     save_figure_atomic(f, "fig3_realignment_bootstrap_comparison", FIGSFOLDER)




# %%
np.isnan(bootstrapped_control_da).sum(), np.isnan(bootstrapped_expt_da).sum(), np.isnan(bootstrapped_control_behav).sum(), np.isnan(bootstrapped_expt_behav).sum()

# %%
## To make bootstrapped distribution of k values from original (unaligned) data
## AND bootstrapped distribution of k values from realigned data, to compare against each other and against shuffled controls
## e.g. should use same df and fitting procedure but just resample animals with replacement instead of shuffling x0 values

rats = subset_aligned.id.unique()

import random

num_bootstraps = 1000    # how many random lists you want
list_len = len(rats)     # sample the same number of animals as the active transition source

random_lists = [random.choices(list(rats), k=list_len) for _ in range(num_bootstraps)]

bootstrapped_ks_da_orig = []
for i, random_list in enumerate(random_lists, 1):
    subset_bootstrapped = []
    for rat in random_list:
        rat_trials = subset_aligned[subset_aligned.id == rat]
        subset_bootstrapped.append(rat_trials)
    boot_df = pd.concat(subset_bootstrapped)
    popt_boot = fit_sigmoid(boot_df, "trial", da_column, initial_k=-1)
    k_boot = popt_boot[2] if np.all(np.isfinite(popt_boot)) else np.nan
    bootstrapped_ks_da_orig.append(k_boot)
    # print(f"Bootstrap {i}: k = {k_boot:.3f}")
    
bootstrapped_ks_da_realigned = []
for i, random_list in enumerate(random_lists, 1):
    n_unique = len(set(random_list))   # actual unique animals in this bootstrap draw
    subset_bootstrapped = []
    for rat in random_list:
        rat_trials = subset_aligned[subset_aligned.id == rat]
        subset_bootstrapped.append(rat_trials)
    boot_df = pd.concat(subset_bootstrapped)
    boot_df_realigned = make_realigned_trials(boot_df, fits_df, verbose=False)
    boot_df_realigned_complete, _ = filter_balanced_trials(
        boot_df_realigned,
        target_count=n_unique,
        fallback_to_max=True,
    )
    popt_boot = fit_sigmoid(boot_df_realigned_complete, "trial_aligned", da_column, initial_k=-1)
    k_boot = popt_boot[2] if np.all(np.isfinite(popt_boot)) else np.nan
    bootstrapped_ks_da_realigned.append(k_boot)
    # print(f"Bootstrap {i}: k = {k_boot:.3f}")

bootstrapped_ks_da_orig = np.array([k for k in bootstrapped_ks_da_orig if not np.isnan(k)])
bootstrapped_ks_da_realigned = np.array([k for k in bootstrapped_ks_da_realigned if not np.isnan(k)])


# %%

# %%

# %%
da_column = "auc_snips"
da_column = "auc_snips_early"
da_column = "auc_snips_late"
da_column = "cluster_photo"

# # Compare original vs realigned AUC
f, ax = plt.subplots(ncols=2, figsize=(6, 2.5), sharey=True)

# Left plot: Original trial order
popt_orig = plot_auc_and_sigmoid(subset_aligned, "trial", da_column, ax=ax[0], include_steepness=True)
ax[0].set_xlabel('Trial Number', fontsize=10)
ax[0].set_title('Original Trial Order', fontsize=11)
ax[0].set_ylabel('Dopamine AUC', fontsize=10)

# Right plot: Realigned to transition
popt_realigned = plot_auc_and_sigmoid(include_only_complete_trials(subset_aligned_8), "trial_aligned", da_column, ax=ax[1], 
                                      first_trial=int(include_only_complete_trials(subset_aligned_8).trial_aligned.min()),
                                      include_steepness=True)
ax[1].set_xlabel('Trial Relative to Transition', fontsize=10)
ax[1].set_title('Realigned to Transition Point', fontsize=11)


# Mark transition point
ax[1].axvline(0, color='gray', linestyle='--', alpha=0.5, linewidth=1)

sns.despine(ax=ax[0], offset=5)
sns.despine(ax=ax[1], offset=5)

if SAVE_FIGS:
    save_figure_atomic(f, "fig4b_realigned_transitions", FIGSFOLDER)


print(f"Original order - transition midpoint: {popt_orig[1]:.1f}, steepness (k): {popt_orig[2]:.3f}")
print(f"Realigned - transition midpoint: {popt_realigned[1]:.1f}, steepness (k): {popt_realigned[2]:.3f}")

# %% [markdown]
# ## Figure 4-2: Dopamine and Behavior Co-Aligned to Transition

# %% [markdown]
# ## Realigning to dopamine transitions

# %%
# plot single parameters

DA_COLUMN = "auc_snips"
SIMBA_COLUMN = "simba_median_balance"
REALIGNMENT_COLUMN = "realigned_trials_da"


subset_aligned = get_realigned_data(x_array, REALIGNMENT_COLUMN,
                                    # rats_to_exclude=["PB27", "PB71"],
                                    #verbose=False
                                    )
data_for_figure = prepare_shared_data(subset_aligned, REALIGNMENT_COLUMN, DA_COLUMN, SIMBA_COLUMN)

fig, ax = make_realignment_panel_only_one_parameter(data_for_figure, REALIGNMENT_COLUMN, DA_COLUMN, "da",
                                          color=DA_COLOR)

DA_COLUMN = "cluster_photo"

data_for_figure = prepare_shared_data(subset_aligned, REALIGNMENT_COLUMN, DA_COLUMN, SIMBA_COLUMN)
fig, ax = make_realignment_panel_only_one_parameter(data_for_figure, REALIGNMENT_COLUMN, DA_COLUMN, "da",
                                          color=DA_COLOR)
ax[0].set_ylabel("Cluster 1 Probability")


DA_COLUMN = "auc_snips_early"
data_for_figure = prepare_shared_data(subset_aligned, REALIGNMENT_COLUMN, DA_COLUMN, SIMBA_COLUMN)
fig, ax = make_realignment_panel_only_one_parameter(data_for_figure, REALIGNMENT_COLUMN, DA_COLUMN, "da",
                                          color=DA_COLOR)
ax[0].set_ylabel("Dopamine (Early) AUC")

DA_COLUMN = "auc_snips_late"
data_for_figure = prepare_shared_data(subset_aligned, REALIGNMENT_COLUMN, DA_COLUMN, SIMBA_COLUMN)
fig, ax = make_realignment_panel_only_one_parameter(data_for_figure, REALIGNMENT_COLUMN, DA_COLUMN, "da",
                                          color=DA_COLOR)
ax[0].set_ylabel("Dopamine (Late) AUC")
ax[0].set_ylim(-5, 11)


fig, ax = make_realignment_panel_only_one_parameter(data_for_figure, REALIGNMENT_COLUMN, SIMBA_COLUMN, "behav",
                                          color=BEHAV_COLOR)

ax[0].set_yticks([-1, 0, 1])


# %%


df_complete = fig4_shared["df_complete"]

f, ax = plt.subplots(figsize=(2.3, 2.1),
                     gridspec_kw={"left": 0.3, "right": 0.7, "top": 0.8, "bottom": 0.25},)

popt_orig = plot_auc_and_sigmoid(subset_aligned, "trial", da_column, ax=ax, include_steepness=False)
ax.set_xlabel('Trial', fontsize=10)
ax.set_ylabel('Dopamine AUC', fontsize=10, color=da_color)
ax.tick_params(axis='y', labelcolor=da_color)
ax.text(0.5, 1.1,
        f"Dopamine k = {popt_orig[2]:.3f}\nBehavior k = {fig4_shared['k_behav_orig']:.3f}",
        ha="center", transform=ax.transAxes)

ax2 = ax.twinx()
plot_realigned_behaviour(subset_aligned, "trial", ax=ax2, include_steepness=False, simba_col=simba_col)
# ax2.set_ylabel('Time Moving', fontsize=10, color=behav_color)
ax2.tick_params(axis='y', labelcolor=behav_color)

if simba_col == "simba_raw_mean":
    ax2.set_ylim(0.1, 0.4)  # Set y-limits for raw metric
    ax2.set_yticks((0.2, 0.4), labels=["0.2", "0.4"])
else:
    ax2.set_ylim([-0.70, 1.1])
    ax2.set_yticks([-0.5, 0, 0.5, 1.0], labels=["","","",""])
sns.despine(ax=ax2, offset=5, left=True, right=False)
ax2.spines['right'].set_position(('outward', 5))

sns.despine(ax=ax, offset=5)

if SAVE_FIGS:
    save_figure_atomic(f, "fig4a_realigned_transitions_with_behavior_panel1", FIGSFOLDER)

print(f"Panel 1 (original) — Dopamine k: {popt_orig[2]:.3f}, Behavior k: {fig4_shared['k_behav_orig']:.3f}")

# %%
# Panel 2: Realigned to transition (dopamine + behavior)
if "fig4_shared" not in locals():
    fig4_shared = prepare_fig4_shared_data(subset_aligned)
df_complete = fig4_shared["df_complete"]

f, ax = plt.subplots(figsize=(2.3, 2.1),
                     gridspec_kw={"left": 0.3, "right": 0.7, "top": 0.8, "bottom": 0.25},)

popt_realigned = plot_auc_and_sigmoid(
    df_complete,
    "trial_aligned",
    da_column,
    ax=ax,
    first_trial=int(df_complete.trial_aligned.min()),
    include_steepness=False,
 )
ax.set_xlabel('Trial after realignment', fontsize=10)
# ax.set_ylabel('Dopamine AUC', fontsize=10, color=colors[3])
ax.tick_params(axis='y', labelcolor=da_color)
ax.axvline(0, color='gray', linestyle='--', alpha=0.5, linewidth=1)
ax.set_xticks([-10, 0, 10, 20])
ax.set_yticks([-50,0,50,100], labels=["","","",""])
ax.text(
    0.5,
    1.1,
    f"Dopamine k = {popt_realigned[2]:.3f}\nBehavior k = {fig4_shared['k_behav_realigned']:.3f}",
    ha="center",
    transform=ax.transAxes
)

ax2 = ax.twinx()
plot_realigned_behaviour(
    df_complete,
    "trial_aligned",
    ax=ax2,
    first_trial=int(df_complete.trial_aligned.min()),
    include_steepness=False,
 )
ax2.set_ylabel('Appetitive probability', fontsize=10, color=behav_color, rotation=270, labelpad=10)
ax2.tick_params(axis='y', labelcolor=behav_color)

if simba_col == "simba_raw_mean":
    ax2.set_ylim(0.1, 0.4)  # Set y-limits for raw metric
    ax2.set_yticks((0.2, 0.4), labels=["0.2", "0.4"])
else:
    ax2.set_ylim([-0.70, 1.1])
    ax2.set_yticks([-0.5, 0, 0.5, 1.0], labels=["","","",""])
sns.despine(ax=ax2, offset=5, left=True, right=False)
ax2.spines['right'].set_position(('outward', 5))

sns.despine(ax=ax, offset=5)

if SAVE_FIGS:
    save_figure_atomic(f, "fig4b_realigned_transitions_with_behavior", FIGSFOLDER)

print(f"Panel 2 (realigned) — Dopamine k: {popt_realigned[2]:.3f}, Behavior k: {fig4_shared['k_behav_realigned']:.3f}")

# %%
data_for_figure.keys()


# %%

# %% [markdown]
# ## Figure 4D: Shuffled data

# %%
# Legacy exploratory block removed.
# Canonical fitting utilities are defined above in the helper section:
# - sigmoid (model function)
# - fit_sigmoid (bounded + validated fitter)
# Shuffled analyses should use get_shuffled_k below.

# %%
def shuffle_without_fixed_values(values, max_attempts=10000, rng=None):
    values = np.asarray(values)
    rng = np.random.default_rng() if rng is None else rng

    for _ in range(max_attempts):
        shuffled = rng.permutation(values)
        if np.all(shuffled != values):
            return shuffled

    raise ValueError(
        "Could not find a valid shuffle where every value changes position. "
        "This usually means duplicate values make the constraint impossible or very unlikely."
    )


def filter_balanced_trials(df, group_col="trial_aligned", target_count=None, fallback_to_max=False):
    """Keep bins with equal representation across animals.

    If target_count bins are unavailable and fallback_to_max=True,
    fall back to the largest achievable equal-count subset."""
    if df.empty:
        return df.copy(), 0

    counts = df.groupby(group_col).size()
    if counts.empty:
        return df.copy(), 0

    if target_count is None:
        target_count = int(counts.max())

    filtered = df.groupby(group_col, group_keys=False).filter(lambda g: len(g) == target_count).copy()
    if (not filtered.empty) or (not fallback_to_max):
        return filtered, int(target_count)

    max_count = int(counts.max())
    filtered = df.groupby(group_col, group_keys=False).filter(lambda g: len(g) == max_count).copy()
    return filtered, max_count


def get_shuffled_k(
    df,
    parameter,
    num_repeats=1000,
    num_runs=1,
    initial_k=-1,
    k_min=-5,
    k_max=5,
    x0_min=None,
    x0_max=None,
    asymptote_tol_frac=0.25,
    enforce_complete_trials=True,
    validate_real_fit=True,
    validate_shuffled_fits=True,  # reject boundary-stuck / degenerate shuffled fits
    verbose=False,
):
    base_df = df.copy()
    rat_ids = np.asarray(sorted(base_df["id"].unique()))
    n_required = len(rat_ids)

    if n_required < 2:
        raise ValueError("Need at least two animals to compute shuffled transition controls.")

    z = base_df.copy()
    real_n_used = n_required
    if enforce_complete_trials:
        z, real_n_used = filter_balanced_trials(
            z,
            target_count=n_required,
            fallback_to_max=False,
        )

    real_k = fit_sigmoid(
        z,
        "trial_aligned",
        parameter,
        initial_k=initial_k,
        k_min=k_min,
        k_max=k_max,
        x0_min=x0_min,
        x0_max=x0_max,
        asymptote_tol_frac=asymptote_tol_frac,
        validate_fit=validate_real_fit,
    )[2]

    fits_selected = (
        fits_df
        .loc[fits_df["id"].isin(rat_ids), ["id", "x0_orig"]]
        .drop_duplicates(subset=["id"])
        .sort_values("id")
        .reset_index(drop=True)
    )

    if len(fits_selected) != n_required:
        raise ValueError(
            f"Expected fits for {n_required} animals from the active transition source, "
            f"but found {len(fits_selected)}. Re-run the load/prep cells for the current toggle setting."
        )

    original_x0s = fits_selected["x0_orig"].to_numpy(dtype=float)
    shuffled_k = np.full(num_repeats, np.nan, dtype=float)
    shuffled_x0_abs_distances = np.full((num_repeats, original_x0s.size), np.nan, dtype=float)

    for i in range(num_repeats):
        shuffled_x0s = shuffle_without_fixed_values(original_x0s)
        shuffled_x0_abs_distances[i] = np.abs(original_x0s - shuffled_x0s)

        df_fitted_params_shuffled = fits_selected.assign(x0_orig=shuffled_x0s)
        z_temp = make_realigned_trials(base_df.copy(), df_fitted_params_shuffled, verbose=False)

        shuffled_n_used = n_required
        if enforce_complete_trials:
            z_temp, shuffled_n_used = filter_balanced_trials(
                z_temp,
                target_count=n_required,
                fallback_to_max=True,
            )

        if z_temp.empty:
            if verbose:
                print(f"shuffle {i}: no balanced bins available")
            continue

        tmp_ks = np.full(num_runs, np.nan, dtype=float)
        for j in range(num_runs):
            tmp_ks[j] = fit_sigmoid(
                z_temp,
                "trial_aligned",
                parameter,
                initial_k=initial_k,
                k_min=k_min,
                k_max=k_max,
                x0_min=x0_min,
                x0_max=x0_max,
                asymptote_tol_frac=asymptote_tol_frac,
                validate_fit=validate_shuffled_fits,
            )[2]

        valid_ks = tmp_ks[np.isfinite(tmp_ks)]
        if valid_ks.size > 0:
            shuffled_k[i] = np.median(valid_ks)

        if verbose:
            print(f"shuffle {i}: k={shuffled_k[i]}, n_used={shuffled_n_used}, real_n_used={real_n_used}")

    return real_k, shuffled_k, shuffled_x0_abs_distances


# Helper definitions only. Run the DA/behavior shuffle cells below to compute results.

# %%
fits_df.x0_orig.values

# %%
real_k_da, shuffled_k_da, shuffled_x0_abs_distances_da = get_shuffled_k(
    subset_aligned, da_column, num_repeats=500, initial_k=-1
)


# %%
def plot_shuffled_ks(real_k, shuffled_k, ax=None, color_real=da_color, color_shuf=da_color, kde=False):
    if ax is None:
        f, ax = plt.subplots(figsize=(3, 2.5))

    finite_k = shuffled_k[np.isfinite(shuffled_k)]
    lo = min(np.floor(np.nanmin(finite_k) * 2) / 2, -2) if finite_k.size else -2
    hi = max(np.ceil(np.nanmax(finite_k) * 2) / 2, 0.5) if finite_k.size else 0.5
    bins = np.linspace(lo, hi, 50)
    sns.histplot(finite_k, bins=bins,
                 kde=kde,
                    ax=ax, color=color_shuf, alpha=0.5,
    )

    ax.axvline(real_k, ymax=0.8, linestyle="--", color=color_real)
    sns.despine(ax=ax)

    ax.set_ylabel("Frequency")

    return f, ax

f, ax = plot_shuffled_ks(real_k_da, shuffled_k_da, color_real=da_color, color_shuf=da_color)
ax.set_xlabel("DA k after realignment")

unaligned_k = fig4_shared["k_da_orig"]

real_k_diff = real_k_da - unaligned_k
shuffled_k_diff = shuffled_k_da - unaligned_k

f, ax = plot_shuffled_ks(real_k_diff, shuffled_k_diff, color_real=da_color, color_shuf=da_color)
ax.set_xlabel(" $\Delta$ k after realignment")

if shuffled_k_da.size > 0:
    print("p = {}".format(1 - np.sum(shuffled_k_da > real_k_da) / shuffled_k_da.size))
else:
    print("p = nan (no valid shuffled fits)")

if SAVE_FIGS:
    save_figure_atomic(f, "fig4d_realignment_shuffle_control_dopamine", FIGSFOLDER)

# %%
real_k_behav, shuffled_k_behav, shuffled_x0_abs_distances_behav = get_shuffled_k(
    subset_aligned, simba_col, num_repeats=500, initial_k=-1
)

# %%
f, ax = plot_shuffled_ks(real_k_behav, shuffled_k_behav, color_real=behav_color, color_shuf=behav_color)
ax.set_xlabel("Behavior k after realignment")

unaligned_k = fig4_shared["k_behav_orig"]

f, ax = plot_shuffled_ks(real_k_behav - unaligned_k, shuffled_k_behav - unaligned_k, color_real=behav_color, color_shuf=behav_color)
ax.set_xlabel(" $\Delta$ k after realignment")

if shuffled_k_behav.size > 0:
    print("p = {}".format(1 - np.sum(shuffled_k_behav > real_k_behav) / shuffled_k_behav.size))
else:
    print("p = nan (no valid shuffled fits)")

if SAVE_FIGS:
    save_figure_atomic(f, "fig4d_realignment_shuffle_control_behavior", FIGSFOLDER)

print(real_k_behav)

# %%
shuffled_k_behav

# %%
## To make bootstrapped distribution of k values from original (unaligned) data
## AND bootstrapped distribution of k values from realigned data, to compare against each other and against shuffled controls
## e.g. should use same df and fitting procedure but just resample animals with replacement instead of shuffling x0 values

rats = subset_aligned.id.unique()

import random

num_bootstraps = 1000    # how many random lists you want
list_len = len(rats)     # sample the same number of animals as the active transition source

random_lists = [random.choices(list(rats), k=list_len) for _ in range(num_bootstraps)]

bootstrapped_ks_da_orig = []
for i, random_list in enumerate(random_lists, 1):
    subset_bootstrapped = []
    for rat in random_list:
        rat_trials = subset_aligned[subset_aligned.id == rat]
        subset_bootstrapped.append(rat_trials)
    boot_df = pd.concat(subset_bootstrapped)
    popt_boot = fit_sigmoid(boot_df, "trial", da_column, initial_k=-1)
    k_boot = popt_boot[2] if np.all(np.isfinite(popt_boot)) else np.nan
    bootstrapped_ks_da_orig.append(k_boot)
    # print(f"Bootstrap {i}: k = {k_boot:.3f}")
    
bootstrapped_ks_da_realigned = []
for i, random_list in enumerate(random_lists, 1):
    n_unique = len(set(random_list))   # actual unique animals in this bootstrap draw
    subset_bootstrapped = []
    for rat in random_list:
        rat_trials = subset_aligned[subset_aligned.id == rat]
        subset_bootstrapped.append(rat_trials)
    boot_df = pd.concat(subset_bootstrapped)
    boot_df_realigned = make_realigned_trials(boot_df, fits_df, verbose=False)
    boot_df_realigned_complete, _ = filter_balanced_trials(
        boot_df_realigned,
        target_count=n_unique,
        fallback_to_max=True,
    )
    popt_boot = fit_sigmoid(boot_df_realigned_complete, "trial_aligned", da_column, initial_k=-1)
    k_boot = popt_boot[2] if np.all(np.isfinite(popt_boot)) else np.nan
    bootstrapped_ks_da_realigned.append(k_boot)
    # print(f"Bootstrap {i}: k = {k_boot:.3f}")

bootstrapped_ks_da_orig = np.array([k for k in bootstrapped_ks_da_orig if not np.isnan(k)])
bootstrapped_ks_da_realigned = np.array([k for k in bootstrapped_ks_da_realigned if not np.isnan(k)])


# %%
## To make bootstrapped distribution of k values from original (unaligned) data
## AND bootstrapped distribution of k values from realigned data, to compare against each other and against shuffled controls
## e.g. should use same df and fitting procedure but just resample animals with replacement instead of shuffling x0 values

rats = subset_aligned.id.unique()

import random

num_bootstraps = 100    # how many random lists you want
list_len = len(rats)     # sample the same number of animals as the active transition source

random_lists = [random.choices(list(rats), k=list_len) for _ in range(num_bootstraps)]

bootstrapped_ks_behav_orig = []
for i, random_list in enumerate(random_lists, 1):
    subset_bootstrapped = []
    for rat in random_list:
        rat_trials = subset_aligned[subset_aligned.id == rat]
        subset_bootstrapped.append(rat_trials)
    boot_df = pd.concat(subset_bootstrapped)
    popt_boot = fit_sigmoid(boot_df, "trial", simba_col, initial_k=-1)
    k_boot = popt_boot[2] if np.all(np.isfinite(popt_boot)) else np.nan
    bootstrapped_ks_behav_orig.append(k_boot)
    # print(f"Bootstrap {i}: k = {k_boot:.3f}")
    
bootstrapped_ks_behav_realigned = []
for i, random_list in enumerate(random_lists, 1):
    n_unique = len(set(random_list))   # actual unique animals in this bootstrap draw
    subset_bootstrapped = []
    for rat in random_list:
        rat_trials = subset_aligned[subset_aligned.id == rat]
        subset_bootstrapped.append(rat_trials)
    boot_df = pd.concat(subset_bootstrapped)
    boot_df_realigned = make_realigned_trials(boot_df, fits_df, verbose=False)
    boot_df_realigned_complete, _ = filter_balanced_trials(
        boot_df_realigned,
        target_count=n_unique,
        fallback_to_max=True,
    )
    popt_boot = fit_sigmoid(boot_df_realigned_complete, "trial_aligned", simba_col, initial_k=-1)
    k_boot = popt_boot[2] if np.all(np.isfinite(popt_boot)) else np.nan
    bootstrapped_ks_behav_realigned.append(k_boot)
    # print(f"Bootstrap {i}: k = {k_boot:.3f}")

bootstrapped_ks_behav_orig = np.array([k for k in bootstrapped_ks_behav_orig if not np.isnan(k)])
bootstrapped_ks_behav_realigned = np.array([k for k in bootstrapped_ks_behav_realigned if not np.isnan(k)])


# %%
# horizontal box plots


# %%


box_data = [
    bootstrapped_ks_da_orig,
    bootstrapped_ks_da_realigned,
    bootstrapped_ks_behav_orig,
    bootstrapped_ks_behav_realigned,
]
positions = [0, 1, 3, 4]
labels = ["DA Original", "DA Realigned", "Behavior Original", "Behavior Realigned"]
box_colors = ["grey", da_color, "grey", behav_color]

f, ax = plt.subplots(figsize=(1.7, 2.1),
                     gridspec_kw={"left": 0.3, "right": 0.95, "top": 0.8, "bottom": 0.25},)
bp = ax.boxplot(
    box_data,
    positions=positions,
    widths=0.55,
    showfliers=False,
    patch_artist=True,
)

for patch, color in zip(bp["boxes"], box_colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.6)
    patch.set_edgecolor("black")

for key in ["medians", "whiskers", "caps"]:
    for artist in bp[key]:
        artist.set_color("black")

ax.set_xticks([])
# ax.set_xticklabels(labels, rotation=25)
ax.set_ylabel("Steepness of transition (k)")
ax.set_xlabel("")
sns.despine(ax=ax, offset=5, bottom=True)

if SAVE_FIGS:
    save_figure_atomic(f, "fig3_realignment_bootstrap_comparison", FIGSFOLDER)




# %%
f, ax = plt.subplots(figsize=(1.7, 2.3), nrows=2, sharex=True,
                     gridspec_kw={"left": 0.2, "right": 0.95, "top": 0.8, "bottom": 0.25,
                                  },)

sns.kdeplot(bootstrapped_ks_da_orig, color="grey", fill=True, alpha=0.5, label="DA Original", ax=ax[0], cut=0)
sns.kdeplot(bootstrapped_ks_da_realigned, color=da_color, fill=True, alpha=0.5, label="DA Realigned", ax=ax[0], cut=0)


for axis in ax:
    sns.despine(ax=axis, offset=5)
    axis.set_yticks([])
    axis.set_xlim(-2.1, 0.5)
    
ax[1].set_xlabel("Steepness of transition (k)")

if SAVE_FIGS:
    save_figure_atomic(f, "fig3_realignment_bootstrap_comparison_kde", FIGSFOLDER)
    
    
stat, p = mannwhitneyu(bootstrapped_ks_da_orig, bootstrapped_ks_da_realigned, alternative="two-sided")
print(f"Dopamine - Mann-Whitney U: U = {stat:.1f}, p = {p:.4f}")

stat, p = mannwhitneyu(bootstrapped_ks_behav_orig, bootstrapped_ks_behav_realigned, alternative="two-sided")
print(f"Behavior - Mann-Whitney U: U = {stat:.1f}, p = {p:.4f}")

# %%
f, ax = plt.subplots(figsize=(1.7, 2.3), nrows=3, sharex=True,
                     gridspec_kw={"left": 0.2, "right": 0.95, "top": 0.8, "bottom": 0.25,
                                  "height_ratios": [1, 0.2, 1]},)

sns.kdeplot(bootstrapped_ks_da_orig, color="grey", fill=True, alpha=0.5, label="DA Original", ax=ax[0], cut=0)
sns.kdeplot(bootstrapped_ks_da_realigned, color=da_color, fill=True, alpha=0.5, label="DA Realigned", ax=ax[0], cut=0)

sns.kdeplot(bootstrapped_ks_behav_orig, color="grey", fill=True, alpha=0.5, label="Behavior Original", ax=ax[2], cut=0)
sns.kdeplot(bootstrapped_ks_behav_realigned, color=behav_color, fill=True, alpha=0.5, label="Behavior Realigned", ax=ax[2], cut=0)

ax[1].axis("off")

for axis in ax:
    sns.despine(ax=axis, offset=5)
    axis.set_yticks([])
    axis.set_xlim(-2.1, 0.5)
    
ax[2].set_xlabel("Steepness of transition (k)")

if SAVE_FIGS:
    save_figure_atomic(f, "fig3_realignment_bootstrap_comparison_kde", FIGSFOLDER)
    
    
stat, p = mannwhitneyu(bootstrapped_ks_da_orig, bootstrapped_ks_da_realigned, alternative="two-sided")
print(f"Dopamine - Mann-Whitney U: U = {stat:.1f}, p = {p:.4f}")

stat, p = mannwhitneyu(bootstrapped_ks_behav_orig, bootstrapped_ks_behav_realigned, alternative="two-sided")
print(f"Behavior - Mann-Whitney U: U = {stat:.1f}, p = {p:.4f}")

# %%
bootstrapped_ks = np.array(bootstrapped_ks)

bootstrapped_ks = bootstrapped_ks[bootstrapped_ks > -2.0]

if bootstrapped_ks.size > 0:
    print("p = {}".format(1 - np.sum(bootstrapped_ks > real_k_behav) / bootstrapped_ks.size))
else:
    print("p = nan (no valid shuffled fits)")

# %%
# plotted as difference from unaligned k
unaligned_k = fig4_shared["k_behav_orig"]

valid_shuffled_k_behav = shuffled_k_behav[np.isfinite(shuffled_k_behav)]
shuf_k = valid_shuffled_k_behav - unaligned_k
real_k = real_k_behav - unaligned_k

bins = np.linspace(-2, 2, 80)
f, ax = plt.subplots(figsize=(3, 2.5),
                     gridspec_kw={"left": 0.25, "right": 0.95, "top": 0.8, "bottom": 0.3},)

# ax.hist(shuf_k, bins=bins, color=color_shuf)
sns.histplot(shuf_k, bins=bins,
             # kde=True,
             ax=ax, color=behav_color, alpha=0.5,
             #edgecolor='none'
             )
ax.axvline(real_k, ymax=0.8, linestyle="--", color=behav_color)

sns.despine(ax=ax, offset=5)
ax.set_xlim(-1.1, 1.1)
# ax.set_xticks([0, 1, 2])
# ax.text(real_k, 110,
#         "Real k, p = {:2.3f}".format(1 - np.sum(shuf_k < real_k) / shuf_k.size),
#         ha="center",
#         color=color_real)

# ax.text(0.5, 70, "Shuffled ks", color=color_shuf)

ax.set_ylabel("Frequency")
ax.set_xlabel("Steepness of transition (k)")

if shuf_k.size > 0:
    print("p = {}".format(1 - np.sum(shuf_k > real_k) / shuf_k.size))
else:
    print("p = nan (no valid shuffled fits)")

if SAVE_FIGS:
    save_figure_atomic(f, "fig4d_realignment_shuffle_control_behavior", FIGSFOLDER)


# %%
1 - np.sum(shuf_k > real_k) / 1000

# %%
# Plot before/after AUC comparison
f, axes = plt.subplots(1, 2, figsize=(5, 3))

# Plot 1: Individual animal comparisons
ax = axes[0]
x_pos = np.arange(len(auc_df))
width = 0.35

bars1 = ax.bar(x_pos - width/2, auc_df['auc_before'], width, label='Before Transition',
                color=colors[0], alpha=0.7, edgecolor='k', linewidth=1.5)
bars2 = ax.bar(x_pos + width/2, auc_df['auc_after'], width, label='After Transition',
                color=colors[1], alpha=0.7, edgecolor='k', linewidth=1.5)

ax.set_xlabel('Animal', fontsize=10)
ax.set_ylabel('Mean Photometry AUC', fontsize=10)
ax.set_title('Photometry AUC Before/After Transition', fontsize=11)
ax.set_xticks(x_pos)
ax.set_xticklabels(auc_df['animal'], rotation=45, ha='right', fontsize=9)
ax.legend(fontsize=9)
sns.despine(ax=ax)

# Plot 2: Change in AUC
ax = axes[1]
auc_df['auc_change'] = auc_df['auc_after'] - auc_df['auc_before']
colors_change = [colors[1] if x > 0 else colors[0] for x in auc_df['auc_change']]

ax.bar(x_pos, auc_df['auc_change'], color=colors_change, alpha=0.7, edgecolor='k', linewidth=1.5)
ax.axhline(0, color='k', linestyle='-', linewidth=0.5, alpha=0.5)
ax.set_xlabel('Animal', fontsize=10)
ax.set_ylabel('AUC Change (After - Before)', fontsize=10)
ax.set_title('Dopamine Response Change at Transition', fontsize=11)
ax.set_xticks(x_pos)
ax.set_xticklabels(auc_df['animal'], rotation=45, ha='right', fontsize=9)
sns.despine(ax=ax)

plt.tight_layout()
if SAVE_FIGS:
    save_figure(f, "fig4b_auc_before_after_transition", FIGSFOLDER)
plt.show()

# %% [markdown]
# ## Export Results

# %%
# Export AUC summary
auc_df.to_csv(RESULTSFOLDER / "transition_aligned_auc_summary.csv", index=False)
print(f"Exported AUC summary to {RESULTSFOLDER / 'transition_aligned_auc_summary.csv'}")

print(f"\nFigure 4B generation complete!")
print(f"Summary:")
print(f"  - {len(animals)} animals with valid transition fits")
print(f"  - {len(subset_aligned)} total trials used for alignment")
print(f"  - Mean AUC change: {auc_df['auc_change'].mean():.3f} ± {auc_df['auc_change'].std():.3f}")


# %%
# Alternative null distribution: random rat-level realignment points (with replacement)
def get_random_x0_null_distribution(
    df,
    parameter=da_column,
    num_repeats=100,
    x0_low=11,
    x0_high=26,
    initial_k=-1,
    k_min=-2,
    k_max=2,
    enforce_complete_trials=True,
    validate_fit=False,
    rng=None,
):
    """Build null by drawing random x0 per rat each repeat, then realigning and fitting sigmoid."""
    rng = np.random.default_rng() if rng is None else rng

    base_df = df.copy()
    params_template = fits_df.copy()

    if "id" in params_template.columns:
        rat_ids = np.asarray(sorted(params_template["id"].dropna().unique()))
        id_based = True
    else:
        rat_ids = np.arange(len(params_template))
        id_based = False

    sampled_x0 = np.full((num_repeats, len(rat_ids)), np.nan, dtype=float)
    null_k = np.full(num_repeats, np.nan, dtype=float)

    for i in range(num_repeats):
        x0_draw = rng.integers(x0_low, x0_high + 1, size=len(rat_ids))
        sampled_x0[i] = x0_draw

        params_random = params_template.copy()
        if id_based:
            x0_map = dict(zip(rat_ids, x0_draw))
            params_random["x0_orig"] = params_random["id"].map(x0_map).astype(float)
        else:
            params_random["x0_orig"] = x0_draw.astype(float)

        z_temp = make_realigned_trials(base_df.copy(), params_random, verbose=False)
        if enforce_complete_trials:
            z_temp = include_only_complete_trials(z_temp)

        null_k[i] = fit_sigmoid(
            z_temp,
            "trial_aligned",
            parameter,
            initial_k=initial_k,
            k_min=k_min,
            k_max=k_max,
            validate_fit=validate_fit,
        )[2]

    return null_k, sampled_x0, rat_ids


# Example run (dopamine)
null_k_random_behav, sampled_x0_random_behav, rat_ids_random = get_random_x0_null_distribution(
    subset_aligned,
    parameter=simba_col,
    num_repeats=1000,
    x0_low=11,
    x0_high=26,
    initial_k=-1,
    validate_fit=False,
)

print(f"Random-x0 null ks (DA): {np.isfinite(null_k_random_da).sum()} / {null_k_random_da.size} valid")
print(f"Sampled x0 matrix shape: {sampled_x0_random_da.shape} (repeats x rats)")
print(f"x0 draw range observed: {int(np.nanmin(sampled_x0_random_da))} to {int(np.nanmax(sampled_x0_random_da))}")

# %%
f, ax = plot_shuffled_ks(real_k_behav, null_k_random_behav, color_real=behav_color, color_shuf=behav_color)
ax.set_xlim(-2, 2)
ax.set_xticks([-2, -1, 0])
ax.set_xlabel("Behavior k after realignment")

if null_k_random_behav.size > 0:
    p_random = 1 - np.sum(null_k_random_behav > real_k_behav) / null_k_random_behav.size
    print(f"Random-x0 null p-value (behavior): {p_random:.4f}")

# %%

# %% [markdown]
# ## Realign Dopamine and Behavior to Behavior-Derived Transition Points
#
# Estimate per-animal transition points from a behavioral signal, keep steep negative slopes, reject positive slopes, and realign dopamine (`auc_snips`) and behavior to those transition points.

# %%
# Build transition points from behavior, then realign dopamine and behavior
transition_signal_col = "simba_median_balance" if "simba_median_balance" in x_array.columns else simba_col
print(f"Using transition signal column: {transition_signal_col}")

base_subset = (
    x_array
    .query("condition == 'deplete' & infusiontype == '45NaCl'")
    .copy()
)

# Keep behavior-derived transition section aligned with the primary Figure 4 subset.
if "rats_to_exclude" in locals():
    base_subset = base_subset.query("id not in @rats_to_exclude").copy()


def fit_single_animal_transition(animal_df, signal_col):
    animal_df = animal_df.sort_values("trial").copy()
    x = animal_df["trial"].to_numpy(dtype=float)
    y = animal_df[signal_col].to_numpy(dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]

    if len(x) < 6:
        return {"ok": False, "reason": "too_few_points", "x0": np.nan, "k": np.nan, "popt": None}

    if np.nanmax(y) - np.nanmin(y) < 1e-8:
        return {"ok": False, "reason": "no_variation", "x0": np.nan, "k": np.nan, "popt": None}

    fit_df = animal_df[["trial", signal_col]].dropna().copy()
    popt = fit_sigmoid(
        fit_df,
        "trial",
        signal_col,
        fit_to_raw_data=True,
        initial_k=-1,
        x0_min=float(np.nanmin(x)),
        x0_max=float(np.nanmax(x)),
        validate_fit=True,
        maxfev=30000,
        return_param_order="legacy",
    )

    if not np.all(np.isfinite(popt)):
        return {"ok": False, "reason": "fit_failed_or_shape_check_failed", "x0": np.nan, "k": np.nan, "popt": None}

    _, x0, k, _ = [float(v) for v in popt]

    # Allow steep k values, but reject positive k.
    if not np.isfinite(k) or k >= 0:
        return {"ok": False, "reason": "k_positive_or_nonfinite", "x0": x0, "k": k, "popt": popt}

    return {"ok": True, "reason": "ok", "x0": x0, "k": k, "popt": popt}


transition_rows = []
for animal in sorted(base_subset["id"].unique()):
    one = base_subset.loc[base_subset["id"] == animal].copy()
    fit = fit_single_animal_transition(one, transition_signal_col)
    transition_rows.append({
        "id": animal,
        "x0_orig": fit["x0"],
        "k": fit["k"],
        "is_valid": bool(fit["ok"]),
        "reason": fit["reason"],
    })

transition_points_df = pd.DataFrame(transition_rows)
valid_transition_points = transition_points_df.query("is_valid == True").dropna(subset=["x0_orig"]).copy()

print("\nTransition fit summary:")
print(transition_points_df[["id", "x0_orig", "k", "is_valid", "reason"]].to_string(index=False))
print(f"\nValid transition points: {len(valid_transition_points)}/{len(transition_points_df)} animals")

if valid_transition_points.empty:
    print("No valid transition points found, cannot realign.")
else:
    valid_ids = sorted(valid_transition_points["id"].unique())
    aligned_input = base_subset.query("id in @valid_ids").copy()

    # make_realigned_trials expects columns 'id' and 'x0_orig'
    realigned_by_behavior = make_realigned_trials(
        aligned_input.copy(),
        valid_transition_points[["id", "x0_orig"]].copy(),
        verbose=False,
    )

    # Prefer fully complete bins; if none exist, gracefully fall back.
    n_required = len(valid_ids)
    realigned_complete = include_only_complete_trials(
        realigned_by_behavior,
        group_col="trial_aligned",
        n_required=n_required,
    )
    alignment_mode = f"complete bins (n={n_required})"

    if realigned_complete.empty:
        # First fallback: bins represented by at least two animals.
        realigned_complete = (
            realigned_by_behavior
            .groupby("trial_aligned", group_keys=False)
            .filter(lambda g: len(g) >= 2)
            .copy()
        )
        alignment_mode = "bins with >=2 animals"

    if realigned_complete.empty or realigned_complete["trial_aligned"].nunique() < 6:
        # Final fallback: keep all bins so the right panel remains interpretable.
        realigned_complete = realigned_by_behavior.copy()
        alignment_mode = "all available bins (unbalanced)"

    def _group_mean_sem(df, group_col, value_col):
        g = df.groupby(group_col)[value_col]
        mean = g.mean().sort_index()
        sem = g.sem().sort_index().fillna(0)
        return mean.index.to_numpy(dtype=float), mean.to_numpy(dtype=float), sem.to_numpy(dtype=float)

    def _fit_group_k(x_vals, y_vals):
        if len(x_vals) < 6:
            return np.nan
        if np.nanmax(y_vals) - np.nanmin(y_vals) < 1e-8:
            return np.nan
        fit_df = pd.DataFrame({"trial_like": x_vals, "signal": y_vals})
        popt = fit_sigmoid(
            fit_df,
            "trial_like",
            "signal",
            fit_to_raw_data=True,
            initial_k=-1,
            x0_min=float(np.nanmin(x_vals)),
            x0_max=float(np.nanmax(x_vals)),
            validate_fit=True,
            maxfev=30000,
            return_param_order="legacy",
        )
        if not np.all(np.isfinite(popt)):
            return np.nan
        return float(popt[2])

    x_do, y_do, se_do = _group_mean_sem(aligned_input, "trial", da_column)
    x_dr, y_dr, se_dr = _group_mean_sem(realigned_complete, "trial_aligned", da_column)

    x_bo, y_bo, se_bo = _group_mean_sem(aligned_input, "trial", transition_signal_col)
    x_br, y_br, se_br = _group_mean_sem(realigned_complete, "trial_aligned", transition_signal_col)

    k_do = _fit_group_k(x_do, y_do)
    k_dr = _fit_group_k(x_dr, y_dr)
    k_bo = _fit_group_k(x_bo, y_bo)
    k_br = _fit_group_k(x_br, y_br)

    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.8), sharey=False)

    # Panel 1: original
    axes[0].plot(x_do, y_do, linestyle="", marker="o", markersize=4.5, color=da_color, markerfacecolor="white")
    axes[0].fill_between(x_do, y_do - se_do, y_do + se_do, color=da_color, alpha=0.12)
    ax0b = axes[0].twinx()
    ax0b.plot(x_bo, y_bo, linestyle="", marker="o", markersize=4.5, color=behav_color, markerfacecolor="white")
    ax0b.fill_between(x_bo, y_bo - se_bo, y_bo + se_bo, color=behav_color, alpha=0.12)

    axes[0].set_title("Original trial order")
    axes[0].set_xlabel("Trial")
    axes[0].set_ylabel("Dopamine AUC", color=da_color)
    ax0b.set_ylabel(f"Behavior ({transition_signal_col})", color=behav_color, rotation=270, labelpad=13)

    # Panel 2: realigned
    axes[1].plot(x_dr, y_dr, linestyle="", marker="o", markersize=4.5, color=da_color, markerfacecolor="white")
    axes[1].fill_between(x_dr, y_dr - se_dr, y_dr + se_dr, color=da_color, alpha=0.12)
    ax1b = axes[1].twinx()
    ax1b.plot(x_br, y_br, linestyle="", marker="o", markersize=4.5, color=behav_color, markerfacecolor="white")
    ax1b.fill_between(x_br, y_br - se_br, y_br + se_br, color=behav_color, alpha=0.12)
    axes[1].axvline(0, color="0.5", linestyle="--", linewidth=1)

    axes[1].set_title("Realigned to behavior-derived transition")
    axes[1].set_xlabel("Trial relative to transition")
    axes[1].set_ylabel("Dopamine AUC", color=da_color)
    ax1b.set_ylabel(f"Behavior ({transition_signal_col})", color=behav_color, rotation=270, labelpad=13)

    for ax in [axes[0], axes[1]]:
        ax.tick_params(axis="y", labelcolor=da_color)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    for ax in [ax0b, ax1b]:
        ax.tick_params(axis="y", labelcolor=behav_color)
        ax.spines["top"].set_visible(False)

    axes[0].text(0.02, 0.98, f"DA k={k_do:.3f}\nBeh k={k_bo:.3f}", transform=axes[0].transAxes, va="top", fontsize=9)
    axes[1].text(0.02, 0.98, f"DA k={k_dr:.3f}\nBeh k={k_br:.3f}", transform=axes[1].transAxes, va="top", fontsize=9)

    fig.suptitle(
        f"Realignment using behavior-derived transition points ({len(valid_ids)} animals)",
        y=1.03,
    )
    fig.tight_layout()
    plt.show()

    print("\nRealignment summary:")
    print(f"  Included animals: {valid_ids}")
    print(f"  Alignment mode: {alignment_mode}")
    print(f"  Original rows: {len(aligned_input)}")
    print(f"  Realigned rows used: {len(realigned_complete)}")
    print(f"  Dopamine k: original={k_do:.4f}, realigned={k_dr:.4f}, delta={k_dr - k_do:.4f}")
    print(f"  Behavior k: original={k_bo:.4f}, realigned={k_br:.4f}, delta={k_br - k_bo:.4f}")
