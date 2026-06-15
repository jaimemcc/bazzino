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
# # Paper Figures: Figure 3 - Lagged Correlation Analysis
#
# Cross-correlation analysis between dopamine (photometry AUC) and behavior (time moving) across different time lags. Determines optimal lag for each animal where correlation is strongest.

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
from scipy.stats import pearsonr
from scipy.ndimage import gaussian_filter1d

from trompy import save_figure_atomic

from figure_config import (
    configure_matplotlib, COLORS,
    DATAFOLDER, RESULTSFOLDER, FIGSFOLDER,
    SAVE_FIGS
)
from figure_plotting import plot_lag_histogram, plot_lag_peak_sharpness

from utils import recalculate_time_moving

SAVE_FIGS = True

# Configure matplotlib
configure_matplotlib()
colors = COLORS

SIMBA_METRIC = "median" # other options...

# %% [markdown]
# ## Load data

# %%
# Load assembled data from figure_4 preprocessing
assembled_data_path = DATAFOLDER / "assembled_data.pickle"

with open(assembled_data_path, "rb") as f:
    data = dill.load(f)

x_array = data["x_array"]
snips_movement = data["snips_movement"]

if SIMBA_METRIC == "zscore":
    snips_simba = data["snips_simba_zscore"]
    simba_col = "simba_mean_zscore_auc"
elif SIMBA_METRIC == "median":
    snips_simba = data["snips_simba_baseline"]
    simba_col = "simba_median_balance"
else:
    raise ValueError(f"Unknown SIMBA_METRIC: {SIMBA_METRIC}")

da_col = "auc_snips"

#x_array = recalculate_time_moving(x_array, snips_movement, threshold=3)

# Get the deplete + 45NaCl animals
subset_condition = (
    x_array
    .query("condition == 'deplete' & infusiontype == '45NaCl'")
    .reset_index(drop=True)
    .sort_values(['id', 'trial'])
)

subset_aligned = (
    subset_condition
    .dropna(subset=['trial_aligned'])
    .reset_index(drop=True)
    .sort_values(['id', 'trial'])
)


# %% [markdown]
# ## Helper Functions

# %%
def include_only_complete_trials(df, group_col="trial_aligned", n_required=8):
    """Filter df to include only groups with at least n_required entries."""
    return (
        df
        .groupby(group_col, group_keys=False)
        .filter(lambda g: len(g) >= n_required)
        .copy()
    )


def compute_lagged_pearson(dopamine, behavior, lag_range=(-5, 5), min_overlap=10):
    """Compute Pearson r for each lag in the requested window."""
    n_trials = len(dopamine)
    lags = np.arange(-(n_trials - 1), n_trials)

    valid_lag_mask = (lags >= lag_range[0]) & (lags <= lag_range[1])
    valid_lags = lags[valid_lag_mask]

    rows = []
    for lag in valid_lags:
        overlap = len(dopamine) - abs(lag)
        if overlap < min_overlap:
            rows.append({"lag_trials": lag, "pearson_r": np.nan, "n_samples": overlap})
            continue

        if lag < 0:
            dop_lag = dopamine[:lag]
            behav_lag = behavior[-lag:]
        elif lag > 0:
            dop_lag = dopamine[lag:]
            behav_lag = behavior[:-lag]
        else:
            dop_lag = dopamine
            behav_lag = behavior

        try:
            r, _ = pearsonr(dop_lag, behav_lag)
        except:
            r = np.nan

        rows.append({"lag_trials": lag, "pearson_r": r, "n_samples": len(dop_lag)})

    return pd.DataFrame(rows)


def select_optimal_lag(dopamine, behavior, lag_df, direction="both"):
    """Select lag using an r-sign constraint, allowing both positive and negative lags."""
    direction = str(direction).lower()
    if direction not in {"both", "positive", "negative"}:
        raise ValueError("direction must be one of: 'both', 'positive', 'negative'")

    valid_df = lag_df[lag_df["pearson_r"].notna()].copy()

    if direction == "positive":
        valid_df = valid_df[valid_df["pearson_r"] > 0]
    elif direction == "negative":
        valid_df = valid_df[valid_df["pearson_r"] < 0]

    if valid_df.empty:
        return {
            "lag_trials": np.nan,
            "pearson_r": np.nan,
            "p_value": np.nan,
            "n_samples": np.nan,
        }

    if direction == "positive":
        best_idx = valid_df["pearson_r"].idxmax()
    elif direction == "negative":
        best_idx = valid_df["pearson_r"].idxmin()
    else:
        best_idx = valid_df["pearson_r"].abs().idxmax()

    best_lag = int(lag_df.loc[best_idx, "lag_trials"])
    best_r = lag_df.loc[best_idx, "pearson_r"]

    if best_lag < 0:
        dop_lag = dopamine[:best_lag]
        behav_lag = behavior[-best_lag:]
    elif best_lag > 0:
        dop_lag = dopamine[best_lag:]
        behav_lag = behavior[:-best_lag]
    else:
        dop_lag = dopamine
        behav_lag = behavior

    try:
        _, p_value = pearsonr(dop_lag, behav_lag)
    except:
        p_value = np.nan

    return {
        "lag_trials": best_lag,
        "pearson_r": best_r,
        "p_value": p_value,
        "n_samples": len(dop_lag),
    }


def compute_optimal_lag_correlation(dopamine, behavior, lag_range=(-5, 5), min_overlap=10, direction="both"):
    """Compute optimal lag correlation between two signals using Pearson r."""
    lag_df = compute_lagged_pearson(
        dopamine,
        behavior,
        lag_range=lag_range,
        min_overlap=min_overlap,
    )
    return select_optimal_lag(dopamine, behavior, lag_df, direction=direction)


# %%
def compute_lag_peak_sharpness(
    subset_condition,
    animals,
    lag_range=(-5, 5),
    min_overlap=10,
    max_offset=5,
    direction="both",
    trial_col="trial",
    dopamine_col="auc_snips",
    behavior_col="mean_simba",
):
    """Compute peak sharpness as r_best - r_at_offset for each animal."""
    lag_difference_data = []

    for animal_name in sorted(animals):
        animal_data = (
            subset_condition
            .query("id == @animal_name")
            .sort_values(trial_col)
        )

        dopamine_auc = animal_data[dopamine_col].to_numpy()
        behavior_auc = animal_data[behavior_col].to_numpy()

        dopamine_norm = (
            (dopamine_auc - dopamine_auc.min())
            / (dopamine_auc.max() - dopamine_auc.min() + 1e-10)
        )
        behavior_norm = (
            (behavior_auc - behavior_auc.min())
            / (behavior_auc.max() - behavior_auc.min() + 1e-10)
        )

        lag_df = compute_lagged_pearson(
            dopamine_norm,
            behavior_norm,
            lag_range=lag_range,
            min_overlap=min_overlap,
        )

        best_lag_info = select_optimal_lag(dopamine_norm, behavior_norm, lag_df, direction=direction)
        if not np.isfinite(best_lag_info["lag_trials"]):
            continue

        best_lag = int(best_lag_info["lag_trials"])
        best_r = best_lag_info["pearson_r"]

        lag_difference_data.append({
            "animal": animal_name,
            "offset": 0,
            "lag_target": best_lag,
            "r_difference": 0.0,
            "best_lag": best_lag,
            "best_r": best_r,
        })

        lag_df_dict = dict(zip(lag_df["lag_trials"], lag_df["pearson_r"]))
        for offset in range(1, max_offset + 1):
            lag_minus = best_lag - offset
            lag_plus = best_lag + offset

            r_minus = lag_df_dict.get(lag_minus, np.nan)
            r_plus = lag_df_dict.get(lag_plus, np.nan)

            if np.isfinite(r_minus):
                lag_difference_data.append({
                    "animal": animal_name,
                    "offset": -offset,
                    "lag_target": lag_minus,
                    "r_difference": best_r - r_minus,
                    "best_lag": best_lag,
                    "best_r": best_r,
                })

            if np.isfinite(r_plus):
                lag_difference_data.append({
                    "animal": animal_name,
                    "offset": offset,
                    "lag_target": lag_plus,
                    "r_difference": best_r - r_plus,
                    "best_lag": best_lag,
                    "best_r": best_r,
                })

    lag_diff_df = pd.DataFrame(lag_difference_data)
    lag_diff_summary = (
        lag_diff_df.groupby("offset")
        .agg({"r_difference": ["mean", "sem", "std", "count"]})
        .reset_index()
    )
    lag_diff_summary.columns = ["offset", "mean_diff", "sem_diff", "std_diff", "n"]

    return lag_diff_df, lag_diff_summary


def print_lag_peak_sharpness_summary(lag_diff_summary, offsets=(1, 2, 5)):
    """Print summary drop values at selected offset magnitudes."""
    print("\nPeak Sharpness Summary:")
    for offset in offsets:
        row = lag_diff_summary[lag_diff_summary["offset"].abs() == offset]
        print(
            f"  Mean drop at ±{offset} lag: "
            f"{row['mean_diff'].mean():.4f} ± {row['sem_diff'].mean():.4f}"
        )


# %% [markdown]
# ## Figure 3A: Cross-Correlation Analysis (Full Data, -10 to +10 trials)
#
# Compute cross-correlation across all trials and animals to see correlation patterns without transition alignment.

# %%
def get_lags_all_animals(x_array, condition, infusiontype, lag_range=(-10, 10), direction="both", smoothing=False):
    """Compute optimal lag for each animal with an r-sign constraint."""
    cross_corr_results = []

    subset_condition = x_array.query("condition == @condition & infusiontype == @infusiontype")
    animals = subset_condition.id.unique()

    for animal in animals:
        animal_data = (
            subset_condition
            .query("id == @animal")
            .sort_values('trial')
        )

        if len(animal_data) < 5:
            print(f"Warning: {animal} has {len(animal_data)} trials in window - skipping")
            continue

        dopamine_auc = animal_data[da_col].values
        behavior_auc = animal_data[simba_col].values
        
        if smoothing:
            dopamine_auc = gaussian_filter1d(dopamine_auc, sigma=1)
            behavior_auc = gaussian_filter1d(behavior_auc, sigma=1)

        dopamine_norm = (dopamine_auc - dopamine_auc.min()) / (dopamine_auc.max() - dopamine_auc.min() + 1e-10)
        behavior_norm = (behavior_auc - behavior_auc.min()) / (behavior_auc.max() - behavior_auc.min() + 1e-10)

        result = compute_optimal_lag_correlation(
            dopamine_norm,
            behavior_norm,
            lag_range=lag_range,
            direction=direction,
        )
        result['animal'] = animal
        result['n_trials'] = len(animal_data)
        result['analysis_type'] = f'{condition}_{infusiontype}'

        cross_corr_results.append(result)

    return pd.DataFrame(cross_corr_results)

# Compute optimal lag for each animal (full data, -10 to +10 trials)
lag_range = (-5, 5)

direction = "positive"  # r-value constraint: "positive", "negative", "both"

cross_corr_deplete_45NaCl = get_lags_all_animals(
    x_array,
    lag_range=lag_range,
    condition='deplete',
    infusiontype='45NaCl',
    direction=direction,
)

print("\nFull Data Lagged Correlation (Dopamine AUC vs Time Moving)")
print("Window: all trials (lags can be negative or positive)")
print(f"r-value constraint: {direction}")
print(cross_corr_deplete_45NaCl[['animal', 'lag_trials', 'pearson_r', 'p_value', 'n_trials']].to_string(index=False))
print(f"\nLag Summary (in trials):")
print(f"  Mean lag: {cross_corr_deplete_45NaCl['lag_trials'].mean():.2f} trials")
print(f"  Median lag: {cross_corr_deplete_45NaCl['lag_trials'].median():.2f} trials")
print(f"  Std: {cross_corr_deplete_45NaCl['lag_trials'].std():.2f} trials")
print(f"\nMean Pearson correlation at optimal lag: {cross_corr_deplete_45NaCl['pearson_r'].mean():.3f}")

# %%
# Stats: compare optimal lags to zero for all groups and between groups
from itertools import combinations
from scipy.stats import f_oneway, ttest_1samp, ttest_ind, wilcoxon

# Ensure all groups are available in case this cell is run independently.
if "cross_corr_deplete_45NaCl" not in globals():
    cross_corr_deplete_45NaCl = get_lags_all_animals(
        x_array, lag_range=lag_range, condition="deplete", infusiontype="45NaCl", direction=direction
    )
if "cross_corr_deplete_10NaCl" not in globals():
    cross_corr_deplete_10NaCl = get_lags_all_animals(
        x_array, lag_range=lag_range, condition="deplete", infusiontype="10NaCl", direction=direction
    )
if "cross_corr_replete_45NaCl" not in globals():
    cross_corr_replete_45NaCl = get_lags_all_animals(
        x_array, lag_range=lag_range, condition="replete", infusiontype="45NaCl", direction=direction
    )
if "cross_corr_replete_10NaCl" not in globals():
    cross_corr_replete_10NaCl = get_lags_all_animals(
        x_array, lag_range=lag_range, condition="replete", infusiontype="10NaCl", direction=direction
    )

group_lag_data = {
    "replete_10NaCl": cross_corr_replete_10NaCl,
    "replete_45NaCl": cross_corr_replete_45NaCl,
    "deplete_10NaCl": cross_corr_deplete_10NaCl,
    "deplete_45NaCl": cross_corr_deplete_45NaCl,
}

# --- One-sample tests: each group against zero lag ---
rows = []
for group_name, group_df in group_lag_data.items():
    lags = group_df["lag_trials"].dropna().to_numpy(dtype=float)
    n = len(lags)

    t_stat, t_p = (np.nan, np.nan)
    w_stat, w_p = (np.nan, np.nan)

    if n >= 2:
        t_stat, t_p = ttest_1samp(lags, popmean=0.0, alternative="two-sided")

        # If all lags are exactly zero, Wilcoxon is undefined after dropping zero diffs.
        if np.allclose(lags, 0):
            w_stat, w_p = (0.0, 1.0)
        else:
            try:
                w_stat, w_p = wilcoxon(lags, zero_method="wilcox", alternative="two-sided", mode="auto")
            except ValueError:
                w_stat, w_p = (np.nan, np.nan)

    rows.append({
        "group": group_name,
        "n_animals": n,
        "mean_lag": np.mean(lags) if n else np.nan,
        "median_lag": np.median(lags) if n else np.nan,
        "t_stat": t_stat,
        "t_p": t_p,
        "wilcoxon_stat": w_stat,
        "wilcoxon_p": w_p,
    })

lag_vs_zero_stats = pd.DataFrame(rows)

# Bonferroni correction across 4 groups for each test family.
lag_vs_zero_stats["t_p_bonf"] = np.minimum(lag_vs_zero_stats["t_p"] * len(lag_vs_zero_stats), 1.0)
lag_vs_zero_stats["wilcoxon_p_bonf"] = np.minimum(
    lag_vs_zero_stats["wilcoxon_p"] * len(lag_vs_zero_stats), 1.0
)

print("Optimal lag vs zero (two-sided one-sample tests)")
display(
    lag_vs_zero_stats[[
        "group",
        "n_animals",
        "mean_lag",
        "median_lag",
        "t_stat",
        "t_p",
        "t_p_bonf",
        "wilcoxon_stat",
        "wilcoxon_p",
        "wilcoxon_p_bonf",
    ]].round(4)
)

# --- Between-group test: one-way ANOVA ---
lag_group_rows = []
for group_name, group_df in group_lag_data.items():
    for lag in group_df["lag_trials"].dropna().to_numpy(dtype=float):
        lag_group_rows.append({"group": group_name, "lag_trials": lag})

lag_group_df = pd.DataFrame(lag_group_rows)
group_order = list(group_lag_data.keys())
group_arrays = [lag_group_df.loc[lag_group_df["group"] == g, "lag_trials"].to_numpy() for g in group_order]

if all(len(arr) >= 2 for arr in group_arrays):
    f_stat, anova_p = f_oneway(*group_arrays)
    k = len(group_arrays)
    n_total = int(np.sum([len(arr) for arr in group_arrays]))
    df_between = k - 1
    df_within = n_total - k
    eta_sq = (f_stat * df_between) / (f_stat * df_between + df_within) if np.isfinite(f_stat) else np.nan

    lag_group_anova = pd.DataFrame([
        {
            "F": f_stat,
            "p": anova_p,
            "df_between": df_between,
            "df_within": df_within,
            "eta_sq": eta_sq,
        }
    ])

    print("\nOne-way ANOVA on optimal lag by group")
    display(lag_group_anova.round(4))

    if anova_p < 0.05:
        print("ANOVA is significant (p < 0.05); running post hoc tests.")
        try:
            from statsmodels.stats.multicomp import pairwise_tukeyhsd

            tukey = pairwise_tukeyhsd(
                endog=lag_group_df["lag_trials"],
                groups=lag_group_df["group"],
                alpha=0.05,
            )
            lag_group_posthoc = pd.DataFrame(
                tukey.summary().data[1:],
                columns=tukey.summary().data[0],
            )
            for col in ["meandiff", "p-adj", "lower", "upper"]:
                lag_group_posthoc[col] = pd.to_numeric(lag_group_posthoc[col], errors="coerce")

            print("Tukey HSD post hoc results")
            display(lag_group_posthoc.round(4))
        except Exception as exc:
            print(f"Tukey HSD unavailable ({exc}); using Welch t-tests with Bonferroni correction.")
            posthoc_rows = []
            for g1, g2 in combinations(group_order, 2):
                arr1 = lag_group_df.loc[lag_group_df["group"] == g1, "lag_trials"].to_numpy()
                arr2 = lag_group_df.loc[lag_group_df["group"] == g2, "lag_trials"].to_numpy()
                t_stat_pair, p_pair = ttest_ind(arr1, arr2, equal_var=False)
                posthoc_rows.append({
                    "group1": g1,
                    "group2": g2,
                    "mean_diff": np.mean(arr1) - np.mean(arr2),
                    "t_stat": t_stat_pair,
                    "p_uncorrected": p_pair,
                })

            lag_group_posthoc = pd.DataFrame(posthoc_rows)
            lag_group_posthoc["p_bonf"] = np.minimum(
                lag_group_posthoc["p_uncorrected"] * len(lag_group_posthoc),
                1.0,
            )
            lag_group_posthoc["significant_bonf"] = lag_group_posthoc["p_bonf"] < 0.05

            print("Pairwise Welch t-tests post hoc results")
            display(lag_group_posthoc.round(4))
    else:
        print("ANOVA is not significant (p >= 0.05); post hoc tests not run.")
else:
    print("Not enough observations per group for ANOVA (need at least n=2 in each group).")

# %%
# Plot histogram of lags (full data)
f, ax = plt.subplots(figsize=(1.6, 1.8),
                     gridspec_kw={"left": 0.25, "right": 0.95, "top": 0.93, "bottom": 0.25})

plot_lag_histogram(
    ax,
    cross_corr_deplete_45NaCl['lag_trials'],
    lag_range=lag_range,
)

ax.set_yticks([0, 1, 2, 3])
ax.set_ylabel("No. animals")
ax.set_xlabel("Optimal lag")

if SAVE_FIGS:
    save_figure_atomic(f, "fig3_lagged_correlation_full", FIGSFOLDER)


# %%

cross_corr_deplete_45NaCl = get_lags_all_animals(
    x_array,
    lag_range=lag_range,
    condition='deplete',
    infusiontype='45NaCl',
    direction=direction
)

cross_corr_deplete_10NaCl = get_lags_all_animals(
    x_array,
    lag_range=lag_range,
    condition='deplete',
    infusiontype='10NaCl',
    direction=direction
)

cross_corr_replete_45NaCl = get_lags_all_animals(
    x_array,
    lag_range=lag_range,
    condition='replete',
    direction=direction,
    infusiontype='45NaCl'
)

cross_corr_replete_10NaCl = get_lags_all_animals(
    x_array,
    lag_range=lag_range,
    condition='replete',
    infusiontype='10NaCl',
    direction=direction
)

cross_corr_not_dep45 = pd.concat([
    cross_corr_deplete_10NaCl,
    cross_corr_replete_45NaCl,
    cross_corr_replete_10NaCl,])

# %%
f, ax = plt.subplots(figsize=(1.8, 1.8),
                     gridspec_kw={"left": 0.35, "right": 0.95, "top": 0.93, "bottom": 0.25})

# rng = np.random.default_rng(7)
np.random.seed(123)
jitter_width = 0.08

for idx, df in enumerate([
    cross_corr_replete_10NaCl,
    cross_corr_replete_45NaCl,
    cross_corr_deplete_10NaCl,
    cross_corr_deplete_45NaCl,
]):
    x_jitter = idx + np.random.uniform(-jitter_width, jitter_width, size=len(df))
    ax.scatter(x_jitter, df.pearson_r.values, color=colors[idx], alpha=0.7, edgecolor='k', linewidth=0.5, s=30)
    ax.bar(idx, df.pearson_r.mean(), color=colors[idx], alpha=0.5)

ax.set_xticks((0, 1, 2, 3))
ax.set_xticklabels(["","","",""])

ax.set_ylim(0, 0.8)

ax.set_ylabel("Pearson r at optimal lag")

sns.despine(ax=ax)

if SAVE_FIGS:
    save_figure_atomic(f, "fig3_optimal_lag_correlation_all_conds", FIGSFOLDER)

# %%
# Stats: test whether Pearson r differs across groups
from itertools import combinations
from scipy.stats import f_oneway, ttest_ind

r_group_data = {
    "replete_10NaCl": cross_corr_replete_10NaCl["pearson_r"].dropna().to_numpy(dtype=float),
    "replete_45NaCl": cross_corr_replete_45NaCl["pearson_r"].dropna().to_numpy(dtype=float),
    "deplete_10NaCl": cross_corr_deplete_10NaCl["pearson_r"].dropna().to_numpy(dtype=float),
    "deplete_45NaCl": cross_corr_deplete_45NaCl["pearson_r"].dropna().to_numpy(dtype=float),
}

group_order = list(r_group_data.keys())
group_arrays = [r_group_data[g] for g in group_order]

if all(len(arr) >= 2 for arr in group_arrays):
    f_stat, anova_p = f_oneway(*group_arrays)
    k = len(group_arrays)
    n_total = int(np.sum([len(arr) for arr in group_arrays]))
    df_between = k - 1
    df_within = n_total - k
    eta_sq = (f_stat * df_between) / (f_stat * df_between + df_within) if np.isfinite(f_stat) else np.nan

    pearson_group_anova = pd.DataFrame([
        {
            "F": f_stat,
            "p": anova_p,
            "df_between": df_between,
            "df_within": df_within,
            "eta_sq": eta_sq,
        }
    ])

    print("One-way ANOVA on Pearson r by group")
    display(pearson_group_anova.round(4))

    if anova_p < 0.05:
        print("ANOVA is significant (p < 0.05); running post hoc tests.")
        try:
            from statsmodels.stats.multicomp import pairwise_tukeyhsd

            r_group_df = pd.concat(
                [pd.DataFrame({"group": g, "pearson_r": vals}) for g, vals in r_group_data.items()],
                ignore_index=True,
            )
            tukey = pairwise_tukeyhsd(
                endog=r_group_df["pearson_r"],
                groups=r_group_df["group"],
                alpha=0.05,
            )
            pearson_group_posthoc = pd.DataFrame(
                tukey.summary().data[1:],
                columns=tukey.summary().data[0],
            )
            for col in ["meandiff", "p-adj", "lower", "upper"]:
                pearson_group_posthoc[col] = pd.to_numeric(pearson_group_posthoc[col], errors="coerce")

            print("Tukey HSD post hoc results")
            display(pearson_group_posthoc.round(4))
        except Exception as exc:
            print(f"Tukey HSD unavailable ({exc}); using Welch t-tests with Bonferroni correction.")
            posthoc_rows = []
            for g1, g2 in combinations(group_order, 2):
                arr1 = r_group_data[g1]
                arr2 = r_group_data[g2]
                t_stat_pair, p_pair = ttest_ind(arr1, arr2, equal_var=False)
                posthoc_rows.append({
                    "group1": g1,
                    "group2": g2,
                    "mean_diff": np.mean(arr1) - np.mean(arr2),
                    "t_stat": t_stat_pair,
                    "p_uncorrected": p_pair,
                })

            pearson_group_posthoc = pd.DataFrame(posthoc_rows)
            pearson_group_posthoc["p_bonf"] = np.minimum(
                pearson_group_posthoc["p_uncorrected"] * len(pearson_group_posthoc),
                1.0,
            )
            pearson_group_posthoc["significant_bonf"] = pearson_group_posthoc["p_bonf"] < 0.05

            print("Pairwise Welch t-tests post hoc results")
            display(pearson_group_posthoc.round(4))
    else:
        print("ANOVA is not significant (p >= 0.05); post hoc tests not run.")
else:
    print("Not enough observations per group for ANOVA (need at least n=2 in each group).")

# %%
bw=0.75

f, ax = plt.subplots(figsize=(1.8,1.8),
                     gridspec_kw={"left": 0.35, "right": 0.95, "top": 0.93, "bottom": 0.25})

sns.kdeplot(cross_corr_not_dep45['lag_trials'],
            fill=True, ax=ax, color='gray',
            bw_adjust=bw,
            clip=[cross_corr_not_dep45['lag_trials'].min(),
                  cross_corr_not_dep45['lag_trials'].max()],
            common_norm=True,
            label='Other conditions')

sns.kdeplot(cross_corr_deplete_45NaCl['lag_trials'],
            fill=True, ax=ax, color=colors[3],
            bw_adjust=bw,
            clip=[cross_corr_deplete_45NaCl['lag_trials'].min(),
                  cross_corr_deplete_45NaCl['lag_trials'].max()],
            common_norm=True,
            label='Deplete 45NaCl')

ax.set_xlim(lag_range)
ax.set_xlabel("Optimal lag")
ax.set_yticks([0, 0.1])
ax.axvline(0, color='k', linestyle='--', alpha=0.7)
sns.despine(ax=ax, offset=0)


if SAVE_FIGS:
    save_figure_atomic(f, "fig3_lagged_correlation_kde", FIGSFOLDER)


# %%
# Peak sharpness analysis: difference between best lag and adjacent lags
animals_full = x_array.query("condition == 'deplete' & infusiontype == '45NaCl'").id.unique()

lag_range = (-10, 10)
min_overlap = 10

lag_diff_df, lag_diff_summary = compute_lag_peak_sharpness(
    subset_condition=subset_condition,
    animals=animals_full,
    lag_range=lag_range,
    min_overlap=min_overlap,
    max_offset=5,
    trial_col="trial",
    dopamine_col="auc_snips",
    behavior_col=simba_col,
)

f, ax = plt.subplots(figsize=(1.8, 1.8),
                     gridspec_kw={"left": 0.35, "right": 0.95, "top": 0.93, "bottom": 0.25})
plot_lag_peak_sharpness(
    ax=ax,
    lag_diff_df=lag_diff_df,
    lag_diff_summary=lag_diff_summary,
    color=colors[2],
)

ax.axhline(0, color="black", linestyle="--", linewidth=1, alpha=0.5)
ax.axvline(0, color="black", linestyle="--", linewidth=1, alpha=0.5)
ax.set_xlabel("Offset from best lag")
ax.set_ylabel("Diff. in Pearson r")
ax.set_ylim([-0.05, 0.25])
sns.despine(ax=ax)

if SAVE_FIGS:
    save_figure_atomic(f, "fig3_lag_sharpness", FIGSFOLDER)

print_lag_peak_sharpness_summary(lag_diff_summary, offsets=(1, 2, 5))

# %%
# probability density of lags for each condition
bw=0.75

f, ax = plt.subplots(ncols=4, figsize=(4, 3), sharey=True,
                     gridspec_kw={'width_ratios': [1, 1, 1, 1], 'wspace': 0.3})
sns.kdeplot(cross_corr_replete_10NaCl['lag_trials'], color=colors[0], label='replete 10NaCl', fill=True, alpha=0.5, ax=ax[0], bw_adjust=bw, clip=lag_range)
sns.kdeplot(cross_corr_replete_45NaCl['lag_trials'], color=colors[1], label='replete 45NaCl', fill=True, alpha=0.5, ax=ax[1], bw_adjust=bw, clip=lag_range)
sns.kdeplot(cross_corr_deplete_10NaCl['lag_trials'], color=colors[2], label='deplete 10NaCl', fill=True, alpha=0.5, ax=ax[2], bw_adjust=bw, clip=lag_range)
sns.kdeplot(cross_corr_deplete_45NaCl['lag_trials'], color=colors[3], label='deplete 45NaCl', fill=True, alpha=0.5, ax=ax[3], bw_adjust=bw, clip=lag_range)

for axis in ax:
    sns.despine(ax=axis, offset=5, )

    axis.axvline(0, color='gray', linestyle=':', linewidth=1, alpha=0.5)
    axis.set_xlabel('Lag')
    axis.set_xlim(lag_range)


# %% [markdown]
# ## Individual Animals: Lag Profiles

# %%
# Single animal lag plot (match optimal-lag calculations)
animal_to_plot = 'PB48'
print(animal_to_plot)

lag_range = (-10, 10)
min_overlap = 10
direction = "both"  # r-value constraint: "positive", "negative", "both"

animal_data = (
    subset_condition
    .query("id == @animal_to_plot")
    .sort_values("trial")
)

dopamine_auc = animal_data["auc_snips"].to_numpy()
time_moving = animal_data["time_moving"].to_numpy()

dopamine_norm = (dopamine_auc - dopamine_auc.min()) / (dopamine_auc.max() - dopamine_auc.min() + 1e-10)
behavior_norm = (time_moving - time_moving.min()) / (time_moving.max() - time_moving.min() + 1e-10)

lag_r_df = compute_lagged_pearson(
    dopamine_norm,
    behavior_norm,
    lag_range=lag_range,
    min_overlap=min_overlap,
)

best = select_optimal_lag(dopamine_norm, behavior_norm, lag_r_df, direction=direction)
best_lag = best["lag_trials"]

fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.5), gridspec_kw={"wspace": 0.35})

axes[0].plot(lag_r_df["lag_trials"], lag_r_df["pearson_r"], marker="o")
axes[0].axhline(0, color="0.5", linewidth=1)
axes[0].axvline(0, color="0.5", linewidth=1, linestyle="--")
if np.isfinite(best_lag):
    axes[0].axvline(best_lag, color="0.25", linewidth=1.5, linestyle="--",
                    label=f"Best lag: {int(best_lag)}")
    axes[0].legend(frameon=False)
axes[0].set_xlabel("Lag (trials)")
axes[0].set_ylabel("Pearson r")
axes[0].set_title(f"Lagged correlation for {animal_to_plot} (r: {direction})")

trial_idx = np.arange(len(dopamine_norm))
axes[1].plot(trial_idx, dopamine_norm, label="Dopamine (norm)", color=colors[0])
axes[1].plot(trial_idx, behavior_norm, label="Time moving (norm)", color=colors[1])
axes[1].set_xlabel("Trial")
axes[1].set_ylabel("Normalized value")
axes[1].set_title("Normalized signals")
axes[1].legend(frameon=False)

if np.isfinite(best_lag):
    if best_lag < 0:
        dop_lag = dopamine_norm[:best_lag]
        behav_lag = behavior_norm[-best_lag:]
    elif best_lag > 0:
        dop_lag = dopamine_norm[best_lag:]
        behav_lag = behavior_norm[:-best_lag]
    else:
        dop_lag = dopamine_norm
        behav_lag = behavior_norm

    aligned_idx = np.arange(len(dop_lag))
    axes[2].plot(aligned_idx, dop_lag, label=f"Dopamine (lag {int(best_lag)})", color=colors[0])
    axes[2].plot(aligned_idx, behav_lag, label="Time moving", color=colors[1])
    axes[2].set_xlabel("Aligned trial index")
    axes[2].set_ylabel("Normalized value")
    axes[2].set_title("Signals aligned by best lag")
    axes[2].legend(frameon=False)
else:
    axes[2].axis("off")
    axes[2].set_title("No valid lag")

sns.despine(fig=fig)
plt.tight_layout()
plt.show()

# lag_r_df

# %%
# Multi-animal overlay: lagged Pearson r for all animals + mean
lag_range = (-10, 10)
min_overlap = 10

all_lag_dfs = []
for animal_name in sorted(animals_full):
    animal_data = (
        subset_condition
        .query("id == @animal_name")
        .sort_values("trial")
    )

    dopamine_auc = animal_data["auc_snips"].to_numpy()
    time_moving = animal_data["time_moving"].to_numpy()

    dopamine_norm = (dopamine_auc - dopamine_auc.min()) / (dopamine_auc.max() - dopamine_auc.min() + 1e-10)
    behavior_norm = (time_moving - time_moving.min()) / (time_moving.max() - time_moving.min() + 1e-10)

    lag_df = compute_lagged_pearson(
        dopamine_norm,
        behavior_norm,
        lag_range=lag_range,
        min_overlap=min_overlap,
    )
    lag_df["animal"] = animal_name

    all_lag_dfs.append(lag_df)

lag_r_all_df = pd.concat(all_lag_dfs, ignore_index=True)
lag_r_mean_df = lag_r_all_df.groupby("lag_trials", as_index=False)["pearson_r"].mean()

plt.figure(figsize=(6.5, 4.5))
for animal_name, df_animal in lag_r_all_df.groupby("animal"):
    plt.plot(df_animal["lag_trials"], df_animal["pearson_r"], color="0.75", linewidth=1, alpha=0.8)

plt.plot(
    lag_r_mean_df["lag_trials"],
    lag_r_mean_df["pearson_r"],
    color="black",
    linewidth=2.5,
    marker="o",
    label="Mean across animals",
)

plt.axhline(0, color="0.5", linewidth=1)
plt.axvline(0, color="0.5", linewidth=1, linestyle="--")
plt.xlabel("Lag (trials)")
plt.ylabel("Pearson r")
plt.title("Lagged correlation across all animals")
plt.legend(frameon=False)
sns.despine()
plt.tight_layout()

if SAVE_FIGS:
    save_figure_atomic(plt.gcf(), "fig5_lagged_correlation_overlay", FIGSFOLDER)

plt.show()

lag_r_all_df.head()

# %% [markdown]
# ## Diagnostic: Individual Animal Lag Selection

# %% [markdown]
# ## Grid Search: Optimal Threshold and Lag Exploration
#
# Explore parameter space to find optimal combination of time_moving threshold and lag for maximizing correlation between dopamine and behavior.

# %%
plt.hist(snips_movement.flatten(), bins=range(0, 50, 1))

# %%
# Import the utility function for recalculating time_moving
from utils import recalculate_time_moving

# Load movement snips from the assembled data
snips_movement = data["snips_movement"]

# Define parameter grids
lag_range = (-10, 10)
direction = "both"  # r-value constraint: "positive", "negative", "both"
threshold_range = np.linspace(0.005, 0.10, 20)  # 20 different thresholds from 0.005 to 0.10 (for normalized movement)
threshold_range = np.arange(0.5, 4, 0.25) # in pixels
min_overlap = 10

# for testing
# lag_range = (-2, 2) # for testing
# threshold_range = np.array([0.02, 0.2]) # for testing

print(f"Grid search parameters:")
print(f"  Lag range: {lag_range}")
print(f"  r-value constraint: {direction}")
print(f"  Number of threshold values: {len(threshold_range)}")
print(f"  Threshold range: [{threshold_range.min():.4f}, {threshold_range.max():.4f}]")
print(f"  Minimum overlap: {min_overlap} trials")
print(f"\nMovement snips shape: {snips_movement.shape}")
print(f"Full data shape: {subset_condition.shape}")

# %%
# Perform grid search: vary threshold and lag, compute optimal correlation for each combination
grid_results = []

for threshold in threshold_range:

    # Recalculate time_moving with this threshold
    x_array_with_new_threshold = recalculate_time_moving(
        x_array,
        snips_movement,
        threshold=threshold
    )

    subset_with_new_threshold = (x_array_with_new_threshold
                                 .query("condition == 'deplete' and infusiontype == '45NaCl'")
                                 )

    # For each animal, compute lagged correlations and find best lag
    for animal in animals_full:
        animal_data = (
            subset_with_new_threshold
            .query("id == @animal")
            .reset_index(drop=True)
        )

        if len(animal_data) < min_overlap * 2:
            continue

        # Extract signals
        dopamine_auc = animal_data['auc_snips'].values
        time_moving = animal_data['time_moving'].values

        # Normalize signals
        dopamine_norm = (dopamine_auc - dopamine_auc.min()) / (dopamine_auc.max() - dopamine_auc.min() + 1e-10)
        behavior_norm = (time_moving - time_moving.min()) / (time_moving.max() - time_moving.min() + 1e-10)

        # Compute lagged correlations
        lag_df = compute_lagged_pearson(dopamine_norm, behavior_norm, lag_range=lag_range, min_overlap=min_overlap)

        # Find best lag using r-value sign constraint (lags remain bidirectional)
        if len(lag_df) > 0:
            best_lag_info = select_optimal_lag(dopamine_norm, behavior_norm, lag_df, direction=direction)
            best_lag = best_lag_info['lag_trials']
            best_r = best_lag_info['pearson_r']
            best_p = best_lag_info['p_value']

            grid_results.append({
                'animal': animal,
                'threshold': threshold,
                'best_lag': best_lag,
                'best_r': best_r,
                'best_abs_r': abs(best_r),
                'pvalue': best_p,
                'direction': direction,
            })

# Convert to DataFrame
grid_df = pd.DataFrame(grid_results)

print(f"Grid search complete!")
print(f"  Total combinations explored: {len(threshold_range)} thresholds × {len(animals_full)} animals = {len(grid_df)} results")
print(f"  r-value constraint: {direction}")
print(f"\nBest overall correlation:")
best_overall = grid_df.loc[grid_df['best_abs_r'].idxmax()]
print(f"  Animal: {best_overall['animal']}")
print(f"  Threshold: {best_overall['threshold']:.4f}")
print(f"  Lag: {int(best_overall['best_lag'])} trials")
print(f"  Pearson r: {best_overall['best_r']:.3f}")
print(f"  p-value: {best_overall['pvalue']:.4f}")

# %%
# Create pivot tables for heatmap visualization
# For each (threshold, lag) combination, show the mean |r| across animals

# Get all unique lags that were found as "best" lags (excluding NaN)
all_lags_found = sorted([lag for lag in grid_df['best_lag'].unique() if not np.isnan(lag)])

# Create a summary: for each threshold, what's the distribution of best lags?
threshold_lag_summary = []
for threshold in threshold_range:
    thresh_data = grid_df.query("threshold == @threshold")
    for lag in all_lags_found:
        lag_count = len(thresh_data.query("best_lag == @lag"))
        mean_r = thresh_data.query("best_lag == @lag")['best_abs_r'].mean() if lag_count > 0 else np.nan
        threshold_lag_summary.append({
            'threshold': threshold,
            'lag': int(lag),
            'count': lag_count,
            'mean_abs_r': mean_r
        })

threshold_lag_df = pd.DataFrame(threshold_lag_summary)

# Pivot for heatmap: rows = threshold, columns = lag, values = mean |r|
heatmap_data = threshold_lag_df.pivot(index='threshold', columns='lag', values='mean_abs_r')

print(f"Heatmap shape: {heatmap_data.shape}")
print(f"Threshold range: {len(heatmap_data.index)} values")
print(f"Lag range covered: {heatmap_data.columns.min()} to {heatmap_data.columns.max()} trials")

# %%
# Plot heatmap: mean |r| as a function of threshold and lag
fig, ax = plt.subplots(1, 1, figsize=(6, 5))

# Create heatmap
im = ax.imshow(
    heatmap_data.values,
    aspect='auto',
    cmap='viridis',
    origin='lower',
    interpolation='nearest'
)

# Set axis labels and ticks
ax.set_xlabel('Lag (trials)', fontsize=12)
ax.set_ylabel('Movement Threshold', fontsize=12)
ax.set_title('Mean |Pearson r| across animals\n(Dopamine vs. Time Moving)', fontsize=14, pad=15)

# Set x-axis ticks to show lag values
lag_ticks = np.arange(len(heatmap_data.columns))
lag_labels = [int(x) for x in heatmap_data.columns]
ax.set_xticks(lag_ticks[::2])  # Show every other lag for readability
ax.set_xticklabels(lag_labels[::2])

# Set y-axis ticks to show threshold values
threshold_ticks = np.arange(len(heatmap_data.index))
threshold_labels = [f'{x:.3f}' for x in heatmap_data.index]
ax.set_yticks(threshold_ticks[::2])  # Show every other threshold for readability
ax.set_yticklabels(threshold_labels[::2])

# Add colorbar
cbar = plt.colorbar(im, ax=ax)
cbar.set_label('Mean |Pearson r|', rotation=270, labelpad=20, fontsize=11)

# Add grid for better readability
ax.set_xticks(lag_ticks - 0.5, minor=True)
ax.set_yticks(threshold_ticks - 0.5, minor=True)
ax.grid(which='minor', color='w', linestyle='-', linewidth=0.5, alpha=0.3)

plt.tight_layout()
plt.show()

# Find and highlight the global maximum
max_r_idx = np.unravel_index(np.nanargmax(heatmap_data.values), heatmap_data.shape)
max_threshold = heatmap_data.index[max_r_idx[0]]
max_lag = heatmap_data.columns[max_r_idx[1]]
max_r = heatmap_data.values[max_r_idx]

print(f"\nGlobal maximum (mean |r| across animals):")
print(f"  Threshold: {max_threshold:.4f}")
print(f"  Lag: {int(max_lag)} trials")
print(f"  Mean |r|: {max_r:.3f}")

# %%
# Additional visualization: Best correlation vs threshold (aggregated across lags)
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Panel 1: Mean and max |r| across animals for each threshold (taking best lag per animal)
threshold_summary = (
    grid_df
    .groupby('threshold')
    .agg({
        'best_abs_r': ['mean', 'std', 'max'],
        'animal': 'count'
    })
    .reset_index()
)
threshold_summary.columns = ['threshold', 'mean_abs_r', 'std_abs_r', 'max_abs_r', 'n_animals']

ax = axes[0]
ax.plot(threshold_summary['threshold'], threshold_summary['mean_abs_r'], marker='o', label='Mean |r|', linewidth=2)
ax.fill_between(
    threshold_summary['threshold'],
    threshold_summary['mean_abs_r'] - threshold_summary['std_abs_r'],
    threshold_summary['mean_abs_r'] + threshold_summary['std_abs_r'],
    alpha=0.3
)
ax.plot(threshold_summary['threshold'], threshold_summary['max_abs_r'], marker='s', linestyle='--', 
        label='Max |r|', linewidth=1.5, alpha=0.7)
ax.axvline(0.02, color='red', linestyle=':', alpha=0.5, label='Default (0.02)')
ax.set_xlabel('Movement Threshold', fontsize=12)
ax.set_ylabel('|Pearson r|', fontsize=12)
ax.set_title('Correlation vs. Threshold\n(Best lag per animal)', fontsize=13)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

# Panel 2: Distribution of best lags across all thresholds
ax = axes[1]
lag_counts = grid_df['best_lag'].value_counts().sort_index()
ax.bar(lag_counts.index, lag_counts.values, color='steelblue', edgecolor='black', alpha=0.7)
ax.set_xlabel('Best Lag (trials)', fontsize=12)
ax.set_ylabel('Frequency', fontsize=12)
ax.set_title('Distribution of Optimal Lags\n(All thresholds & animals)', fontsize=13)
ax.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.show()

# Print summary statistics
print(f"\nSummary Statistics:")
print(f"  Default threshold (0.02) mean |r|: {threshold_summary.loc[threshold_summary['threshold'].sub(0.02).abs().idxmin(), 'mean_abs_r']:.3f}")
best_threshold_row = threshold_summary.loc[threshold_summary['mean_abs_r'].idxmax()]
print(f"  Optimal threshold: {best_threshold_row['threshold']:.4f}")
print(f"  Optimal threshold mean |r|: {best_threshold_row['mean_abs_r']:.3f}")
print(f"\nMost common best lag: {lag_counts.idxmax()} trials (occurred {lag_counts.max()} times)")

# %%

# %% [markdown]
# ## Per-Animal Optimal Parameters
#
# Identify the optimal threshold and lag for each individual animal based on maximum |Pearson r|.

# %%
# Find optimal threshold and lag for each animal
per_animal_optimal = []

for animal in animals_full:
    animal_grid_data = grid_df.query("animal == @animal")
    
    if len(animal_grid_data) == 0:
        continue
    
    # Find the row with maximum |r| for this animal
    best_idx = animal_grid_data['best_abs_r'].idxmax()
    best_row = animal_grid_data.loc[best_idx]
    
    per_animal_optimal.append({
        'animal': animal,
        'optimal_threshold': best_row['threshold'],
        'optimal_lag': best_row['best_lag'],
        'best_r': best_row['best_r'],
        'best_abs_r': best_row['best_abs_r'],
        'pvalue': best_row['pvalue']
    })

per_animal_df = pd.DataFrame(per_animal_optimal).sort_values('best_abs_r', ascending=False)

print("Optimal parameters for each animal (sorted by |r|):")
print("=" * 85)
print(per_animal_df.to_string(index=False))
print("\n" + "=" * 85)
print(f"\nSummary Statistics:")
print(f"  Mean optimal threshold: {per_animal_df['optimal_threshold'].mean():.4f}")
print(f"  Median optimal threshold: {per_animal_df['optimal_threshold'].median():.4f}")
print(f"  Std optimal threshold: {per_animal_df['optimal_threshold'].std():.4f}")
print(f"  Range: [{per_animal_df['optimal_threshold'].min():.4f}, {per_animal_df['optimal_threshold'].max():.4f}]")
print(f"\n  Mean optimal lag: {per_animal_df['optimal_lag'].mean():.2f} trials")
print(f"  Median optimal lag: {per_animal_df['optimal_lag'].median():.2f} trials")
print(f"  Std optimal lag: {per_animal_df['optimal_lag'].std():.2f} trials")
print(f"\n  Mean best |r|: {per_animal_df['best_abs_r'].mean():.3f}")
print(f"  Median best |r|: {per_animal_df['best_abs_r'].median():.3f}")

# %%
# Visualize per-animal optimal parameters
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Panel 1: Optimal threshold per animal
ax = axes[0, 0]
animal_sorted = per_animal_df.sort_values('optimal_threshold')
y_pos = np.arange(len(animal_sorted))
bars = ax.barh(y_pos, animal_sorted['optimal_threshold'], color='steelblue', alpha=0.7, edgecolor='black')
ax.axvline(0.02, color='red', linestyle=':', linewidth=2, alpha=0.5, label='Default (0.02)')
ax.set_yticks(y_pos)
ax.set_yticklabels(animal_sorted['animal'], fontsize=9)
ax.set_xlabel('Optimal Movement Threshold', fontsize=11)
ax.set_ylabel('Animal', fontsize=11)
ax.set_title('Optimal Threshold per Animal', fontsize=12, fontweight='bold')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3, axis='x')

# Panel 2: Optimal lag per animal
ax = axes[0, 1]
animal_sorted_lag = per_animal_df.sort_values('optimal_lag')
y_pos = np.arange(len(animal_sorted_lag))
colors_lag = ['#d34d4d' if lag > 0 else '#3b6fb6' if lag < 0 else '#9a9a9a' 
              for lag in animal_sorted_lag['optimal_lag']]
bars = ax.barh(y_pos, animal_sorted_lag['optimal_lag'], color=colors_lag, alpha=0.7, edgecolor='black')
ax.axvline(0, color='gray', linestyle='-', linewidth=1, alpha=0.5)
ax.set_yticks(y_pos)
ax.set_yticklabels(animal_sorted_lag['animal'], fontsize=9)
ax.set_xlabel('Optimal Lag (trials)', fontsize=11)
ax.set_ylabel('Animal', fontsize=11)
ax.set_title('Optimal Lag per Animal', fontsize=12, fontweight='bold')
ax.grid(True, alpha=0.3, axis='x')

# Panel 3: Best |r| per animal
ax = axes[1, 0]
animal_sorted_r = per_animal_df.sort_values('best_abs_r')
y_pos = np.arange(len(animal_sorted_r))
bars = ax.barh(y_pos, animal_sorted_r['best_abs_r'], color='forestgreen', alpha=0.7, edgecolor='black')
ax.set_yticks(y_pos)
ax.set_yticklabels(animal_sorted_r['animal'], fontsize=9)
ax.set_xlabel('Best |Pearson r|', fontsize=11)
ax.set_ylabel('Animal', fontsize=11)
ax.set_title('Maximum Correlation per Animal', fontsize=12, fontweight='bold')
ax.grid(True, alpha=0.3, axis='x')

# Panel 4: Scatter plot - threshold vs lag colored by |r|
ax = axes[1, 1]
scatter = ax.scatter(
    per_animal_df['optimal_threshold'],
    per_animal_df['optimal_lag'],
    c=per_animal_df['best_abs_r'],
    s=150,
    cmap='viridis',
    edgecolor='black',
    linewidth=1.5,
    alpha=0.8
)
ax.axvline(0.02, color='red', linestyle=':', linewidth=2, alpha=0.5, label='Default threshold')
ax.axhline(0, color='gray', linestyle=':', linewidth=1, alpha=0.5)
ax.set_xlabel('Optimal Threshold', fontsize=11)
ax.set_ylabel('Optimal Lag (trials)', fontsize=11)
ax.set_title('Optimal Parameters per Animal', fontsize=12, fontweight='bold')
cbar = plt.colorbar(scatter, ax=ax)
cbar.set_label('|Pearson r|', rotation=270, labelpad=15, fontsize=10)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Add animal labels to scatter plot
for idx, row in per_animal_df.iterrows():
    ax.annotate(
        row['animal'],
        (row['optimal_threshold'], row['optimal_lag']),
        xytext=(5, 5),
        textcoords='offset points',
        fontsize=7,
        alpha=0.7
    )

plt.tight_layout()
plt.show()

# %%
# Create individual heatmaps for selected animals showing correlation landscape
# Select top 6 animals by best |r| for visualization
top_animals = per_animal_df.head(6)['animal'].values

fig, axes = plt.subplots(2, 3, figsize=(16, 10))
axes = axes.flatten()

for idx, animal in enumerate(top_animals):
    ax = axes[idx]
    
    # Get all grid search results for this animal
    animal_data = grid_df.query("animal == @animal")
    
    # Create pivot table: threshold x lag -> |r|
    # For each threshold, we have one "best_lag", so need to expand to all lags
    # Actually, we need to recompute all lag correlations for each threshold
    # But that's stored in grid_df - grid_df has one row per (animal, threshold) with best_lag
    
    # For visualization, let's show threshold vs best_lag with color = |r|
    scatter = ax.scatter(
        animal_data['threshold'],
        animal_data['best_lag'],
        c=animal_data['best_abs_r'],
        s=60,
        cmap='viridis',
        edgecolor='black',
        linewidth=0.5,
        alpha=0.8
    )
    
    # Mark the global optimum for this animal
    animal_opt = per_animal_df.query("animal == @animal").iloc[0]
    ax.scatter(
        [animal_opt['optimal_threshold']],
        [animal_opt['optimal_lag']],
        s=200,
        marker='*',
        c='red',
        edgecolor='black',
        linewidth=1.5,
        zorder=10,
        label=f"Optimum: r={animal_opt['best_r']:.3f}"
    )
    
    ax.axvline(0.02, color='red', linestyle=':', linewidth=1, alpha=0.3)
    ax.axhline(0, color='gray', linestyle=':', linewidth=1, alpha=0.3)
    ax.set_xlabel('Threshold', fontsize=10)
    ax.set_ylabel('Best Lag (trials)', fontsize=10)
    ax.set_title(f'{animal}', fontsize=11, fontweight='bold')
    ax.legend(fontsize=8, loc='best')
    ax.grid(True, alpha=0.2)
    
    # Add colorbar for each subplot
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('|r|', rotation=270, labelpad=12, fontsize=9)

plt.tight_layout()
plt.show()

print(f"Showing correlation landscapes for top {len(top_animals)} animals by best |r|")
print("Red star indicates the optimal (threshold, lag) combination for each animal")

# %%
# Compare default threshold (0.02) vs optimal threshold for each animal
comparison_data = []

for animal in animals_full:
    # Get optimal parameters
    opt_row = per_animal_df.query("animal == @animal")
    if len(opt_row) == 0:
        continue
    opt_row = opt_row.iloc[0]
    
    # Get performance at default threshold (0.02)
    default_data = grid_df.query("animal == @animal and threshold == 0.02")
    
    if len(default_data) > 0:
        default_r = default_data.iloc[0]['best_abs_r']
        default_lag = default_data.iloc[0]['best_lag']
    else:
        # Find closest threshold to 0.02
        animal_thresholds = grid_df.query("animal == @animal")
        closest_idx = (animal_thresholds['threshold'] - 0.02).abs().idxmin()
        default_r = animal_thresholds.loc[closest_idx, 'best_abs_r']
        default_lag = animal_thresholds.loc[closest_idx, 'best_lag']
    
    comparison_data.append({
        'animal': animal,
        'default_r': default_r,
        'default_lag': default_lag,
        'optimal_r': opt_row['best_abs_r'],
        'optimal_threshold': opt_row['optimal_threshold'],
        'optimal_lag': opt_row['optimal_lag'],
        'improvement': opt_row['best_abs_r'] - default_r,
        'pct_improvement': ((opt_row['best_abs_r'] - default_r) / default_r * 100) if default_r > 0 else 0
    })

comparison_df = pd.DataFrame(comparison_data).sort_values('improvement', ascending=False)

# Visualization
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Panel 1: Comparison of default vs optimal |r|
ax = axes[0]
x_pos = np.arange(len(comparison_df))
width = 0.35

bars1 = ax.bar(x_pos - width/2, comparison_df['default_r'], width, 
               label='Default (0.02)', color='lightcoral', alpha=0.8, edgecolor='black')
bars2 = ax.bar(x_pos + width/2, comparison_df['optimal_r'], width,
               label='Optimal', color='lightgreen', alpha=0.8, edgecolor='black')

ax.set_xlabel('Animal', fontsize=11)
ax.set_ylabel('|Pearson r|', fontsize=11)
ax.set_title('Default vs Optimal Threshold Performance', fontsize=12, fontweight='bold')
ax.set_xticks(x_pos)
ax.set_xticklabels(comparison_df['animal'], rotation=45, ha='right', fontsize=9)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3, axis='y')

# Panel 2: Improvement distribution
ax = axes[1]
improvements = comparison_df['improvement'].values
colors = ['green' if x > 0 else 'red' if x < 0 else 'gray' for x in improvements]
ax.bar(x_pos, improvements, color=colors, alpha=0.7, edgecolor='black')
ax.axhline(0, color='black', linestyle='-', linewidth=1)
ax.set_xlabel('Animal', fontsize=11)
ax.set_ylabel('Improvement in |r|', fontsize=11)
ax.set_title('Improvement from Optimizing Threshold', fontsize=12, fontweight='bold')
ax.set_xticks(x_pos)
ax.set_xticklabels(comparison_df['animal'], rotation=45, ha='right', fontsize=9)
ax.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.show()

# Print summary
print("Comparison of Default (0.02) vs Optimal Threshold:")
print("=" * 100)
print(comparison_df.to_string(index=False))
print("\n" + "=" * 100)
print(f"\nSummary:")
print(f"  Mean improvement: {comparison_df['improvement'].mean():.4f} (|r| units)")
print(f"  Median improvement: {comparison_df['improvement'].median():.4f}")
print(f"  Mean % improvement: {comparison_df['pct_improvement'].mean():.2f}%")
print(f"  Animals improved: {(comparison_df['improvement'] > 0).sum()} / {len(comparison_df)}")
print(f"  Animals worsened: {(comparison_df['improvement'] < 0).sum()} / {len(comparison_df)}")
print(f"  Max improvement: {comparison_df['improvement'].max():.4f} ({comparison_df.loc[comparison_df['improvement'].idxmax(), 'animal']})")

# %%

# %% [markdown]
# ## Time Moving Trajectories with Optimal Thresholds
#
# Visualize time_moving vs trial for each animal using their individualized optimal thresholds.

# %%
# Recalculate time_moving for each animal using their optimal threshold
optimal_time_moving_data = []

for idx, animal_opt in per_animal_df.iterrows():
    animal = animal_opt['animal']
    optimal_threshold = animal_opt['optimal_threshold']
    
    # Recalculate time_moving for entire x_array with optimal threshold
    x_array_with_optimal = recalculate_time_moving(
        x_array,
        snips_movement,
        threshold=optimal_threshold
    )
    
    # Get subset and animal's data
    subset_with_optimal = x_array_with_optimal.query("condition == 'deplete' & infusiontype == '45NaCl'")
    animal_data_optimal = subset_with_optimal.query("id == @animal").sort_values('trial').copy()
    
    # Get default data for comparison
    animal_data_default = subset_condition.query("id == @animal").sort_values('trial').copy()
    
    # Store for plotting
    optimal_time_moving_data.append({
        'animal': animal,
        'trials': animal_data_optimal['trial'].values,
        'time_moving_optimal': animal_data_optimal['time_moving'].values,
        'time_moving_default': animal_data_default['time_moving'].values,  # Original
        'optimal_threshold': optimal_threshold,
        'best_r': animal_opt['best_r']
    })

print(f"Recalculated time_moving for {len(optimal_time_moving_data)} animals using individualized thresholds")
print(f"Threshold range: [{per_animal_df['optimal_threshold'].min():.4f}, {per_animal_df['optimal_threshold'].max():.4f}]")

# %%
# Plot time_moving vs trial for each animal using optimal thresholds
n_animals = len(optimal_time_moving_data)
n_cols = 3
n_rows = int(np.ceil(n_animals / n_cols))

fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 4 * n_rows), sharex=False, sharey=True)
axes = axes.flatten() if n_animals > 1 else [axes]

for idx, animal_info in enumerate(optimal_time_moving_data):
    ax = axes[idx]
    
    animal = animal_info['animal']
    trials = animal_info['trials']
    time_moving_optimal = animal_info['time_moving_optimal']
    time_moving_default = animal_info['time_moving_default']
    optimal_threshold = animal_info['optimal_threshold']
    best_r = animal_info['best_r']
    
    # Plot both optimal and default for comparison
    ax.plot(trials, time_moving_default, 'o-', alpha=0.4, color='gray', 
            linewidth=1.5, markersize=4, label=f'Default (0.02)')
    ax.plot(trials, time_moving_optimal, 'o-', alpha=0.8, color=colors[1], 
            linewidth=2, markersize=5, label=f'Optimal ({optimal_threshold:.4f})')
    
    ax.set_xlabel('Trial', fontsize=10)
    ax.set_ylabel('Time Moving (proportion)', fontsize=10)
    ax.set_title(f'{animal} | r={best_r:.3f}', fontsize=11, fontweight='bold')
    ax.legend(fontsize=8, loc='best', framealpha=0.9)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-0.05, 1.05)

# Hide unused subplots
for idx in range(n_animals, len(axes)):
    axes[idx].axis('off')

plt.tight_layout()
plt.show()

print(f"\nShowing time_moving trajectories for all {n_animals} animals")
print("Gray line: default threshold (0.02)")
print("Colored line: individualized optimal threshold")

# %%
# Overlay plot: all animals' time_moving with optimal thresholds on one plot
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Panel 1: All animals with optimal thresholds
ax = axes[0]
for animal_info in sorted(optimal_time_moving_data, key=lambda x: x['best_r'], reverse=True):
    trials = animal_info['trials']
    time_moving_optimal = animal_info['time_moving_optimal']
    animal = animal_info['animal']
    best_r = animal_info['best_r']
    
    # Color by correlation strength
    alpha = 0.3 + 0.5 * (abs(best_r) / per_animal_df['best_abs_r'].max())
    ax.plot(trials, time_moving_optimal, 'o-', alpha=alpha, linewidth=1.5, 
            markersize=3, label=f'{animal} (r={best_r:.2f})')

ax.set_xlabel('Trial', fontsize=11)
ax.set_ylabel('Time Moving (proportion)', fontsize=11)
ax.set_title('All Animals: Time Moving with Optimal Thresholds', fontsize=12, fontweight='bold')
ax.grid(True, alpha=0.3)
ax.set_ylim(-0.05, 1.05)
ax.legend(fontsize=7, loc='center left', bbox_to_anchor=(1, 0.5), framealpha=0.9)

# Panel 2: Mean ± SD across animals
ax = axes[1]

# Align all animals to common trial grid (use min and max trials)
all_trials = []
for animal_info in optimal_time_moving_data:
    all_trials.extend(animal_info['trials'])
trial_min, trial_max = min(all_trials), max(all_trials)
common_trials = np.arange(trial_min, trial_max + 1)

# Interpolate each animal's data to common trial grid
time_moving_matrix = []
for animal_info in optimal_time_moving_data:
    trials = animal_info['trials']
    time_moving = animal_info['time_moving_optimal']
    
    # Interpolate to common grid
    time_moving_interp = np.interp(common_trials, trials, time_moving, 
                                   left=np.nan, right=np.nan)
    time_moving_matrix.append(time_moving_interp)

time_moving_matrix = np.array(time_moving_matrix)

# Calculate mean and std across animals (ignoring NaNs)
mean_time_moving = np.nanmean(time_moving_matrix, axis=0)
std_time_moving = np.nanstd(time_moving_matrix, axis=0)

ax.plot(common_trials, mean_time_moving, 'o-', color=colors[0], linewidth=2.5, 
        markersize=5, label='Mean across animals')
ax.fill_between(common_trials, 
                mean_time_moving - std_time_moving,
                mean_time_moving + std_time_moving,
                color=colors[0], alpha=0.3, label='±1 SD')

ax.set_xlabel('Trial', fontsize=11)
ax.set_ylabel('Time Moving (proportion)', fontsize=11)
ax.set_title('Mean Time Moving Across Animals\n(with optimal thresholds)', fontsize=12, fontweight='bold')
ax.legend(fontsize=10, loc='best', framealpha=0.9)
ax.grid(True, alpha=0.3)
ax.set_ylim(-0.05, 1.05)

plt.tight_layout()
plt.show()

print(f"Left panel: Individual animals (opacity scaled by |r|)")
print(f"Right panel: Population average ± SD")

# %%

# %% [markdown]
# ## Time Moving Trajectories (Threshold = 0.005)
#
# Visualize time_moving vs trial for each animal using a fixed threshold of 0.005.

# %%
# Recalculate time_moving using a fixed threshold for all animals
fixed_threshold = 3

x_array_fixed = recalculate_time_moving(
    x_array,
    snips_movement,
    threshold=fixed_threshold
)

subset_fixed = x_array_fixed.query("condition == 'deplete' & infusiontype == '45NaCl'")

best_r_by_animal = {}
if 'per_animal_df' in globals() and not per_animal_df.empty:
    best_r_by_animal = per_animal_df.set_index('animal')['best_r'].to_dict()

fixed_time_moving_data = []
for animal in animals_full:
    animal_data_fixed = subset_fixed.query("id == @animal").sort_values('trial').copy()
    if animal_data_fixed.empty:
        continue

    animal_data_default = subset_condition.query("id == @animal").sort_values('trial').copy()

    fixed_time_moving_data.append({
        'animal': animal,
        'trials': animal_data_fixed['trial'].values,
        'time_moving_fixed': animal_data_fixed['time_moving'].values,
        'time_moving_default': animal_data_default['time_moving'].values,
        'fixed_threshold': fixed_threshold,
        'best_r': best_r_by_animal.get(animal, np.nan)
    })

print(f"Recalculated time_moving for {len(fixed_time_moving_data)} animals using threshold={fixed_threshold}")

# %%
# Plot time_moving vs trial for each animal using fixed threshold
n_animals_fixed = len(fixed_time_moving_data)
n_cols = 3
n_rows = int(np.ceil(n_animals_fixed / n_cols))

fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 4 * n_rows), sharex=False, sharey=True)
axes = axes.flatten() if n_animals_fixed > 1 else [axes]

for idx, animal_info in enumerate(fixed_time_moving_data):
    ax = axes[idx]

    animal = animal_info['animal']
    trials = animal_info['trials']
    time_moving_fixed = animal_info['time_moving_fixed']
    time_moving_default = animal_info['time_moving_default']
    best_r = animal_info['best_r']

    # Plot both fixed threshold and default for comparison
    ax.plot(trials, time_moving_default, 'o-', alpha=0.4, color='gray',
            linewidth=1.5, markersize=4, label='Default (0.02)')
    ax.plot(trials, time_moving_fixed, 'o-', alpha=0.8, color=colors[1],
            linewidth=2, markersize=5, label=f'Fixed ({fixed_threshold:.3f})')

    title = f"{animal}"
    if np.isfinite(best_r):
        title += f" | r={best_r:.3f}"

    ax.set_xlabel('Trial', fontsize=10)
    ax.set_ylabel('Time Moving (proportion)', fontsize=10)
    ax.set_title(title, fontsize=11, fontweight='bold')
    ax.legend(fontsize=8, loc='best', framealpha=0.9)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-0.05, 1.05)

# Hide unused subplots
for idx in range(n_animals_fixed, len(axes)):
    axes[idx].axis('off')

plt.tight_layout()
plt.show()

print(f"\nShowing time_moving trajectories for all {n_animals_fixed} animals")
print("Gray line: default threshold (0.02)")
print(f"Colored line: fixed threshold ({fixed_threshold:.3f})")

# %%
# Overlay plot: all animals' time_moving with fixed threshold on one plot
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Panel 1: All animals with fixed threshold
ax = axes[0]
for animal_info in sorted(fixed_time_moving_data, key=lambda x: (np.nan_to_num(x['best_r'], nan=0.0)), reverse=True):
    trials = animal_info['trials']
    time_moving_fixed = animal_info['time_moving_fixed']
    animal = animal_info['animal']
    best_r = animal_info['best_r']

    alpha_scale = np.nan_to_num(abs(best_r), nan=0.0)
    alpha = 0.3 + 0.5 * (alpha_scale / max(per_animal_df['best_abs_r'].max(), 1e-6))
    ax.plot(trials, time_moving_fixed, 'o-', alpha=alpha, linewidth=1.5,
            markersize=3, label=f'{animal} (r={best_r:.2f})' if np.isfinite(best_r) else animal)

ax.set_xlabel('Trial', fontsize=11)
ax.set_ylabel('Time Moving (proportion)', fontsize=11)
ax.set_title(f'All Animals: Time Moving (threshold={fixed_threshold:.3f})', fontsize=12, fontweight='bold')
ax.grid(True, alpha=0.3)
ax.set_ylim(-0.05, 1.05)
ax.legend(fontsize=7, loc='center left', bbox_to_anchor=(1, 0.5), framealpha=0.9)

# Panel 2: Mean ± SD across animals
ax = axes[1]

# Align all animals to common trial grid (use min and max trials)
all_trials = []
for animal_info in fixed_time_moving_data:
    all_trials.extend(animal_info['trials'])
trial_min, trial_max = min(all_trials), max(all_trials)
common_trials = np.arange(trial_min, trial_max + 1)

# Interpolate each animal's data to common trial grid
time_moving_matrix = []
for animal_info in fixed_time_moving_data:
    trials = animal_info['trials']
    time_moving = animal_info['time_moving_fixed']

    # Interpolate to common grid
    time_moving_interp = np.interp(common_trials, trials, time_moving,
                                   left=np.nan, right=np.nan)
    time_moving_matrix.append(time_moving_interp)

time_moving_matrix = np.array(time_moving_matrix)

# Calculate mean and std across animals (ignoring NaNs)
mean_time_moving = np.nanmean(time_moving_matrix, axis=0)
std_time_moving = np.nanstd(time_moving_matrix, axis=0)

ax.plot(common_trials, mean_time_moving, 'o-', color=colors[0], linewidth=2.5,
        markersize=5, label='Mean across animals')
ax.fill_between(common_trials,
                mean_time_moving - std_time_moving,
                mean_time_moving + std_time_moving,
                color=colors[0], alpha=0.3, label='±1 SD')

ax.set_xlabel('Trial', fontsize=11)
ax.set_ylabel('Time Moving (proportion)', fontsize=11)
ax.set_title('Mean Time Moving Across Animals\n(fixed threshold)', fontsize=12, fontweight='bold')
ax.legend(fontsize=10, loc='best', framealpha=0.9)
ax.grid(True, alpha=0.3)
ax.set_ylim(-0.05, 1.05)

plt.tight_layout()
plt.show()

print("Left panel: Individual animals")
print("Right panel: Population average ± SD")

# %%

# %%

# %%

# %%
