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

from scipy.stats import linregress
import pingouin as pg

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

def make_auc_and_sem_data_per_trial(x_array, condition, infusiontype, y_col):
    
    grouped = (
        x_array
        .query("condition == @condition & infusiontype == @infusiontype")
        .groupby("trial")[y_col]
    )
    
    auc_data = (grouped.mean().values)
    sem_data = (grouped.sem().values)

    return auc_data, sem_data

# deplete_100_mean, deplete_100_sem = make_auc_and_sem_data_per_trial(x_array, "deplete", "10NaCl", "simba_median_balance")

# get xvals for diff infusion volumes
scaling_factor_100mM = 2.3376
scaling_factor_450mM = 2.3376 * 4.5

x_vals_100 = np.arange(0,49)*scaling_factor_100mM
x_vals_450 = np.arange(0,49)*scaling_factor_450mM


# %%
# make panels with trial
def make_trial_and_nacl_plot(x_vals_100, deplete_100_mean, deplete_100_sem,
                             x_vals_450, deplete_450_mean, deplete_450_sem,
                             colors=colors, scattersize=15):
    
    f, ax = plt.subplots(ncols=3, figsize=(6, 1.8), sharey=True,
                         gridspec_kw={"wspace": 0.5, "width_ratios": [0.5, 0.5, 1.2],
                                      "bottom": 0.3})
    
    x_vals_trial = np.arange(0,49)
    ax[0].scatter(x_vals_trial, deplete_100_mean,
                  color=colors[2], label="10 mM NaCl", alpha=0.5, s=scattersize)
    ax[0].fill_between(x_vals_trial,
                       deplete_100_mean - deplete_100_sem,
                       deplete_100_mean + deplete_100_sem,
                       color=colors[2], alpha=0.3, zorder=-10)
    
    ax[0].scatter(x_vals_trial, deplete_450_mean,
                  color=colors[3], label="45 mM NaCl", alpha=0.5, s=scattersize)
    ax[0].fill_between(x_vals_trial,
                       deplete_450_mean - deplete_450_sem,
                       deplete_450_mean + deplete_450_sem,
                       color=colors[3], alpha=0.3, zorder=-10)

    ax[0].set_xlabel("Trial")
    ax[0].axhline(0, color="k", linestyle="--",alpha=0.3, zorder=-10)
    
    ax[1].scatter(x_vals_100, deplete_100_mean,
                  color=colors[2], label="10 mM NaCl", alpha=0.3, s=scattersize)
    ax[1].fill_between(x_vals_100,
                       deplete_100_mean - deplete_100_sem,
                       deplete_100_mean + deplete_100_sem,
                       color=colors[2], alpha=0.3, zorder=-10)
    
    ax[1].scatter(x_vals_450, deplete_450_mean,
                  color=colors[3], label="45 mM NaCl", alpha=0.5, s=scattersize)
    ax[1].fill_between(x_vals_450,
                       deplete_450_mean - deplete_450_sem,
                       deplete_450_mean + deplete_450_sem,
                       color=colors[3], alpha=0.3, zorder=-10)

    ax[1].axvspan(np.min(x_vals_100), np.max(x_vals_100), color="k", alpha=0.05, zorder=-10)
    
    ax[1].set_xlabel("NaCl (mg)")
    
    ax[2].errorbar(x_vals_100, deplete_100_mean, yerr=deplete_100_sem, color=colors[2], alpha=0.5, marker="o")
    
    x_vals_450_red, deplete_450_red, deplete_450_sem_red = zip(*[(x, y, err) for x, y, err in zip(x_vals_450, deplete_450_mean, deplete_450_sem) if x <= np.max(x_vals_100)])
    
    
    ax[2].errorbar(x_vals_450_red, deplete_450_red, yerr=deplete_450_sem_red, color=colors[3], alpha=0.7, marker="o")
    
    ax[2].set_xlabel("NaCl (mg)")
    
    
    for axis in ax:
        sns.despine(ax=axis, offset=5)
        axis.axhline(0, color="k", linestyle="--",alpha=0.3, zorder=-10)
        
    return f, ax

deplete_100_mean, deplete_100_sem = make_auc_and_sem_data_per_trial(x_array, "deplete", "10NaCl", "simba_median_balance")
deplete_450_mean, deplete_450_sem = make_auc_and_sem_data_per_trial(x_array, "deplete", "45NaCl", "simba_median_balance")

f, ax = make_trial_and_nacl_plot(x_vals_100, deplete_100_mean, deplete_100_sem, x_vals_450, deplete_450_mean, deplete_450_sem, colors=colors)
ax[0].set_ylabel("Appetitive behaviour")
save_figure_atomic(f, "fig5_nacl_behav", FIGSFOLDER)

deplete_100_mean, deplete_100_sem = make_auc_and_sem_data_per_trial(x_array, "deplete", "10NaCl", "auc_snips")
deplete_450_mean, deplete_450_sem = make_auc_and_sem_data_per_trial(x_array, "deplete", "45NaCl", "auc_snips")

f, ax = make_trial_and_nacl_plot(x_vals_100, deplete_100_mean, deplete_100_sem, x_vals_450, deplete_450_mean, deplete_450_sem, colors=colors)
ax[0].set_ylabel("Dopamine (AUC)")
save_figure_atomic(f, "fig5_nacl_da", FIGSFOLDER)

# %%
# first testing behavioural data with anova (function of trial)

beh_df = (
    x_array
    .query('condition == "deplete" and infusiontype in ["10NaCl", "45NaCl"]')
    [["id", "infusiontype", "trial", "simba_median_balance"]]
    .dropna()
    .groupby([ "id", "infusiontype", "trial"], as_index=False)["simba_median_balance"]
    .mean()
)

beh_anova = pg.mixed_anova(
    data=beh_df,
    dv="simba_median_balance",
    within="trial",
    between="infusiontype",
    subject="id",
)

display(beh_anova)

# %%
# behaviour ANOVA using sodium delivered (mg) as within-subject factor
# binned into 5 levels based on common bins across infusion types

beh_mg_df = (
    x_array
    .query('condition == "deplete" and infusiontype in ["10NaCl", "45NaCl"]')
    [["id", "infusiontype", "trial", "simba_median_balance"]]
    .dropna()
    .groupby(["id", "infusiontype", "trial"], as_index=False)["simba_median_balance"]
    .mean()
)

scale_map = {
    "10NaCl": scaling_factor_100mM,
    "45NaCl": scaling_factor_450mM,
}
beh_mg_df["sodium_mg"] = (beh_mg_df["trial"] + 1) * beh_mg_df["infusiontype"].map(scale_map)

# Use common sodium bins so the within-subject levels are shared across subjects.
mg_max_overlap = beh_mg_df[beh_mg_df["infusiontype"] == "10NaCl"]["sodium_mg"].max()
beh_mg_df = beh_mg_df[beh_mg_df["sodium_mg"] <= mg_max_overlap].copy()
n_bins = 5
beh_mg_df["sodium_bin"] = pd.cut(
    beh_mg_df["sodium_mg"],
    bins=np.linspace(0, mg_max_overlap, n_bins + 1),
    labels=[f"bin_{i}" for i in range(1, n_bins + 1)],
    include_lowest=True,
)

beh_mg_df = (
    beh_mg_df
    .dropna(subset=["sodium_bin"])
    .groupby(["id", "infusiontype", "sodium_bin"], as_index=False)["simba_median_balance"]
    .mean()
)

beh_mg_anova = pg.mixed_anova(
    data=beh_mg_df,
    dv="simba_median_balance",
    within="sodium_bin",
    between="infusiontype",
    subject="id",
)

display(beh_mg_anova)

# %%
# behaviour continuous mixed model using sodium delivered (mg)

import statsmodels.formula.api as smf

beh_cont_df = (
    x_array
    .query('condition == "deplete" and infusiontype in ["10NaCl", "45NaCl"]')
    [["id", "infusiontype", "trial", "simba_median_balance"]]
    .dropna()
    .groupby(["id", "infusiontype", "trial"], as_index=False)["simba_median_balance"]
    .mean()
)

beh_cont_scale_map = {
    "10NaCl": scaling_factor_100mM,
    "45NaCl": scaling_factor_450mM,
}
beh_cont_df["sodium_mg"] = (beh_cont_df["trial"] + 1) * beh_cont_df["infusiontype"].map(beh_cont_scale_map)

beh_cont_model = smf.mixedlm(
    "simba_median_balance ~ sodium_mg * C(infusiontype)",
    data=beh_cont_df,
    groups=beh_cont_df["id"],
).fit(reml=False)

display(pd.DataFrame({"coef": beh_cont_model.params, "p_value": beh_cont_model.pvalues}))

# %%
# dopamine ANOVA as a function of trial

da_df = (
    x_array
    .query('condition == "deplete" and infusiontype in ["10NaCl", "45NaCl"]')
    [["id", "infusiontype", "trial", "auc_snips"]]
    .dropna()
    .groupby(["id", "infusiontype", "trial"], as_index=False)["auc_snips"]
    .mean()
)

da_anova = pg.mixed_anova(
    data=da_df,
    dv="auc_snips",
    within="trial",
    between="infusiontype",
    subject="id",
)

display(da_anova)

# %%
# dopamine ANOVA using sodium delivered (mg) as within-subject factor
# binned into 5 levels based on common bins across infusion types

da_mg_df = (
    x_array
    .query('condition == "deplete" and infusiontype in ["10NaCl", "45NaCl"]')
    [["id", "infusiontype", "trial", "auc_snips"]]
    .dropna()
    .groupby(["id", "infusiontype", "trial"], as_index=False)["auc_snips"]
    .mean()
)

da_scale_map = {
    "10NaCl": scaling_factor_100mM,
    "45NaCl": scaling_factor_450mM,
}
da_mg_df["sodium_mg"] = (da_mg_df["trial"] + 1) * da_mg_df["infusiontype"].map(da_scale_map)

# Use common sodium bins so the within-subject levels are shared across subjects.
da_mg_max_overlap = da_mg_df[da_mg_df["infusiontype"] == "10NaCl"]["sodium_mg"].max()
da_mg_df = da_mg_df[da_mg_df["sodium_mg"] <= da_mg_max_overlap].copy()
n_bins = 5
da_mg_df["sodium_bin"] = pd.cut(
    da_mg_df["sodium_mg"],
    bins=np.linspace(0, da_mg_max_overlap, n_bins + 1),
    labels=[f"bin_{i}" for i in range(1, n_bins + 1)],
    include_lowest=True,
)

da_mg_df = (
    da_mg_df
    .dropna(subset=["sodium_bin"])
    .groupby(["id", "infusiontype", "sodium_bin"], as_index=False)["auc_snips"]
    .mean()
)

da_mg_anova = pg.mixed_anova(
    data=da_mg_df,
    dv="auc_snips",
    within="sodium_bin",
    between="infusiontype",
    subject="id",
)

display(da_mg_anova)

# %%
# dopamine continuous mixed model using sodium delivered (mg)

import statsmodels.formula.api as smf

da_cont_df = (
    x_array
    .query('condition == "deplete" and infusiontype in ["10NaCl", "45NaCl"]')
    [["id", "infusiontype", "trial", "auc_snips"]]
    .dropna()
    .groupby(["id", "infusiontype", "trial"], as_index=False)["auc_snips"]
    .mean()
)

da_cont_scale_map = {
    "10NaCl": scaling_factor_100mM,
    "45NaCl": scaling_factor_450mM,
}
da_cont_df["sodium_mg"] = (da_cont_df["trial"] + 1) * da_cont_df["infusiontype"].map(da_cont_scale_map)

# Stabilize optimization by scaling the continuous predictor.
sodium_std = da_cont_df["sodium_mg"].std(ddof=0)
if sodium_std == 0:
    raise ValueError("sodium_mg has zero variance; cannot fit continuous model.")
da_cont_df["sodium_mg_z"] = (da_cont_df["sodium_mg"] - da_cont_df["sodium_mg"].mean()) / sodium_std

da_cont_formula = "auc_snips ~ sodium_mg_z * C(infusiontype)"
da_cont_mixedlm = smf.mixedlm(
    da_cont_formula,
    data=da_cont_df,
    groups=da_cont_df["id"],
)

da_cont_model = da_cont_mixedlm.fit(method="lbfgs", maxiter=2000, reml=True)
if not da_cont_model.converged:
    da_cont_model = da_cont_mixedlm.fit(method="powell", maxiter=2000, reml=True)

display(pd.DataFrame({"coef": da_cont_model.params, "p_value": da_cont_model.pvalues}))
print(f"converged: {da_cont_model.converged}")

# %%
# dopamine continuous fallback: OLS with subject-clustered SEs

da_cont_ols = smf.ols(da_cont_formula, data=da_cont_df).fit(
    cov_type="cluster",
    cov_kwds={"groups": da_cont_df["id"]},
)

display(pd.DataFrame({"coef": da_cont_ols.params, "p_value": da_cont_ols.pvalues}))

# %%
# dopamine continuous fallback with early-trial trimming by infusion group
# remove first 13 trials in 10NaCl and first 3 trials in 45NaCl

da_trim_df = da_cont_df.copy()
da_trim_df = da_trim_df[
    ((da_trim_df["infusiontype"] == "10NaCl") & (da_trim_df["trial"] >= 13))
    | ((da_trim_df["infusiontype"] == "45NaCl") & (da_trim_df["trial"] >= 3))
].copy()

sodium_std_trim = da_trim_df["sodium_mg"].std(ddof=0)
if sodium_std_trim == 0:
    raise ValueError("sodium_mg has zero variance after trimming; cannot fit model.")
da_trim_df["sodium_mg_z"] = (
    da_trim_df["sodium_mg"] - da_trim_df["sodium_mg"].mean()
) / sodium_std_trim

da_trim_formula = "auc_snips ~ sodium_mg_z * C(infusiontype)"
da_trim_ols = smf.ols(da_trim_formula, data=da_trim_df).fit(
    cov_type="cluster",
    cov_kwds={"groups": da_trim_df["id"]},
)

display(pd.DataFrame({"coef": da_trim_ols.params, "p_value": da_trim_ols.pvalues}))
print(f"n rows after trimming: {len(da_trim_df)}")
print(f"n subjects after trimming: {da_trim_df['id'].nunique()}")

# %%
# behaviour continuous fallback: OLS with subject-clustered SEs

beh_ols_df = beh_cont_df.copy()
beh_sodium_std = beh_ols_df["sodium_mg"].std(ddof=0)
if beh_sodium_std == 0:
    raise ValueError("sodium_mg has zero variance; cannot fit behaviour OLS model.")
beh_ols_df["sodium_mg_z"] = (
    beh_ols_df["sodium_mg"] - beh_ols_df["sodium_mg"].mean()
) / beh_sodium_std

beh_ols_formula = "simba_median_balance ~ sodium_mg_z * C(infusiontype)"
beh_cont_ols = smf.ols(beh_ols_formula, data=beh_ols_df).fit(
    cov_type="cluster",
    cov_kwds={"groups": beh_ols_df["id"]},
)

display(pd.DataFrame({"coef": beh_cont_ols.params, "p_value": beh_cont_ols.pvalues}))

# %%
# behaviour continuous fallback with early-trial trimming by infusion group
# remove first 13 trials in 10NaCl and first 3 trials in 45NaCl

beh_trim_df = beh_ols_df.copy()
beh_trim_df = beh_trim_df[
    ((beh_trim_df["infusiontype"] == "10NaCl") & (beh_trim_df["trial"] >= 13))
    | ((beh_trim_df["infusiontype"] == "45NaCl") & (beh_trim_df["trial"] >= 3))
].copy()

beh_sodium_std_trim = beh_trim_df["sodium_mg"].std(ddof=0)
if beh_sodium_std_trim == 0:
    raise ValueError("sodium_mg has zero variance after trimming; cannot fit behaviour model.")
beh_trim_df["sodium_mg_z"] = (
    beh_trim_df["sodium_mg"] - beh_trim_df["sodium_mg"].mean()
) / beh_sodium_std_trim

beh_trim_formula = "simba_median_balance ~ sodium_mg_z * C(infusiontype)"
beh_trim_ols = smf.ols(beh_trim_formula, data=beh_trim_df).fit(
    cov_type="cluster",
    cov_kwds={"groups": beh_trim_df["id"]},
)

display(pd.DataFrame({"coef": beh_trim_ols.params, "p_value": beh_trim_ols.pvalues}))
print(f"n rows after trimming: {len(beh_trim_df)}")
print(f"n subjects after trimming: {beh_trim_df['id'].nunique()}")

# %%
da_trim_df


# %%
def make_auc_and_sem_data_per_trial_for_trimmed_arrays(trimmed_df, infusiontype, y_col):
    
    grouped = (
        trimmed_df
        .query("infusiontype == @infusiontype")
        .groupby("trial")[y_col]
    )
    
    auc_data = (grouped.mean().values)
    sem_data = (grouped.sem().values)

    return auc_data, sem_data


deplete_100_mean, deplete_100_sem = make_auc_and_sem_data_per_trial_for_trimmed_arrays(da_trim_df, "10NaCl", "auc_snips")
deplete_450_mean, deplete_450_sem = make_auc_and_sem_data_per_trial_for_trimmed_arrays(da_trim_df, "45NaCl", "auc_snips")


# %%

# get xvals for diff infusion volumes
scaling_factor_100mM = 2.3376
scaling_factor_450mM = 2.3376 * 4.5

x_vals_100 = np.arange(0,36)*scaling_factor_100mM
x_vals_450 = np.arange(0,36)*scaling_factor_450mM

# %%
# supplemental figure showing dopamine with first x mg of NaCl removed

f, ax = plt.subplots(ncols=1, figsize=(3, 2))



ax.errorbar(x_vals_100, deplete_100_mean, yerr=deplete_100_sem, color=colors[2], alpha=0.5, marker="o")

x_vals_450_red, deplete_450_red, deplete_450_sem_red = zip(*[(x, y, err) for x, y, err in zip(x_vals_450, deplete_450_mean, deplete_450_sem) if x <= np.max(x_vals_100)])


ax.errorbar(x_vals_450_red, deplete_450_red, yerr=deplete_450_sem_red, color=colors[3], alpha=0.7, marker="o")

ax.set_xlabel("NaCl (mg)")
ax.set_ylabel("Dopamine AUC")

sns.despine(ax=ax, offset=5)
ax.axhline(0, color="k", linestyle="--",alpha=0.3, zorder=-10)

# %%
