# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.2
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
# next testing dopamine data with anova (function of trial)

beh_df = (
    x_array
    .query('condition == "deplete" and infusiontype in ["10NaCl", "45NaCl"]')
    [["id", "infusiontype", "trial", "auc_snips"]]
    .dropna()
    .groupby([ "id", "infusiontype", "trial"], as_index=False)["auc_snips"]
    .mean()
)

beh_anova = pg.mixed_anova(
    data=beh_df,
    dv="auc_snips",
    within="trial",
    between="infusiontype",
    subject="id",
)

display(beh_anova)

# %%
x_array.columns

# %%
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf


def _pick_subject_col(df):
    candidates = [
        "id", "rat", "rat_id", "ratid", "subject", "subject_id", "animal", "animal_id", "rat_n"
    ]
    for c in candidates:
        if c in df.columns:
            return c
    raise ValueError(
        "Could not find a subject identifier column. "
        "Expected one of: id, rat, rat_id, ratid, subject, subject_id, animal, animal_id, rat_n"
    )


def compare_trial_vs_mg_interactions(
    x_array,
    y_col="simba_median_balance",
    condition="deplete",
    infusion_levels=("10NaCl", "45NaCl"),
    scale_10=2.3376,
    scale_45=2.3376 * 4.5,
):
    cols_needed = ["condition", "infusiontype", "trial", y_col]
    missing = [c for c in cols_needed if c not in x_array.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    df = (
        x_array
        .query("condition == @condition and infusiontype in @infusion_levels")
        .copy()
    )

    subj_col = _pick_subject_col(df)

    # Keep one row per subject/condition/trial for this outcome.
    df = (
        df[[subj_col, "infusiontype", "trial", y_col]]
        .dropna(subset=[y_col, "trial", "infusiontype", subj_col])
        .groupby([subj_col, "infusiontype", "trial"], as_index=False)[y_col]
        .mean()
    )

    if df.empty:
        raise ValueError("No rows left after filtering/aggregation.")

    scale_map = {"10NaCl": scale_10, "45NaCl": scale_45}
    df["nacl_mg"] = (df["trial"] + 1) * df["infusiontype"].map(scale_map)

    # Center for numeric stability.
    df["trial_c"] = df["trial"] - df["trial"].mean()
    df["nacl_mg_c"] = df["nacl_mg"] - df["nacl_mg"].mean()

    formula_trial = f"{y_col} ~ trial_c * C(infusiontype)"
    formula_mg = f"{y_col} ~ nacl_mg_c * C(infusiontype)"

    # Subject-clustered SEs account for repeated measures within subject.
    model_trial = smf.ols(formula_trial, data=df).fit(
        cov_type="cluster", cov_kwds={"groups": df[subj_col]}
    )
    model_mg = smf.ols(formula_mg, data=df).fit(
        cov_type="cluster", cov_kwds={"groups": df[subj_col]}
    )

    inter_trial_name = [n for n in model_trial.pvalues.index if "trial_c:C(infusiontype)" in n][0]
    inter_mg_name = [n for n in model_mg.pvalues.index if "nacl_mg_c:C(infusiontype)" in n][0]

    out = pd.DataFrame([
        {
            "analysis": "Outcome ~ Trial * Infusion",
            "term": inter_trial_name,
            "coef": model_trial.params[inter_trial_name],
            "p_value": model_trial.pvalues[inter_trial_name],
            "n_rows": len(df),
            "n_subjects": df[subj_col].nunique(),
        },
        {
            "analysis": "Outcome ~ NaCl(mg) * Infusion",
            "term": inter_mg_name,
            "coef": model_mg.params[inter_mg_name],
            "p_value": model_mg.pvalues[inter_mg_name],
            "n_rows": len(df),
            "n_subjects": df[subj_col].nunique(),
        },
    ])

    return out, model_trial, model_mg, df


# Primary comparison requested: behaviour
behav_results, behav_trial_model, behav_mg_model, behav_df = compare_trial_vs_mg_interactions(
    x_array, y_col="simba_median_balance"
)

print("Behaviour interaction tests")
display(behav_results)

p_trial = behav_results.loc[behav_results["analysis"] == "Outcome ~ Trial * Infusion", "p_value"].iloc[0]
p_mg = behav_results.loc[behav_results["analysis"] == "Outcome ~ NaCl(mg) * Infusion", "p_value"].iloc[0]

print("\nInterpretation:")
print(f"- Trial x infusion interaction p = {p_trial:.3g}")
print(f"- NaCl(mg) x infusion interaction p = {p_mg:.3g}")
if (p_trial < 0.05) and (p_mg >= 0.05):
    print("- Pattern matches your hypothesis: significant by trial, not significant by NaCl consumed.")
else:
    print("- Pattern does not exactly match the expected trial-vs-NaCl contrast in this run.")

# %%
# Dopamine-only interaction comparison (trial vs NaCl consumed)
da_results, da_trial_model, da_mg_model, da_df = compare_trial_vs_mg_interactions(
    x_array, y_col="auc_snips"
)

print("Dopamine interaction tests")
display(da_results)

p_trial_da = da_results.loc[
    da_results["analysis"] == "Outcome ~ Trial * Infusion", "p_value"
].iloc[0]
p_mg_da = da_results.loc[
    da_results["analysis"] == "Outcome ~ NaCl(mg) * Infusion", "p_value"
].iloc[0]

print("\nInterpretation:")
print(f"- Trial x infusion interaction p = {p_trial_da:.3g}")
print(f"- NaCl(mg) x infusion interaction p = {p_mg_da:.3g}")
if (p_trial_da < 0.05) and (p_mg_da >= 0.05):
    print("- Pattern matches your hypothesis: significant by trial, not significant by NaCl consumed.")
else:
    print("- Pattern does not exactly match the expected trial-vs-NaCl contrast in this run.")

# %%
import statsmodels.formula.api as smf
from scipy.stats import chi2


def _lr_test(full_model, reduced_model):
    lr = 2 * (full_model.llf - reduced_model.llf)
    df_diff = int(full_model.df_modelwc - reduced_model.df_modelwc)
    p = chi2.sf(lr, df_diff)
    return lr, df_diff, p


def _format_wald_terms(model_result):
    wt = model_result.wald_test_terms(skip_single=False)
    table = wt.table.reset_index().rename(columns={"index": "term"})

    # Normalize possible column names across statsmodels versions.
    col_map = {}
    for c in table.columns:
        cl = str(c).lower()
        if "chi2" in cl or "stat" in cl:
            col_map[c] = "chi2"
        elif "p>" in cl or "pvalue" in cl or cl == "p":
            col_map[c] = "p_value"
        elif "df" in cl:
            col_map[c] = "df"
    table = table.rename(columns=col_map)

    keep_cols = [c for c in ["term", "chi2", "df", "p_value"] if c in table.columns]
    table = table[keep_cols].copy()
    return table


# --------------------------------------------
# Dopamine mixed-effects models (easy to read)
# Infusion = between-subject factor
# --------------------------------------------

subject_col = _pick_subject_col(x_array)

da_long = (
    x_array
    .query('condition == "deplete" and infusiontype in ["10NaCl", "45NaCl"]')
    [[subject_col, "infusiontype", "trial", "auc_snips"]]
    .dropna()
    .groupby([subject_col, "infusiontype", "trial"], as_index=False)["auc_snips"]
    .mean()
)

# 1) Trial-bin model (within-subject bins)
da_trial = da_long.copy()
da_trial["trial_bin"] = pd.cut(
    da_trial["trial"],
    bins=[-0.1, 15.5, 31.5, 48.5],
    labels=["early", "mid", "late"],
)

da_trial = (
    da_trial
    .dropna(subset=["trial_bin"])
    .groupby([subject_col, "infusiontype", "trial_bin"], as_index=False)["auc_snips"]
    .mean()
)

trial_full = smf.mixedlm(
    "auc_snips ~ C(infusiontype) * C(trial_bin)",
    data=da_trial,
    groups=da_trial[subject_col],
).fit(reml=False)

trial_reduced = smf.mixedlm(
    "auc_snips ~ C(infusiontype) + C(trial_bin)",
    data=da_trial,
    groups=da_trial[subject_col],
).fit(reml=False)

trial_lr, trial_df, trial_p_interaction = _lr_test(trial_full, trial_reduced)
trial_terms = _format_wald_terms(trial_full)

# 2) NaCl(mg)-bin model (within-subject bins in shared mg range)
scale_10 = 2.3376
scale_45 = 2.3376 * 4.5
scale_map = {"10NaCl": scale_10, "45NaCl": scale_45}

da_mg = da_long.copy()
da_mg["nacl_mg"] = (da_mg["trial"] + 1) * da_mg["infusiontype"].map(scale_map)

# Restrict NaCl-bin analyses to the 0.1M coverage requested by user.
mg_max_overlap = 120.0
da_mg = da_mg[da_mg["nacl_mg"] <= mg_max_overlap].copy()

da_mg["mg_bin"] = pd.cut(
    da_mg["nacl_mg"],
    bins=np.linspace(0, mg_max_overlap, 4),
    labels=["low", "mid", "high"],
    include_lowest=True,
)

da_mg = (
    da_mg
    .dropna(subset=["mg_bin"])
    .groupby([subject_col, "infusiontype", "mg_bin"], as_index=False)["auc_snips"]
    .mean()
)

mg_full = smf.mixedlm(
    "auc_snips ~ C(infusiontype) * C(mg_bin)",
    data=da_mg,
    groups=da_mg[subject_col],
).fit(reml=False)

mg_reduced = smf.mixedlm(
    "auc_snips ~ C(infusiontype) + C(mg_bin)",
    data=da_mg,
    groups=da_mg[subject_col],
).fit(reml=False)

mg_lr, mg_df, mg_p_interaction = _lr_test(mg_full, mg_reduced)
mg_terms = _format_wald_terms(mg_full)

interaction_summary = pd.DataFrame([
    {
        "model": "Dopamine ~ infusion * trial_bin",
        "interaction_test": "LR(full vs no interaction)",
        "chi2": trial_lr,
        "df": trial_df,
        "p_value": trial_p_interaction,
        "n_subjects": da_trial[subject_col].nunique(),
    },
    {
        "model": "Dopamine ~ infusion * mg_bin",
        "interaction_test": "LR(full vs no interaction)",
        "chi2": mg_lr,
        "df": mg_df,
        "p_value": mg_p_interaction,
        "n_subjects": da_mg[subject_col].nunique(),
    },
])

print("Dopamine mixed-model interaction tests (infusion as between-subject):")
display(interaction_summary)

print("\nMain effects + interaction (trial-bin model):")
display(trial_terms)

print("Main effects + interaction (mg-bin model):")
display(mg_terms)

print("\nHow to read this:")
print("- In each term table, look for C(infusiontype), C(trial_bin)/C(mg_bin), and their interaction.")
print("- C(infusiontype): between-subject main effect.")
print("- C(trial_bin) or C(mg_bin): within-subject main effect.")
print("- Interaction term: whether infusion differences depend on trial-bin or mg-bin.")

# %%
import numpy as np


def _to_scalar(x):
    arr = np.asarray(x)
    if arr.size == 0:
        return np.nan
    return float(arr.ravel()[0])


def _get_term_stats(term_table, term_name):
    row = term_table.loc[term_table["term"] == term_name].iloc[0]
    return {
        "chi2": _to_scalar(row["chi2"]),
        "df": int(_to_scalar(row["df"])),
        "p": _to_scalar(row["p_value"]),
    }


def _fmt_p(p):
    if p < 0.001:
        return "p < 0.001"
    return f"p = {p:.3f}"


def _effect_phrase(p, positive_label, null_label, trend_label=None):
    if p < 0.05:
        return positive_label
    if (trend_label is not None) and (p < 0.10):
        return trend_label
    return null_label


# Pull key stats from the mixed-model term tables.
trial_inf = _get_term_stats(trial_terms, "C(infusiontype)")
trial_bin = _get_term_stats(trial_terms, "C(trial_bin)")
trial_int = _get_term_stats(trial_terms, "C(infusiontype):C(trial_bin)")

mg_inf = _get_term_stats(mg_terms, "C(infusiontype)")
mg_bin = _get_term_stats(mg_terms, "C(mg_bin)")
mg_int = _get_term_stats(mg_terms, "C(infusiontype):C(mg_bin)")

trial_sentence = (
    "Trial-bin model: dopamine showed "
    f"{_effect_phrase(trial_inf['p'], 'a main effect of infusion', 'no main effect of infusion')} "
    f"(chi2({trial_inf['df']}) = {trial_inf['chi2']:.2f}, {_fmt_p(trial_inf['p'])}), "
    f"{_effect_phrase(trial_bin['p'], 'a main effect of trial bin', 'no main effect of trial bin')} "
    f"(chi2({trial_bin['df']}) = {trial_bin['chi2']:.2f}, {_fmt_p(trial_bin['p'])}), "
    f"and {_effect_phrase(trial_int['p'], 'an infusion x trial-bin interaction', 'no infusion x trial-bin interaction')} "
    f"(chi2({trial_int['df']}) = {trial_int['chi2']:.2f}, {_fmt_p(trial_int['p'])})."
)

mg_sentence = (
    "NaCl-bin model: dopamine showed "
    f"{_effect_phrase(mg_inf['p'], 'a main effect of infusion', 'no main effect of infusion', trend_label='a trend toward a main effect of infusion')} "
    f"(chi2({mg_inf['df']}) = {mg_inf['chi2']:.2f}, {_fmt_p(mg_inf['p'])}), "
    f"{_effect_phrase(mg_bin['p'], 'a main effect of NaCl bin', 'no main effect of NaCl bin')} "
    f"(chi2({mg_bin['df']}) = {mg_bin['chi2']:.2f}, {_fmt_p(mg_bin['p'])}), "
    f"and {_effect_phrase(mg_int['p'], 'an infusion x NaCl-bin interaction', 'no infusion x NaCl-bin interaction', trend_label='a trend toward an infusion x NaCl-bin interaction')} "
    f"(chi2({mg_int['df']}) = {mg_int['chi2']:.2f}, {_fmt_p(mg_int['p'])})."
)

print("Manuscript-ready summary (dopamine):\n")
print(trial_sentence)
print(mg_sentence)

# %%
# Behaviour mixed-model tests and manuscript-ready summary

beh_long = (
    x_array
    .query('condition == "deplete" and infusiontype in ["10NaCl", "45NaCl"]')
    [[subject_col, "infusiontype", "trial", "simba_median_balance"]]
    .dropna()
    .groupby([subject_col, "infusiontype", "trial"], as_index=False)["simba_median_balance"]
    .mean()
)

# 1) Trial-bin model
beh_trial = beh_long.copy()
beh_trial["trial_bin"] = pd.cut(
    beh_trial["trial"],
    bins=[-0.1, 15.5, 31.5, 48.5],
    labels=["early", "mid", "late"],
)

beh_trial = (
    beh_trial
    .dropna(subset=["trial_bin"])
    .groupby([subject_col, "infusiontype", "trial_bin"], as_index=False)["simba_median_balance"]
    .mean()
)

beh_trial_full = smf.mixedlm(
    "simba_median_balance ~ C(infusiontype) * C(trial_bin)",
    data=beh_trial,
    groups=beh_trial[subject_col],
).fit(reml=False)

beh_trial_reduced = smf.mixedlm(
    "simba_median_balance ~ C(infusiontype) + C(trial_bin)",
    data=beh_trial,
    groups=beh_trial[subject_col],
).fit(reml=False)

beh_trial_lr, beh_trial_df, beh_trial_p_interaction = _lr_test(beh_trial_full, beh_trial_reduced)
beh_trial_terms = _format_wald_terms(beh_trial_full)

# 2) NaCl(mg)-bin model
beh_mg = beh_long.copy()
beh_mg["nacl_mg"] = (beh_mg["trial"] + 1) * beh_mg["infusiontype"].map(scale_map)
beh_mg = beh_mg[beh_mg["nacl_mg"] <= mg_max_overlap].copy()

beh_mg["mg_bin"] = pd.cut(
    beh_mg["nacl_mg"],
    bins=np.linspace(0, mg_max_overlap, 4),
    labels=["low", "mid", "high"],
    include_lowest=True,
)

beh_mg = (
    beh_mg
    .dropna(subset=["mg_bin"])
    .groupby([subject_col, "infusiontype", "mg_bin"], as_index=False)["simba_median_balance"]
    .mean()
)

beh_mg_full = smf.mixedlm(
    "simba_median_balance ~ C(infusiontype) * C(mg_bin)",
    data=beh_mg,
    groups=beh_mg[subject_col],
).fit(reml=False)

beh_mg_reduced = smf.mixedlm(
    "simba_median_balance ~ C(infusiontype) + C(mg_bin)",
    data=beh_mg,
    groups=beh_mg[subject_col],
).fit(reml=False)

beh_mg_lr, beh_mg_df, beh_mg_p_interaction = _lr_test(beh_mg_full, beh_mg_reduced)
beh_mg_terms = _format_wald_terms(beh_mg_full)

beh_interaction_summary = pd.DataFrame([
    {
        "model": "Behaviour ~ infusion * trial_bin",
        "interaction_test": "LR(full vs no interaction)",
        "chi2": beh_trial_lr,
        "df": beh_trial_df,
        "p_value": beh_trial_p_interaction,
        "n_subjects": beh_trial[subject_col].nunique(),
    },
    {
        "model": "Behaviour ~ infusion * mg_bin",
        "interaction_test": "LR(full vs no interaction)",
        "chi2": beh_mg_lr,
        "df": beh_mg_df,
        "p_value": beh_mg_p_interaction,
        "n_subjects": beh_mg[subject_col].nunique(),
    },
])

print("Behaviour mixed-model interaction tests (infusion as between-subject):")
display(beh_interaction_summary)

print("\nMain effects + interaction (trial-bin model):")
display(beh_trial_terms)

print("Main effects + interaction (mg-bin model):")
display(beh_mg_terms)

# Manuscript-ready text
beh_trial_inf = _get_term_stats(beh_trial_terms, "C(infusiontype)")
beh_trial_bin = _get_term_stats(beh_trial_terms, "C(trial_bin)")
beh_trial_int = _get_term_stats(beh_trial_terms, "C(infusiontype):C(trial_bin)")

beh_mg_inf = _get_term_stats(beh_mg_terms, "C(infusiontype)")
beh_mg_bin = _get_term_stats(beh_mg_terms, "C(mg_bin)")
beh_mg_int = _get_term_stats(beh_mg_terms, "C(infusiontype):C(mg_bin)")

beh_trial_sentence = (
    "Trial-bin model: behaviour showed "
    f"{_effect_phrase(beh_trial_inf['p'], 'a main effect of infusion', 'no main effect of infusion')} "
    f"(chi2({beh_trial_inf['df']}) = {beh_trial_inf['chi2']:.2f}, {_fmt_p(beh_trial_inf['p'])}), "
    f"{_effect_phrase(beh_trial_bin['p'], 'a main effect of trial bin', 'no main effect of trial bin')} "
    f"(chi2({beh_trial_bin['df']}) = {beh_trial_bin['chi2']:.2f}, {_fmt_p(beh_trial_bin['p'])}), "
    f"and {_effect_phrase(beh_trial_int['p'], 'an infusion x trial-bin interaction', 'no infusion x trial-bin interaction')} "
    f"(chi2({beh_trial_int['df']}) = {beh_trial_int['chi2']:.2f}, {_fmt_p(beh_trial_int['p'])})."
)

beh_mg_sentence = (
    "NaCl-bin model: behaviour showed "
    f"{_effect_phrase(beh_mg_inf['p'], 'a main effect of infusion', 'no main effect of infusion')} "
    f"(chi2({beh_mg_inf['df']}) = {beh_mg_inf['chi2']:.2f}, {_fmt_p(beh_mg_inf['p'])}), "
    f"{_effect_phrase(beh_mg_bin['p'], 'a main effect of NaCl bin', 'no main effect of NaCl bin')} "
    f"(chi2({beh_mg_bin['df']}) = {beh_mg_bin['chi2']:.2f}, {_fmt_p(beh_mg_bin['p'])}), "
    f"and {_effect_phrase(beh_mg_int['p'], 'an infusion x NaCl-bin interaction', 'no infusion x NaCl-bin interaction')} "
    f"(chi2({beh_mg_int['df']}) = {beh_mg_int['chi2']:.2f}, {_fmt_p(beh_mg_int['p'])})."
)

print("\nManuscript-ready summary (behaviour):\n")
print(beh_trial_sentence)
print(beh_mg_sentence)

# %%
