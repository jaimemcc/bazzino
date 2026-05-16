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
# # Notebook for finding sigmoidal transitions in behaviour and dopamine

# %%
import sys
from pathlib import Path

import dill
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from IPython.display import display

from trompy import save_figure_atomic

ROOT = Path.cwd()
if ROOT.name == "notebooks":
    ROOT = ROOT.parent

SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from model_fit_helpers import (
    fit_continuous_sigmoid_transitions_by_id,
    fit_logistic_transitions_by_id,
    build_realigned_trials,
    fit_curve_series,
    binarize_series_threshold,
    quantize_series_threshold,
    fit_logistic_per_series,
    fit_continuous_sigmoid_per_series,
)

from figure_config import COLORS, FIGSFOLDER, DATAFOLDER

# Use absolute paths so notebook behavior is robust to kernel cwd.
DATAFOLDER = (ROOT / "data").resolve()
FIGSFOLDER = (ROOT / "paper" / "figs" / "panels").resolve()
FIGSFOLDER.mkdir(parents=True, exist_ok=True)

SAVEFIGS = True

# %%
with open(DATAFOLDER / "assembled_data.pickle", "rb") as f:
    assembled = dill.load(f)

x_array = assembled["x_array"].copy()
print(f"Loaded x_array with {len(x_array)} rows and {x_array.shape[1]} columns")
print(f"Columns include: {', '.join(sorted(x_array.columns)[:12])} ...")

ids_to_exclude = ["PB27"]


# %%
def get_transition_subset(df, condition="deplete", infusion="45NaCl"):
    return (
        df.query("condition == @condition & infusiontype == @infusion")
        .sort_values(["id", "trial"])
        .copy()
    )


def summarize_transition_fits(fits_df, id_col="id"):
    if fits_df is None or fits_df.empty:
        return pd.DataFrame()

    summary = fits_df.copy()
    if "x0_orig" in summary.columns:
        summary["x0_orig"] = summary["x0_orig"].round(2)
    if "k" in summary.columns:
        summary["k"] = summary["k"].round(3)

    cols = [c for c in [id_col, "model", "x0_orig", "k", "r_squared", "is_valid", "note"] if c in summary.columns]
    return summary[cols].sort_values(id_col).reset_index(drop=True)


def apply_manual_exclusions(fits_df, fit_traces, ids_to_exclude, id_col="id"):
    if fits_df is None:
        fits_df = pd.DataFrame()
    if fit_traces is None:
        fit_traces = {}

    cleaned_ids = [
        str(v).strip()
        for v in (ids_to_exclude or [])
        if pd.notna(v) and str(v).strip()
    ]
    if not cleaned_ids:
        return fits_df, fit_traces

    exclude_set = set(cleaned_ids)
    fits_out = fits_df.copy()
    if not fits_out.empty and id_col in fits_out.columns:
        fits_out = fits_out[~fits_out[id_col].astype(str).isin(exclude_set)].copy()

    traces_out = {
        subject_id: trace
        for subject_id, trace in fit_traces.items()
        if str(subject_id) not in exclude_set
    }

    removed_from_fits = sorted(set(fits_df.get(id_col, pd.Series(dtype=object)).astype(str)) & exclude_set) if not fits_df.empty else []
    removed_from_traces = sorted(set(str(k) for k in fit_traces.keys()) & exclude_set)
    if removed_from_fits or removed_from_traces:
        print(
            f"Applied manual exclusions {sorted(exclude_set)} | "
            f"removed from fits: {removed_from_fits} | removed from traces: {removed_from_traces}"
        )
    else:
        print(f"Manual exclusions requested {sorted(exclude_set)} but no matching ids were found.")

    return fits_out, traces_out


# %%
def plot_subject_fit_traces(fit_traces, y_label, title_prefix, ncols=4, figsize_per_panel=(3.0, 2.4)):
    if not fit_traces:
        print("No fit traces to plot.")
        return None

    subject_ids = sorted(fit_traces.keys())
    n = len(subject_ids)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(figsize_per_panel[0] * ncols, figsize_per_panel[1] * nrows),
        squeeze=False,
    )

    for idx, subject_id in enumerate(subject_ids):
        ax = axes[idx // ncols][idx % ncols]
        trace = fit_traces[subject_id]

        x = trace.get("x", trace.get("x_fit"))
        y = trace.get("y", trace.get("y_fit_raw"))
        y_hat = trace.get("y_hat", trace.get("y_hat_raw"))

        if x is None or y is None:
            ax.axis("off")
            continue

        ax.scatter(x, y, s=16, alpha=0.6, color="#4C5B61", label="data")
        if y_hat is not None:
            ax.plot(x, y_hat, color="#D1495B", lw=1.8, label=trace.get("model", "fit"))
        if np.isfinite(trace.get("x0_orig", np.nan)):
            ax.axvline(trace["x0_orig"], ls="--", lw=1, color="#2F6690", alpha=0.7)

        ax.set_title(str(subject_id), fontsize=9)
        ax.set_xlabel("Trial")
        ax.set_ylabel(y_label)
        sns.despine(ax=ax, offset=3)

    for idx in range(n, nrows * ncols):
        axes[idx // ncols][idx % ncols].axis("off")

    handles, labels = axes[0][0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper right", frameon=False)
    fig.suptitle(title_prefix, y=1.02)
    fig.tight_layout()
    return fig


def _clean_ids(ids):
    return [str(v).strip() for v in (ids or []) if pd.notna(v) and str(v).strip()]


def summarize_failed_subjects(
    df,
    fits_df,
    signal_col,
    id_col="id",
    trial_col="trial",
    excluded_ids=None,
    value_transform=None,
    min_x0=0,
    max_x0=50,
    fit_kind="continuous",
):
    excluded_set = set(_clean_ids(excluded_ids))
    all_ids = [str(v) for v in sorted(df[id_col].dropna().unique())]
    valid_ids = set() if fits_df is None or fits_df.empty else set(fits_df[id_col].astype(str).tolist())

    rows = []
    for subject_id in all_ids:
        if subject_id in valid_ids or subject_id in excluded_set:
            continue

        sdf = df.loc[df[id_col].astype(str) == subject_id].sort_values(trial_col)
        y_raw = sdf[signal_col].to_numpy(dtype=float)
        y_for_fit = np.asarray(value_transform(y_raw), dtype=float) if value_transform is not None else y_raw
        x = np.arange(len(y_for_fit), dtype=float) if fit_kind == "logistic" else np.arange(1, len(y_for_fit) + 1, dtype=float)

        if fit_kind == "logistic":
            fit = fit_logistic_per_series(
                y_for_fit,
                x=x,
                prefer_4p=True,
                direction="decreasing",
                maxfev=60000,
            )
            success = bool(fit.get("success", False))
            note = str(fit.get("note", ""))
            x0 = fit.get("x0_orig", np.nan)
            if success and np.isfinite(x0):
                if (min_x0 is not None and x0 <= min_x0) or (max_x0 is not None and x0 >= max_x0):
                    note = (note + ";x0_out_of_bounds").strip(";")
        else:
            fit = fit_continuous_sigmoid_per_series(
                y_for_fit,
                x=x,
                maxfev=30000,
                require_negative_k=True,
            )
            success = bool(fit.get("success", False) and fit.get("is_valid", False))
            note = str(fit.get("note", ""))
            x0 = fit.get("x0_orig", np.nan)
            if fit.get("success", False) and np.isfinite(x0):
                if (min_x0 is not None and x0 <= min_x0) or (max_x0 is not None and x0 >= max_x0):
                    note = (note + ";x0_out_of_bounds").strip(";")

        rows.append({
            "id": subject_id,
            "reason": note if note else "not_selected",
            "n_trials": int(len(y_for_fit)),
            "n_finite": int(np.isfinite(y_for_fit).sum()),
        })

    failed_df = pd.DataFrame(rows).sort_values("id").reset_index(drop=True) if rows else pd.DataFrame(columns=["id", "reason", "n_trials", "n_finite"])
    return failed_df


def plot_failed_subject_series(
    df,
    failed_df,
    signal_col,
    title_prefix,
    id_col="id",
    trial_col="trial",
    value_transform=None,
    ncols=4,
    figsize_per_panel=(3.0, 2.2),
):
    if failed_df is None or failed_df.empty:
        print("No failed subjects to plot.")
        return None

    subject_ids = failed_df["id"].tolist()
    n = len(subject_ids)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(figsize_per_panel[0] * ncols, figsize_per_panel[1] * nrows),
        squeeze=False,
        sharex=False,
    )

    reason_map = failed_df.set_index("id")["reason"].to_dict()

    for idx, subject_id in enumerate(subject_ids):
        ax = axes[idx // ncols][idx % ncols]
        sdf = df.loc[df[id_col].astype(str) == str(subject_id)].sort_values(trial_col)
        x = np.arange(len(sdf), dtype=float)
        y = sdf[signal_col].to_numpy(dtype=float)
        if value_transform is not None:
            y = np.asarray(value_transform(y), dtype=float)

        ax.scatter(x, y, s=15, alpha=0.75, color="#6C757D")
        ax.plot(x, y, lw=1.2, alpha=0.7, color="#6C757D")
        ax.set_title(f"{subject_id}: {reason_map.get(subject_id, '')}", fontsize=8)
        ax.set_xlabel("Trial")
        ax.set_ylabel(signal_col)
        sns.despine(ax=ax, offset=3)

    for idx in range(n, nrows * ncols):
        axes[idx // ncols][idx % ncols].axis("off")

    fig.suptitle(title_prefix, y=1.02)
    fig.tight_layout()
    return fig


def _keep_complete_aligned_trials(subset, align_col, id_col="id"):
    if subset.empty:
        return subset

    work = subset.copy()
    work["_trial_aligned"] = work[align_col].round().astype(int)

    bounds = work.groupby(id_col)["_trial_aligned"].agg(["min", "max"])
    common_min = int(bounds["min"].max())
    common_max = int(bounds["max"].min())
    if common_min > common_max:
        return work.iloc[0:0].copy()

    expected = np.arange(common_min, common_max + 1, dtype=int)
    work = work.query("_trial_aligned >= @common_min and _trial_aligned <= @common_max").copy()

    valid_ids = []
    for rat_id, rat_df in work.groupby(id_col):
        observed = np.sort(rat_df["_trial_aligned"].unique())
        if np.array_equal(observed, expected):
            valid_ids.append(rat_id)

    return work[work[id_col].isin(valid_ids)].copy()


def _sigmoid_eval(x, params):
    L, k, x0, b = [float(v) for v in params]
    z = np.clip(-k * (x - x0), -60, 60)
    return L / (1 + np.exp(z)) + b


def plot_realignment_overlay(df, align_col, y_cols=("simba_median_balance", "auc_snips"), condition="deplete", infusion="45NaCl"):
    subset = df.query("condition == @condition & infusiontype == @infusion").dropna(subset=[align_col]).copy()
    subset = _keep_complete_aligned_trials(subset, align_col=align_col, id_col="id")

    if subset.empty:
        print(f"No complete aligned rows available for {align_col}")
        return None

    subset["trial_aligned"] = subset[align_col].round().astype(int)
    print(f"Using {subset['id'].nunique()} rats with complete aligned trial windows ({subset['trial_aligned'].min()} to {subset['trial_aligned'].max()})")

    fig, axes = plt.subplots(1, len(y_cols), figsize=(4.5 * len(y_cols), 3.2), sharex=True)
    if len(y_cols) == 1:
        axes = [axes]

    for ax, y_col in zip(axes, y_cols):
        grouped = subset.groupby("trial_aligned")[y_col].agg(["mean", "sem"]).reset_index()
        ax.plot(grouped["trial_aligned"], grouped["mean"], color="#D1495B", lw=2, label="mean")
        ax.fill_between(
            grouped["trial_aligned"],
            grouped["mean"] - grouped["sem"],
            grouped["mean"] + grouped["sem"],
            color="#D1495B",
            alpha=0.25,
        )

        x = grouped["trial_aligned"].to_numpy(dtype=float)
        y = grouped["mean"].to_numpy(dtype=float)
        sig_fit = fit_curve_series(x, y, model_name="sigmoidal", maxfev=30000)
        if sig_fit["success"] and np.all(np.isfinite(sig_fit["params"])):
            x_fit = np.linspace(x.min(), x.max(), 200)
            y_fit = _sigmoid_eval(x_fit, sig_fit["params"])
            k_val = float(sig_fit["params"][1])
            ax.plot(x_fit, y_fit, color="#2F6690", lw=2, ls="--", label="sigmoid fit")
            ax.text(
                0.02,
                0.95,
                f"k={k_val:.3f}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=9,
                color="#2F6690",
            )
            print(f"{align_col} | {y_col}: sigmoid k={k_val:.6f}")
        else:
            print(f"{align_col} | {y_col}: sigmoid fit failed")

        ax.axvline(0, color="#2F6690", ls="--", lw=1)
        ax.axhline(0, color="black", ls=":", lw=1, alpha=0.6)
        ax.set_title(y_col)
        ax.set_xlabel("Trials relative to transition")
        ax.set_ylabel("Mean ± SEM")
        sns.despine(ax=ax, offset=3)

    fig.suptitle(f"Realignment using {align_col}", y=1.03)
    fig.tight_layout()
    return fig


# %% [markdown]
# ### Raw median balance sigmoids

# %%
df_dep45 = get_transition_subset(x_array)

fits_raw, traces_raw = fit_continuous_sigmoid_transitions_by_id(
    df_dep45,
    signal_col="simba_median_balance",
    maxfev=30000,
    require_negative_k=True,
    min_x0=0,
    max_x0=50,
)

print(f"Raw simba_median_balance fits: {len(fits_raw)} valid subjects")
display(summarize_transition_fits(fits_raw))

fig_raw = plot_subject_fit_traces(
    traces_raw,
    y_label="simba_median_balance",
    title_prefix="Raw simba_median_balance sigmoid fits",
)
if savefigs and fig_raw is not None:
    fig_raw.savefig(FIGSFOLDER / "sigmoid_raw_simba_median_balance_subjects.png", dpi=300)

failed_raw = summarize_failed_subjects(
    df_dep45,
    fits_raw,
    signal_col="simba_median_balance",
    excluded_ids=ids_to_exclude if "ids_to_exclude" in globals() else None,
    min_x0=0,
    max_x0=50,
    fit_kind="continuous",
)
print(f"Raw simba_median_balance failed subjects: {len(failed_raw)}")
display(failed_raw)

fig_raw_failed = plot_failed_subject_series(
    df_dep45,
    failed_raw,
    signal_col="simba_median_balance",
    title_prefix="Raw simba_median_balance failed subjects",
)
if savefigs and fig_raw_failed is not None:
    fig_raw_failed.savefig(FIGSFOLDER / "sigmoid_raw_simba_median_balance_failed_subjects.png", dpi=300)

# %%
x_raw_aligned = x_array.copy()
x_raw_aligned["realigned_trials_raw"] = build_realigned_trials(
    x_raw_aligned,
    fits_raw,
    output_col="realigned_trials_raw",
)
x_raw_aligned["trial_aligned"] = x_raw_aligned["realigned_trials_raw"]

fig_raw_align = plot_realignment_overlay(
    x_raw_aligned,
    align_col="realigned_trials_raw",
    y_cols=("simba_median_balance", "auc_snips"),
)
if savefigs and fig_raw_align is not None:
    fig_raw_align.savefig(FIGSFOLDER / "realignment_raw_simba_median_balance.png", dpi=300)

# %%
# Sigmoidal overlays for the same 7 rats selected by median_balance fits
from matplotlib.dates import SA


if "fits_raw" not in globals() or fits_raw is None or fits_raw.empty:
    raise ValueError("fits_raw is not available. Run the raw median_balance fit cell first.")

required_rats = ["PB26", "PB30", "PB31", "PB46", "PB48", "PB73", "PB75"]

median_balance_fits = fits_raw.copy()
median_balance_fits["id"] = median_balance_fits["id"].astype(str)
median_balance_fits = median_balance_fits[median_balance_fits["id"].isin(required_rats)].copy()

included_ids = sorted(median_balance_fits["id"].unique())
missing_ids = sorted(set(required_rats) - set(included_ids))
if missing_ids:
    print(f"Warning: missing required rats in fits_raw: {missing_ids}")
if not included_ids:
    raise ValueError("No required rats found in fits_raw.")

x0_lookup = median_balance_fits.set_index("id")["x0_orig"].to_dict()
output_folder = (ROOT / "paper" / "figs" / "panels").resolve()
output_folder.mkdir(parents=True, exist_ok=True)

df_dep_45 = (
    x_array
    .query("condition == 'deplete' & infusiontype == '45NaCl'")
    .copy()
)
df_dep_45 = df_dep_45[df_dep_45["id"].astype(str).isin(included_ids)].copy()

all_fits = []
f, ax = plt.subplots(figsize=(6, 4),
                     gridspec_kw={"left": 0.2, "right": 0.95, "top": 0.9, "bottom": 0.2})

fit_line_color = COLORS[3]
transition_color = COLORS[3]

for rat in included_ids:
    sig = df_dep_45.loc[df_dep_45["id"].astype(str) == rat, "simba_median_balance"].to_numpy(dtype=float)
    y = sig
    # Match the trial indexing convention used by the validated raw-fit workflow above.
    x = np.arange(1, len(y) + 1, dtype=float)

    fit = fit_continuous_sigmoid_per_series(
        y,
        x=x,
        maxfev=30000,
        require_negative_k=True,
    )

    params = fit.get("params", {}) or {}
    x0_fit = float(fit.get("x0_orig", np.nan)) if np.isfinite(fit.get("x0_orig", np.nan)) else np.nan
    k_val = float(params.get("k", np.nan)) if np.isfinite(params.get("k", np.nan)) else np.nan
    success = bool(fit.get("success", False) and fit.get("is_valid", False))

    all_fits.append({
        "id": rat,
        "x0_median_balance": x0_lookup.get(rat, np.nan),
        "x0_sigmoidal": x0_fit,
        "k_sigmoidal": k_val,
        "model": "sigmoidal",
        "success": success,
        "note": fit.get("note", "") if success else (fit.get("note", "fit_failed_or_invalid") or "fit_failed_or_invalid"),
    })

    trace = fit.get("fit_trace", None)
    if success and trace is not None:
        x_fit = np.asarray(trace.get("x_fit", []), dtype=float)
        y_hat_raw = np.asarray(trace.get("y_hat_raw", []), dtype=float)
        if x_fit.size and y_hat_raw.size:
            ax.plot(x_fit, y_hat_raw, color=fit_line_color, linestyle="--", alpha=0.6)

fits_df = pd.DataFrame(all_fits).sort_values("id").reset_index(drop=True)

x0 = np.array([x0_lookup[r] for r in included_ids if np.isfinite(x0_lookup.get(r, np.nan))], dtype=float)
if x0.size:
    ax.plot(x0, [1.1] * len(x0), marker="o", linestyle="None", color=transition_color, alpha=0.7, clip_on=False)
    ax.text(np.max(x0) + 2, 1.1, "Median-balance transition points", ha="left", va="center", fontsize=10, color=transition_color)

    ax.plot([np.mean(x0), np.mean(x0)], [1.05, 1.15], color=transition_color, linestyle="--", alpha=0.7, clip_on=False)
    ax.text(np.mean(x0), 1.16, f"Mean=trial {int(np.mean(x0))}", ha="center", va="bottom", fontsize=10, color=transition_color)

sns.despine(ax=ax, offset=5)
ax.set_xlabel("Trial Number")
ax.set_ylabel("Apptetive probability")
ax.set_ylim(-1.1, 1.1)
ax.set_yticks([-1, 0, 1])
# ax.set_title(f"Sigmoidal fits using median_balance-selected rats (n={len(included_ids)})")

if SAVEFIGS:
    save_figure_atomic(
        f,
        "figSx_transitions_sigmoidal_fits_45NaCl_median_balance_selected_rats",
        output_folder,
        dpi=600,
    )

fits_df_euclidean = fits_df
print(f"Included rats ({len(included_ids)}): {included_ids}")
display(fits_df_euclidean)

# %%
df_dep_45.columns

# %% [markdown]
# ### Binarized median balance sigmoids
# 1 is median balance > 0.7, 0 is median balance <-0.7, values between -0.7 and 0.7 are ignored/not fitted

# %%
df_dep45 = get_transition_subset(x_array)

fits_bin, traces_bin = fit_logistic_transitions_by_id(
    df_dep45,
    signal_col="simba_median_balance",
    value_transform=lambda values: binarize_series_threshold(values, low=-0.7, high=0.7),
    direction="decreasing",
    maxfev=60000,
    min_x0=0,
    max_x0=50,
)

print(f"Binarized simba_median_balance logistic fits: {len(fits_bin)} valid subjects")
display(summarize_transition_fits(fits_bin))

fig_bin = plot_subject_fit_traces(
    traces_bin,
    y_label="binarized simba_median_balance",
    title_prefix="Binarized simba_median_balance logistic fits",
)
if savefigs and fig_bin is not None:
    fig_bin.savefig(FIGSFOLDER / "logistic_binarized_simba_median_balance_subjects.png", dpi=300)

failed_bin = summarize_failed_subjects(
    df_dep45,
    fits_bin,
    signal_col="simba_median_balance",
    value_transform=lambda values: binarize_series_threshold(values, low=-0.7, high=0.7),
    excluded_ids=ids_to_exclude if "ids_to_exclude" in globals() else None,
    min_x0=0,
    max_x0=50,
    fit_kind="logistic",
)
print(f"Binarized simba_median_balance failed subjects: {len(failed_bin)}")
display(failed_bin)

fig_bin_failed = plot_failed_subject_series(
    df_dep45,
    failed_bin,
    signal_col="simba_median_balance",
    value_transform=lambda values: binarize_series_threshold(values, low=-0.7, high=0.7),
    title_prefix="Binarized simba_median_balance failed subjects",
)
if savefigs and fig_bin_failed is not None:
    fig_bin_failed.savefig(FIGSFOLDER / "logistic_binarized_simba_median_balance_failed_subjects.png", dpi=300)

# %%
x_bin_aligned = x_array.copy()
x_bin_aligned["realigned_trials_bin"] = build_realigned_trials(
    x_bin_aligned,
    fits_bin,
    output_col="realigned_trials_bin",
)
x_bin_aligned["trial_aligned"] = x_bin_aligned["realigned_trials_bin"]

fig_bin_align = plot_realignment_overlay(
    x_bin_aligned,
    align_col="realigned_trials_bin",
    y_cols=("simba_median_balance", "auc_snips"),
)
if savefigs and fig_bin_align is not None:
    fig_bin_align.savefig(FIGSFOLDER / "realignment_binarized_simba_median_balance.png", dpi=300)

# %% [markdown]
# ### Quantized median balance (-1, 0, and +1)
# 1 is median balance > 0.7, -1 is median balance <-0.7, values between -0.7 and 0.7 are 0

# %%
df_dep45_quant = get_transition_subset(x_array).assign(
    simba_median_balance_quantized=lambda frame: quantize_series_threshold(frame["simba_median_balance"].to_numpy(), low=-0.7, high=0.7)
)

fits_quant, traces_quant = fit_continuous_sigmoid_transitions_by_id(
    df_dep45_quant,
    signal_col="simba_median_balance_quantized",
    maxfev=30000,
    require_negative_k=True,
    min_x0=0,
    max_x0=50,
)

print(f"Quantized simba_median_balance sigmoid fits: {len(fits_quant)} valid subjects")
display(summarize_transition_fits(fits_quant))

fig_quant = plot_subject_fit_traces(
    traces_quant,
    y_label="quantized simba_median_balance",
    title_prefix="Quantized (-1/0/1) simba_median_balance sigmoid fits",
)
if savefigs and fig_quant is not None:
    fig_quant.savefig(FIGSFOLDER / "sigmoid_quantized_simba_median_balance_subjects.png", dpi=300)

failed_quant = summarize_failed_subjects(
    df_dep45_quant,
    fits_quant,
    signal_col="simba_median_balance_quantized",
    excluded_ids=ids_to_exclude if "ids_to_exclude" in globals() else None,
    min_x0=0,
    max_x0=50,
    fit_kind="continuous",
)
print(f"Quantized simba_median_balance failed subjects: {len(failed_quant)}")
display(failed_quant)

fig_quant_failed = plot_failed_subject_series(
    df_dep45_quant,
    failed_quant,
    signal_col="simba_median_balance_quantized",
    title_prefix="Quantized simba_median_balance failed subjects",
)
if savefigs and fig_quant_failed is not None:
    fig_quant_failed.savefig(FIGSFOLDER / "sigmoid_quantized_simba_median_balance_failed_subjects.png", dpi=300)

# %%
x_quant_aligned = x_array.copy()
x_quant_aligned["realigned_trials_quant"] = build_realigned_trials(
    x_quant_aligned,
    fits_quant,
    output_col="realigned_trials_quant",
)
x_quant_aligned["trial_aligned"] = x_quant_aligned["realigned_trials_quant"]

fig_quant_align = plot_realignment_overlay(
    x_quant_aligned,
    align_col="realigned_trials_quant",
    y_cols=("simba_median_balance", "auc_snips"),
)
if savefigs and fig_quant_align is not None:
    fig_quant_align.savefig(FIGSFOLDER / "realignment_quantized_simba_median_balance.png", dpi=300)

# %%
transition_comparison = (
    fits_raw[["id", "x0_orig"]]
    .rename(columns={"x0_orig": "x0_raw_sigmoid"})
    .merge(
        fits_bin[["id", "x0_orig"]].rename(columns={"x0_orig": "x0_binarized_logistic"}),
        on="id",
        how="outer",
    )
    .merge(
        fits_quant[["id", "x0_orig"]].rename(columns={"x0_orig": "x0_quantized_sigmoid"}),
        on="id",
        how="outer",
    )
    .sort_values("id")
)

display(transition_comparison.round(2))

if not transition_comparison.empty:
    fig, ax = plt.subplots(figsize=(6.5, 4.0))
    for col, color in [
        ("x0_raw_sigmoid", "#2F6690"),
        ("x0_binarized_logistic", "#D1495B"),
        ("x0_quantized_sigmoid", "#EDAE49"),
    ]:
        vals = transition_comparison[col].dropna()
        ax.scatter(vals, np.repeat(col, len(vals)), alpha=0.8, color=color)

    ax.set_xlabel("Transition trial (x0)")
    ax.set_ylabel("Fit type")
    sns.despine(ax=ax, offset=3)
    fig.tight_layout()

# %%
# Optional: manually remove rats after inspecting fits (no refit needed).
# Example: ids_to_exclude = ["PB71", "PB48"]
ids_to_exclude = ["PB27"]

print(f"ids_to_exclude: {ids_to_exclude}")

if ids_to_exclude:
    if "fits_raw" in globals():
        fits_raw, traces_raw = apply_manual_exclusions(fits_raw, traces_raw, ids_to_exclude)
    if "fits_bin" in globals():
        fits_bin, traces_bin = apply_manual_exclusions(fits_bin, traces_bin, ids_to_exclude)
    if "fits_quant" in globals():
        fits_quant, traces_quant = apply_manual_exclusions(fits_quant, traces_quant, ids_to_exclude)

    print("Manual exclusions applied. Re-run realignment and comparison cells to refresh outputs.")
else:
    print("No manual exclusions applied.")

# %%
