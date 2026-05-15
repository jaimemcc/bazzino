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
# # Figure SX: Within-Trial Cross-Correlations
#
# This notebook computes per-trial cross-correlations between two snip arrays and returns one row per snip with:
# - `best_lag` (in bins)
# - `r` (Pearson correlation at that lag)
#
# Default lag window is `-30` to `+30` bins (for `0.1 s/bin`, that is `-3` to `+3 s`).

# %%
# %load_ext autoreload
# %autoreload 2

from pathlib import Path
import sys
sys.path.insert(0, str(Path("../src").resolve()))

from pickle_compat import enable_dill_pathlib_compat
enable_dill_pathlib_compat()

import numpy as np
import pandas as pd
import dill

from figure_config import DATAFOLDER, RESULTSFOLDER

assembled_data_path = DATAFOLDER / "assembled_data.pickle"
with open(assembled_data_path, "rb") as f:
    data = dill.load(f)

snips_photo = np.asarray(data["snips_photo"], dtype=float)
snips_simba_shifted_baseline = np.asarray(data["snips_simba_shifted_baseline"], dtype=float)

print(f"Loaded: {assembled_data_path}")
print(f"snips_photo shape: {snips_photo.shape}")
print(f"snips_simba_shifted_baseline shape: {snips_simba_shifted_baseline.shape}")


# %%
def _lagged_pearson(x, y, lag, min_overlap=10, min_std=1e-12):
    """Pearson r between 1D vectors at a given lag.

    Convention: positive lag means y is shifted right relative to x.
    """
    if lag > 0:
        x_use = x[:-lag]
        y_use = y[lag:]
    elif lag < 0:
        k = -lag
        x_use = x[k:]
        y_use = y[:-k]
    else:
        x_use = x
        y_use = y

    valid = np.isfinite(x_use) & np.isfinite(y_use)
    x_use = x_use[valid]
    y_use = y_use[valid]

    if x_use.size < min_overlap:
        return np.nan
    if np.std(x_use) <= min_std or np.std(y_use) <= min_std:
        return np.nan

    return float(np.corrcoef(x_use, y_use)[0, 1])


def crosscorr_best_lag_df(
    snips_x,
    snips_y,
    lag_min=-30,
    lag_max=30,
    min_overlap=10,
    min_std=1e-12,
    select_by="max_abs",
):
    """Return dataframe with per-snip best_lag and r.

    Inputs must be arrays shaped (n_snips, n_bins), with matched shape.
    """
    x = np.asarray(snips_x, dtype=float)
    y = np.asarray(snips_y, dtype=float)

    if x.ndim != 2 or y.ndim != 2:
        raise ValueError("Both inputs must be 2D arrays: (n_snips, n_bins).")
    if x.shape != y.shape:
        raise ValueError(f"Shape mismatch: x={x.shape}, y={y.shape}.")
    if not (isinstance(lag_min, int) and isinstance(lag_max, int) and lag_min <= lag_max):
        raise ValueError("lag_min and lag_max must be ints with lag_min <= lag_max.")
    if min_overlap < 3:
        raise ValueError("min_overlap must be >= 3.")
    if select_by not in {"max_abs", "max_pos"}:
        raise ValueError("select_by must be 'max_abs' or 'max_pos'.")

    n_snips, n_bins = x.shape
    lags = np.arange(lag_min, lag_max + 1, dtype=int)

    if np.max(np.abs(lags)) >= n_bins:
        raise ValueError(
            f"Lag range too wide for snip length {n_bins}. "
            f"Use absolute lag <= {n_bins - 1}."
        )

    best_lag = np.full(n_snips, np.nan)
    best_r = np.full(n_snips, np.nan)

    for i in range(n_snips):
        rs = np.array(
            [_lagged_pearson(x[i], y[i], lag, min_overlap=min_overlap, min_std=min_std) for lag in lags],
            dtype=float,
        )
        valid = np.isfinite(rs)
        if not np.any(valid):
            continue

        if select_by == "max_abs":
            score = np.abs(rs)
        else:
            score = rs

        max_score = np.nanmax(score)
        candidates = np.flatnonzero(np.isfinite(score) & np.isclose(score, max_score, rtol=1e-12, atol=1e-15))
        if candidates.size > 1:
            chosen = candidates[np.argmin(np.abs(lags[candidates]))]
        else:
            chosen = int(candidates[0])

        best_lag[i] = lags[chosen]
        best_r[i] = rs[chosen]

    out = pd.DataFrame({"best_lag": best_lag, "r": best_r})
    out["best_lag"] = out["best_lag"].astype("Int64")
    return out


# %%
# Controls
lag_min = -30
lag_max = 30
bin_size_s = 0.1
min_overlap = 10
select_by = "max_abs"  # 'max_abs' or 'max_pos'

xcorr_df = crosscorr_best_lag_df(
    snips_x=snips_photo,
    snips_y=snips_simba_shifted_baseline,
    lag_min=lag_min,
    lag_max=lag_max,
    min_overlap=min_overlap,
    select_by=select_by,
)

print(f"Result shape: {xcorr_df.shape}")
print(f"Rows with valid results: {xcorr_df['r'].notna().sum()} / {len(xcorr_df)}")
print(
    f"Lag window: [{lag_min}, {lag_max}] bins "
    f"(~[{lag_min * bin_size_s:.1f}, {lag_max * bin_size_s:.1f}] s at {bin_size_s:.3f} s/bin)"
)

xcorr_df.head()

# %%
# Build plotting dataframe by attaching trial metadata to xcorr_df
import matplotlib.pyplot as plt


def _resolve_column(df, candidates, required=True):
    for c in candidates:
        if c in df.columns:
            return c
    if required:
        raise KeyError(f"None of these columns were found: {candidates}")
    return None


if "x_array" not in data:
    raise KeyError("Expected key 'x_array' in assembled data for trial metadata.")

meta_df = data["x_array"].copy().reset_index(drop=True)

rat_col = _resolve_column(meta_df, ["id", "rat", "animal", "subject"])
condition_col = _resolve_column(meta_df, ["condition"])
infusion_col = _resolve_column(meta_df, ["infusiontype", "infusion_type", "infusion"])
trial_col = _resolve_column(meta_df, ["trial", "trial_num", "trial_number", "trial_index"])

if len(meta_df) != len(xcorr_df):
    raise ValueError(
        f"Length mismatch between metadata ({len(meta_df)}) and xcorr_df ({len(xcorr_df)})."
    )

plot_df = pd.concat([meta_df[[rat_col, condition_col, infusion_col, trial_col]], xcorr_df], axis=1)
plot_df = plot_df.rename(
    columns={
        rat_col: "rat",
        condition_col: "condition",
        infusion_col: "infusiontype",
        trial_col: "trial",
    }
)

# Numeric trial index (robust if trial labels are strings)
plot_df["trial"] = pd.to_numeric(plot_df["trial"], errors="coerce")
plot_df["best_lag_s"] = plot_df["best_lag"].astype(float) * float(bin_size_s)
plot_df["combo"] = plot_df["condition"].astype(str) + " | " + plot_df["infusiontype"].astype(str)

print(f"plot_df shape: {plot_df.shape}")
print(f"Unique combos: {plot_df['combo'].nunique()}")
plot_df[["rat", "condition", "infusiontype", "trial", "best_lag", "best_lag_s", "r"]].head()

# %%
# (1) Group by rat (average across trials), separate subplots by condition+infusion
rat_avg = (
    plot_df
    .dropna(subset=["r", "best_lag_s"])
    .groupby(["condition", "infusiontype", "combo", "rat"], as_index=False)
    .agg(r=("r", "mean"), best_lag_s=("best_lag_s", "mean"))
)

combos = sorted(rat_avg["combo"].unique())
if len(combos) == 0:
    raise ValueError("No valid rows available to plot rat averages.")

fig, axes = plt.subplots(2, len(combos), figsize=(4.2 * len(combos), 6), sharey="row")
if len(combos) == 1:
    axes = np.array(axes).reshape(2, 1)

for j, combo in enumerate(combos):
    sub = rat_avg[rat_avg["combo"] == combo].sort_values("rat")

    axes[0, j].plot(sub["rat"].astype(str), sub["r"], "o-", alpha=0.9)
    axes[0, j].axhline(0, color="k", lw=1, ls="--", alpha=0.5)
    axes[0, j].set_title(combo)
    axes[0, j].tick_params(axis="x", rotation=45)

    axes[1, j].plot(sub["rat"].astype(str), sub["best_lag_s"], "o-", alpha=0.9, color="tab:orange")
    axes[1, j].axhline(0, color="k", lw=1, ls="--", alpha=0.5)
    axes[1, j].tick_params(axis="x", rotation=45)
    axes[1, j].set_xlabel("Rat")

axes[0, 0].set_ylabel("Mean r (across trials)")
axes[1, 0].set_ylabel("Mean best lag (s)")
fig.suptitle("Cross-correlation summary grouped by rat", y=1.02)
fig.tight_layout()
plt.show()

# %%
# (2) Group by trial (average across rats), separate subplots by condition+infusion
trial_avg = (
    plot_df
    .dropna(subset=["trial", "r", "best_lag_s"])
    .groupby(["condition", "infusiontype", "combo", "trial"], as_index=False)
    .agg(r=("r", "mean"), best_lag_s=("best_lag_s", "mean"))
    .sort_values(["combo", "trial"])
)

combos = sorted(trial_avg["combo"].unique())
if len(combos) == 0:
    raise ValueError("No valid rows available to plot trial averages.")

fig, axes = plt.subplots(2, len(combos), figsize=(4.2 * len(combos), 6), sharey="row")
if len(combos) == 1:
    axes = np.array(axes).reshape(2, 1)

for j, combo in enumerate(combos):
    sub = trial_avg[trial_avg["combo"] == combo]

    axes[0, j].plot(sub["trial"], sub["r"], "o-", alpha=0.9)
    axes[0, j].axhline(0, color="k", lw=1, ls="--", alpha=0.5)
    axes[0, j].set_title(combo)

    axes[1, j].plot(sub["trial"], sub["best_lag_s"], "o-", alpha=0.9, color="tab:orange")
    axes[1, j].axhline(0, color="k", lw=1, ls="--", alpha=0.5)
    axes[1, j].set_xlabel("Trial")

axes[0, 0].set_ylabel("Mean r (across rats)")
axes[1, 0].set_ylabel("Mean best lag (s)")
fig.suptitle("Cross-correlation summary grouped by trial", y=1.02)
fig.tight_layout()
plt.show()

# %%
# Optional: save to CSV
save_csv = False
csv_path = RESULTSFOLDER / "figure_Sx_withintrial_crosscorr_best_lag_r.csv"

if save_csv:
    xcorr_df.to_csv(csv_path, index=False)
    print(f"Saved: {csv_path}")
else:
    print("CSV not saved (set save_csv = True to save).")
