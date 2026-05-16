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
import sys

# Register dill/pathlib compatibility shim BEFORE importing dill
sys.path.insert(0, str(Path("../src").resolve()))
from pickle_compat import enable_dill_pathlib_compat
enable_dill_pathlib_compat()

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
import seaborn as sns
import trompy as tp

from scipy import stats
from scipy.spatial.distance import cdist
from scipy.optimize import curve_fit

import dill

rcParams['font.family'] = 'Arial'
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['savefig.transparent'] = True
colors = ["#67AFD2", "#016895", "#F4795B", "#C74632"]

savefigs = False

DATAFOLDER = Path("..//data")
RESULTSFOLDER = Path("..//results")
FIGSFOLDER = Path("C:/Users/jmc010/Dropbox/Publications in Progress/Bazzino Roitman_sodium/figs")

# %% [markdown]
# ## Load Assembled Data
#
# Load the complete dataset from assembled_data.pickle which includes trial metadata, photometry snips, and cluster assignments.

# %%
assembled_data_path = DATAFOLDER / "assembled_data.pickle"

with open(assembled_data_path, "rb") as f:
    data = dill.load(f)

# Extract main components
x_array = data["x_array"]
snips_photo = data["snips_photo"]
snips_simba = data["snips_simba"]

fits_df = data["fits_df"]
params = data.get("params", {})

print(f"Loaded assembled data from {assembled_data_path}")
print(f"  - x_array shape: {x_array.shape}")
print(f"  - snips_photo shape: {snips_photo.shape}")
print(f"  - Number of unique animals: {x_array['id'].nunique()}")

# %%
# Load PCA-transformed photometry data (optional - if needed for PC-based analysis)
pcafile = RESULTSFOLDER / "transformed_data_photo.pickle"

if pcafile.exists():
    with open(pcafile, 'rb') as f:
        pca = dill.load(f)
    pc1 = pca[:,0]
    pca_data = pca[:, :3]
    print(f"Loaded PCA data: {pca_data.shape}")
else:
    print(f"Warning: PCA file not found at {pcafile}")
    pca_data = None

# %% [markdown]
# ## Calculate trial distances based on clusters and PCs

# %%
# Note: assembled_data.pickle should already contain cluster_photo assignments
# This cell calculates continuous distance metrics based on PCA space

if pca_data is None:
    print("Skipping clusterness calculation - PCA data not available")
elif 'cluster_photo' not in x_array.columns:
    print("Warning: cluster_photo column not found in x_array")
    print("Cannot calculate clusterness without cluster assignments")
else:
    # calculate centroids
    cluster_0_centroid = pca_data[x_array.cluster_photo == 0].mean(axis=0)
    cluster_1_centroid = pca_data[x_array.cluster_photo == 1].mean(axis=0)

    ## First, to work out projections
    # Step 2: Define the cluster separation vector
    cluster_vector = cluster_0_centroid - cluster_1_centroid

    # Step 3: Project each observation onto the cluster vector
    # Normalize the cluster vector
    cluster_vector_norm = cluster_vector / np.linalg.norm(cluster_vector)
    # Compute projections
    projections = np.dot(pca_data - cluster_1_centroid, cluster_vector_norm)

    # Step 4: Normalize the projections to range between 0 and 1
    min_projection = projections.min()
    max_projection = projections.max()
    normalized_projections = (projections - min_projection) / (max_projection - min_projection)

    x_array = x_array.assign(clusterness_photo=normalized_projections)


    ## Second to work out Euclidian distances
    # Stack centroids into a matrix
    centroids = np.vstack([cluster_0_centroid, cluster_1_centroid])

    # Calculate all distances at once using cdist
    # This creates a matrix where each row is an observation and each column is a centroid
    distances = cdist(pca_data, centroids, metric='euclidean')

    distances_diff = distances[:, 1] - distances[:, 0]

    x_array = x_array.assign(euclidean_diff=distances_diff)
    
    print("Added clusterness_photo and euclidean_diff columns to x_array")

# %%
## to save data with clutserness and euclidian diff

with open(DATAFOLDER / "bazzino_data_with_clusters_and_dists.pickle", "wb") as f:
    dill.dump({
        "x_array": x_array,
        "snips_photo": snips_photo,
        "snips_vel": snips_vel,
        "pca": pca
    }, f)

# %%
# to show eligible rats for this fitting
x_array.query("condition == 'deplete' & infusiontype == '45NaCl'").id.unique()


# %%
# Robust logistic fits for binary/near-binary data (handles optional offset)

# 3-parameter logistic (no baseline offset)
def logistic3(x, L, x0, k):
    return L / (1 + np.exp(-k * (x - x0)))

# 4-parameter logistic (baseline offset A)
def logistic4(x, A, L, x0, k):
    return A + (L - A) / (1 + np.exp(-k * (x - x0)))


def _normalize_x(x):
    x = np.asarray(x, dtype=float)
    m, s = float(np.mean(x)), float(np.std(x))
    if not np.isfinite(s) or s == 0:
        s = 1.0
    return (x - m) / s, m, s


def _clip_y(y):
    y = np.asarray(y, dtype=float)
    return np.clip(y, 1e-4, 1 - 1e-4)


def fit_logistic_per_series(y, x=None, prefer_4p=True, direction=None, maxfev=60000):
    """
    Fit a logistic curve to binary/near-binary data with robust inits and bounds.
    - prefer_4p: try 4-parameter (with baseline) first, then fallback to 3-parameter
    - direction: None to infer from corr(x,y); 'increasing' or 'decreasing' to enforce k sign
    Returns dict with keys: {'model','params','y_hat','x0_orig','success','note'}
    """
    if x is None:
        x = np.arange(len(y), dtype=float)
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    # Normalize x for stabler k/x0 estimation
    x_norm, x_mean, x_std = _normalize_x(x)
    y_clip = _clip_y(y)

    # Initial guesses from data
    y_min, y_max = float(np.min(y_clip)), float(np.max(y_clip))
    A_init = y_min
    L_init = y_max
    x0_init = 0.0

    # k sign from direction or correlation
    if direction is None:
        try:
            c = float(np.corrcoef(x, y_clip)[0, 1])
        except Exception:
            c = 0.0
        if not np.isfinite(c):
            c = 0.0
        sign = 1.0 if c >= 0 else -1.0
    else:
        sign = 1.0 if direction == 'increasing' else -1.0

    k_mags = [0.5, 1.0, 2.0]

    def try_fit(func, p0_list, bounds, n_params):
        best = None
        best_rss = np.inf
        for p0 in p0_list:
            try:
                popt, _ = curve_fit(func, x_norm, y_clip, p0=p0, bounds=bounds, maxfev=maxfev)
                y_hat = func(x_norm, *popt)
                rss = float(np.sum((y_clip - y_hat) ** 2))
                if rss < best_rss:
                    best_rss, best = rss, (popt, y_hat)
            except Exception:
                continue
        return best

    # 4-parameter attempt
    res4 = None
    if prefer_4p:
        p0s_4 = [[A_init, L_init, x0_init, sign * km] for km in k_mags]
        # Keep values sensible for binary-ish data
        bnds_4 = ([ -0.1,  0.4, -3.0, -10.0],
                  [  0.6,  1.6,  3.0,  10.0])
        res4 = try_fit(logistic4, p0s_4, bnds_4, 4)

    # 3-parameter fallback or primary
    res3 = None
    p0s_3 = [[L_init, x0_init, sign * km] for km in k_mags]
    bnds_3 = ([0.4, -3.0, -10.0], [1.6, 3.0, 10.0])
    res3 = try_fit(logistic3, p0s_3, bnds_3, 3)

    # Choose result: prefer successful 4p; otherwise 3p
    if res4 is not None:
        popt, y_hat = res4
        # x0 already in normalized units; convert to original units for reporting
        A, L, x0n, k = map(float, popt)
        x0_orig = x0n * x_std + x_mean
        return {
            'model': 'logistic4',
            'params': {'A': A, 'L': L, 'x0_norm': x0n, 'x0_orig': x0_orig, 'k': k},
            'y_hat': y_hat,
            'x0_orig': x0_orig,
            'success': True,
            'note': ''
        }
    elif res3 is not None:
        popt, y_hat = res3
        L, x0n, k = map(float, popt)
        x0_orig = x0n * x_std + x_mean
        return {
            'model': 'logistic3',
            'params': {'L': L, 'x0_norm': x0n, 'x0_orig': x0_orig, 'k': k},
            'y_hat': y_hat,
            'x0_orig': x0_orig,
            'success': True,
            'note': '4p failed; used 3p'
        }
    else:
        return {'model': None, 'params': {}, 'y_hat': None, 'x0_orig': np.nan, 'success': False, 'note': 'fit failed'}
    


# %%
# Logistic fits for raw cluster assignments (binary/inverted cluster_photo)
df_dep_45 = x_array.query("condition == 'deplete' & infusiontype == '45NaCl'").copy()
all_fits = []

f, ax = plt.subplots(figsize=(6,4))

for rat in df_dep_45.id.unique():
    sig = df_dep_45.loc[df_dep_45.id == rat, 'cluster_photo'].to_numpy()
    # Make binary: invert if needed (original code inverted)
    y = np.logical_not(sig).astype(int)
    x = np.arange(len(y), dtype=float)

    fit = fit_logistic_per_series(y, x=x, prefer_4p=True, direction='decreasing')
    all_fits.append({ 'id': rat, **fit['params'], 'model': fit['model'], 'x0_orig': fit['x0_orig'], 'success': fit['success'], 'note': fit['note'] })

    if fit["success"] and fit['x0_orig'] > 0 and fit['x0_orig'] < len(y):
        ax.plot(x, fit["y_hat"], color=colors[2], alpha=0.5, linestyle="--")

fits_df = pd.DataFrame(all_fits)
fits_df = fits_df.query("success == True and x0_orig > 0").copy()

x0 = fits_df['x0_orig'].to_list()
ax.plot(x0, [1.1]*len(x0), marker="o", linestyle="None", color=colors[2], alpha=0.5, clip_on=False)
ax.text(np.max(x0)+2, 1.1, "Transition points", ha="left", va="center", fontsize=10, color=colors[2])

ax.plot([np.mean(x0), np.mean(x0)], [1.05, 1.15], color=colors[2], linestyle="--", alpha=0.5, clip_on=False)
ax.text(np.mean(x0), 1.16, f"Mean=trial {int(np.mean(x0))}", ha="center", va="bottom", fontsize=10, color=colors[2])

sns.despine(ax=ax, offset=5)
ax.set_xlabel("Trial Number")
ax.set_ylabel("Probability of Cluster 1")

ax.set_yticks([0, 0.5, 1])
ax.set_ylim([-0.02, 1.1])

if savefigs:
    f.savefig(FIGSFOLDER / "logistic_fits_45NaCl.pdf", dpi=600, transparent=True)
    
fits_df_cluster_raw = fits_df
        
# print(fits_df)

# %%
# Logistic fits for clusterness (continuous between 0 and 1, based on projection onto cluster vector)
df_dep_45 = x_array.query("condition == 'deplete' & infusiontype == '45NaCl'").copy()
all_fits = []

f, ax = plt.subplots(figsize=(6,4))

for rat in df_dep_45.id.unique():
    sig = df_dep_45.loc[df_dep_45.id == rat, 'clusterness_photo'].to_numpy()
    # Make binary: invert if needed (original code inverted)
    y = sig
    x = np.arange(len(y), dtype=float)

    fit = fit_logistic_per_series(y, x=x, prefer_4p=True, direction='decreasing')
    all_fits.append({ 'id': rat, **fit['params'], 'model': fit['model'], 'x0_orig': fit['x0_orig'], 'success': fit['success'], 'note': fit['note'] })

    if fit["success"] and fit['x0_orig'] > 0 and fit['x0_orig'] < len(y):
        ax.plot(x, fit["y_hat"], color=colors[2], alpha=0.5, linestyle="--")

fits_df = pd.DataFrame(all_fits)
fits_df = fits_df.query("success == True and x0_orig > 0 and x0_orig < 50").copy()

x0 = fits_df['x0_orig'].to_list()
ax.plot(x0, [0.75]*len(x0), marker="o", linestyle="None", color=colors[2], alpha=0.5, clip_on=False)
ax.text(np.max(x0)+2, 0.75, "Transition points", ha="left", va="center", fontsize=10, color=colors[2])

ax.plot([np.mean(x0), np.mean(x0)], [0.73, 0.77], color=colors[2], linestyle="--", alpha=0.5, clip_on=False)
ax.text(np.mean(x0), 0.78, f"Mean=trial {int(np.mean(x0))}", ha="center", va="bottom", fontsize=10, color=colors[2])

sns.despine(ax=ax, offset=5)
ax.set_xlabel("Trial Number")
ax.set_ylabel("Cluster 1-ness")

# ax.set_yticks([0, 0.5, 1])
ax.set_ylim([0.31, 0.7])

if savefigs:
    f.savefig(FIGSFOLDER / "logistic_fits_45NaCl.pdf", dpi=600, transparent=True)
    
fits_df_clusterness = fits_df
        
# print(fits_df)


# %%
df_dep_45.columns

# %%
# Logistic fits for euclidean distance difference (continuous between -inf and +inf)
df_dep_45 = x_array.query("condition == 'deplete' & infusiontype == '45NaCl'").copy()
all_fits = []

f, ax = plt.subplots(figsize=(6,4))

for rat in df_dep_45.id.unique():
    sig = df_dep_45.loc[df_dep_45.id == rat, 'euclidean_diff'].to_numpy()
    # Make binary: invert if needed (original code inverted)
    y = sig
    x = np.arange(len(y), dtype=float)

    fit = fit_logistic_per_series(y, x=x, prefer_4p=True, direction='decreasing')
    all_fits.append({ 'id': rat, **fit['params'], 'model': fit['model'], 'x0_orig': fit['x0_orig'], 'success': fit['success'], 'note': fit['note'] })

    if fit["success"] and fit['x0_orig'] > 0 and fit['x0_orig'] < len(y):
        ax.plot(x, fit["y_hat"], color=colors[2], linestyle="--", alpha=0.5)

fits_df = pd.DataFrame(all_fits)
fits_df = fits_df.query("success == True and x0_orig > 0 and x0_orig < 50").copy()

x0 = fits_df['x0_orig'].to_list()
ax.plot(x0, [1.1]*len(x0), marker="o", linestyle="None", color=colors[2], alpha=0.5, clip_on=False)
ax.text(np.max(x0)+2, 1.1, "Transition points", ha="left", va="center", fontsize=10, color=colors[2])

ax.plot([np.mean(x0), np.mean(x0)], [1.05, 1.15], color=colors[2], linestyle="--", alpha=0.5, clip_on=False)
ax.text(np.mean(x0), 1.16, f"Mean=trial {int(np.mean(x0))}", ha="center", va="bottom", fontsize=10, color=colors[2])

sns.despine(ax=ax, offset=5)
ax.set_xlabel("Trial Number")
ax.set_ylabel("Cluster 1-ness")

# ax.set_yticks([0, 0.5, 1])
# ax.set_ylim([0.31, 0.7])

if savefigs:
    f.savefig(FIGSFOLDER / "logistic_fits_45NaCl.pdf", dpi=600, transparent=True)

fits_df_euclidean = fits_df
# print(fits_df)


# %%
with open(DATAFOLDER / "sigmoidal_fits.pickle", "wb") as f:
    dill.dump({
        "fits_df_cluster_raw": fits_df_cluster_raw,
        "fits_df_clusterness": fits_df_clusterness,
        "fits_df_euclidean": fits_df_euclidean
    }, f)

# %%
df_dep_45.id.unique()

# %%
## Quick visualization for one rat
example_rat = df_dep_45.id.unique()[2]
sig = df_dep_45.loc[df_dep_45.id == example_rat, 'cluster_photo'].to_numpy()
y = np.logical_not(sig).astype(int)
x = np.arange(len(y), dtype=float)
res = fit_logistic_per_series(y, x=x, prefer_4p=True, direction='decreasing')

plt.figure(figsize=(6, 3.5))
plt.scatter(x, y, s=24, color='#016895', alpha=0.8, label='Binary data')
if res['y_hat'] is not None:
    plt.plot(x, res['y_hat'], color='#F4795B', lw=2, label=f"{res['model']} fit")
    if np.isfinite(res['x0_orig']):
        plt.axvline(res['x0_orig'], color='#999', ls='--', lw=1)
        plt.text(res['x0_orig'], 1.02, f"x0≈{res['x0_orig']:.1f}", ha='center', va='bottom', fontsize=9)
plt.ylim(-0.1, 1.1)
plt.xlabel('Trial')
plt.ylabel('Binary outcome')
sns.despine()
plt.legend(frameon=False)
plt.tight_layout()
plt.show()

# %%
merged_fits = (
    fits_df_cluster_raw[['id', 'x0_orig']]
    .rename(columns={'x0_orig': 'x0_cluster_raw'})
    .merge(
        fits_df_clusterness[['id', 'x0_orig']].rename(columns={'x0_orig': 'x0_clusterness'}),
        on='id',
        how='outer'
    )
    .merge(
        fits_df_euclidean[['id', 'x0_orig']].rename(columns={'x0_orig': 'x0_euclidean'}),
        on='id',
        how='outer'
    )
    .dropna()
    .assign(
        x0_cluster_raw=lambda df: df['x0_cluster_raw'].round().astype(int),
        x0_clusterness=lambda df: df['x0_clusterness'].round().astype(int),
        x0_euclidean=lambda df: df['x0_euclidean'].round().astype(int)
    )
)

merged_fits

# %%
## Sigmoidal Fits for time_moving (Deplete + 45NaCl)
# Fit sigmoids to normalized time_moving data for each animal with quality checks from figure_1

from scipy.stats import pearsonr
import matplotlib.pyplot as plt

def _sigmoid_model(x, L, k, x0, b):
    """4-parameter sigmoid function."""
    z = np.clip(-k * (x - x0), -60, 60)
    return L / (1 + np.exp(z)) + b

def _safe_pearson(y_true, y_pred):
    """Safe Pearson correlation handling edge cases."""
    if np.allclose(np.std(y_true), 0) or np.allclose(np.std(y_pred), 0):
        return np.nan, np.nan
    return pearsonr(y_true, y_pred)

def _model_metrics(y_true, y_pred, n_params):
    """Calculate model fit metrics (RMSE, AIC, AICc, BIC)."""
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    n = len(y_true)
    residuals = y_true - y_pred
    sse = np.nansum(residuals ** 2)
    sse = max(float(sse), 1e-12)
    rmse = np.sqrt(sse / n) if n > 0 else np.nan
    aic = n * np.log(sse / n) + 2 * n_params if n > 0 else np.nan
    bic = n * np.log(sse / n) + n_params * np.log(n) if n > 1 else np.nan
    aicc = aic + (2 * n_params * (n_params + 1)) / (n - n_params - 1) if n > (n_params + 1) else np.nan
    return rmse, aic, aicc, bic

def _sigmoid_quality_checks(x, params, pcov):
    """
    Check sigmoid fit quality using criteria from figure_1_behaviour notebook:
    1. x0 must be interior (not near edges)
    2. k must be plausible (0.02 <= |k| <= 2.5)
    3. Confidence intervals must be finite
    4. Both asymptotes must be visible in data range
    """
    if params is None or len(params) != 4 or np.any(~np.isfinite(params)):
        return False, "fit_failed", {"x0_interior": False, "k_plausible": False, "ci_finite": False, "asymptotes_covered": False}

    L, k, x0, b = params
    x_min, x_max = float(np.min(x)), float(np.max(x))
    x_range = max(x_max - x_min, 1.0)
    edge_margin = 0.15 * x_range
    x0_interior = (x_min + edge_margin) <= x0 <= (x_max - edge_margin)
    k_plausible = np.isfinite(k) and (0.02 <= abs(k) <= 2.5)

    ci_finite = False
    if pcov is not None:
        diag = np.diag(pcov)
        if np.all(np.isfinite(diag)) and np.all(diag >= 0):
            se = np.sqrt(diag)
            ci_finite = np.all(np.isfinite(se)) and np.all(se > 0)

    amplitude = abs(L)
    if amplitude < 1e-8:
        asymptotes_covered = False
    else:
        lower_asym = min(b, b + L)
        upper_asym = max(b, b + L)
        y_start = _sigmoid_model(np.array([x_min]), L, k, x0, b)[0]
        y_end = _sigmoid_model(np.array([x_max]), L, k, x0, b)[0]

        tol = 0.2
        start_low = abs(y_start - lower_asym) / amplitude <= tol
        start_high = abs(y_start - upper_asym) / amplitude <= tol
        end_low = abs(y_end - lower_asym) / amplitude <= tol
        end_high = abs(y_end - upper_asym) / amplitude <= tol

        start_side = "low" if start_low else ("high" if start_high else None)
        end_side = "low" if end_low else ("high" if end_high else None)
        asymptotes_covered = (start_side is not None) and (end_side is not None) and (start_side != end_side)

    checks = {
        "x0_interior": bool(x0_interior),
        "k_plausible": bool(k_plausible),
        "ci_finite": bool(ci_finite),
        "asymptotes_covered": bool(asymptotes_covered),
    }

    failed = [name for name, ok in checks.items() if not ok]
    is_valid = len(failed) == 0
    reasons = "ok" if is_valid else ";".join(failed)
    return is_valid, reasons, checks

# Get deplete + 45NaCl subset
deplete_45_subset = (
    x_array
    .query("condition == 'deplete' & infusiontype == '45NaCl'")
    .copy()
)

# Get unique animals
animals = sorted(deplete_45_subset['id'].unique())
print(f"Fitting sigmoids for {len(animals)} animals in deplete + 45NaCl condition")
print(f"Animals: {animals}\n")

# Fit sigmoids for each animal
fit_results = []
fit_traces = []

for animal in animals:
    # Get this animal's data
    animal_data = (
        deplete_45_subset
        .query("id == @animal")
        .sort_values('trial')
        .copy()
    )
    
    # Extract auc_simba values
    y_raw = animal_data.median_balance.values
    x = np.arange(1, len(y_raw) + 1, dtype=float)
    
    # Normalize to [0, 1]
    y_min, y_max = y_raw.min(), y_raw.max()
    y_range = y_max - y_min
    if y_range < 1e-8:
        print(f"  {animal}: SKIPPED (no variation in auc_simba)")
        continue
    
    y = (y_raw - y_min) / y_range
    
    # Determine initial k value based on expected direction
    # For time_moving, we expect increase (aversion develops), so positive k
    initial_k = -1.0
    
    # Set up sigmoid parameters for time_moving (expected to increase)
    p0 = [
        1.0,                    # L: amplitude (normalized to 1)
        initial_k,              # k: steepness (positive for increase)
        np.median(x),           # x0: midpoint
        0.0                     # b: baseline (normalized to 0)
    ]
    bounds = (
        [-np.inf, -np.inf, np.min(x), -np.inf],
        [np.inf, np.inf, np.max(x), np.inf]
    )
    
    # Fit sigmoid
    try:
        params, pcov = curve_fit(_sigmoid_model, x, y, p0=p0, bounds=bounds, maxfev=30000)
        yhat = _sigmoid_model(x, *params)
        yhat_raw = yhat * y_range + y_min
        r, p_val = _safe_pearson(y, yhat)
        rmse, aic, aicc, bic = _model_metrics(y, yhat, n_params=4)
        
        # Quality checks
        is_valid, reasons, checks = _sigmoid_quality_checks(x, params, pcov)
        
        L, k, x0, b = params
        
        fit_results.append({
            'animal': animal,
            'n_trials': len(y),
            'L': L,
            'k': k,
            'x0': x0,
            'b': b,
            'r': r,
            'p': p_val,
            'rmse': rmse,
            'aicc': aicc,
            'bic': bic,
            'is_valid': is_valid,
            'reasons': reasons,
            'x0_interior': checks['x0_interior'],
            'k_plausible': checks['k_plausible'],
            'ci_finite': checks['ci_finite'],
            'asymptotes_covered': checks['asymptotes_covered'],
        })

        fit_traces.append({
            'animal': animal,
            'x': x,
            'y_raw': y_raw,
            'y_fit': yhat_raw,
            'is_valid': is_valid,
        })
        
    except Exception as e:
        print(f"  {animal}: FIT FAILED ({e})")
        fit_results.append({
            'animal': animal,
            'n_trials': len(y),
            'L': np.nan,
            'k': np.nan,
            'x0': np.nan,
            'b': np.nan,
            'r': np.nan,
            'p': np.nan,
            'rmse': np.nan,
            'aicc': np.nan,
            'bic': np.nan,
            'is_valid': False,
            'reasons': 'fit_failed',
            'x0_interior': False,
            'k_plausible': False,
            'ci_finite': False,
            'asymptotes_covered': False,
        })

# Create results DataFrame
sigmoid_fit_df = pd.DataFrame(fit_results)

# Round for display
display_df = sigmoid_fit_df.copy()
for col in ['L', 'k', 'x0', 'b', 'r', 'p', 'rmse', 'aicc', 'bic']:
    if col in display_df.columns:
        display_df[col] = display_df[col].round(4)

print("\n" + "="*80)
print("SIGMOID FIT RESULTS FOR time_moving (Deplete + 45NaCl)")
print("="*80)
print("\nFit Parameters:")
print(display_df[['animal', 'n_trials', 'L', 'k', 'x0', 'b', 'r', 'p', 'rmse', 'aicc']].to_string(index=False))

print("\n" + "-"*80)
print("Quality Checks:")
print(display_df[['animal', 'is_valid', 'reasons', 'x0_interior', 'k_plausible', 'ci_finite', 'asymptotes_covered']].to_string(index=False))

# Summary statistics
valid_fits = sigmoid_fit_df[sigmoid_fit_df['is_valid'] == True]
print("\n" + "="*80)
print(f"SUMMARY: {len(valid_fits)}/{len(sigmoid_fit_df)} fits passed all quality checks")
print("="*80)

if len(valid_fits) > 0:
    print(f"\nValid fits statistics:")
    print(f"  k (steepness):  mean = {valid_fits['k'].mean():.4f}, std = {valid_fits['k'].std():.4f}, range = [{valid_fits['k'].min():.4f}, {valid_fits['k'].max():.4f}]")
    print(f"  x0 (transition): mean = {valid_fits['x0'].mean():.2f}, std = {valid_fits['x0'].std():.2f}, range = [{valid_fits['x0'].min():.2f}, {valid_fits['x0'].max():.2f}]")
    print(f"  r (correlation): mean = {valid_fits['r'].mean():.4f}, std = {valid_fits['r'].std():.4f}")
    print(f"  RMSE:            mean = {valid_fits['rmse'].mean():.4f}, std = {valid_fits['rmse'].std():.4f}")

# Store results for potential use
time_moving_sigmoid_fits = sigmoid_fit_df


# %%
fits_df

# %%
## Plot per-animal sigmoid fits (time_moving)
if len(fit_traces) == 0:
    print("No fits to plot. Run the fitting cell first.")
else:
    for trace in fit_traces:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.scatter(trace['x'], trace['y_raw'], s=30, color='black', alpha=0.7, label='data')
        ax.plot(trace['x'], trace['y_fit'], color='tab:red', linewidth=2, label='sigmoid fit')
        ax.set_title(f"{trace['animal']} time_moving sigmoid fit")
        ax.set_xlabel('trial')
        ax.set_ylabel('time_moving')
        ax.legend(frameon=False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        fig.tight_layout()
        if savefigs:
            outpath = FIGSFOLDER / f"sigmoid_fit_time_moving_{trace['animal']}.png"
            fig.savefig(outpath, dpi=300)


# %%
## Re-fit sigmoids after binarizing auc_simba (above vs below zero)
# Compare binary fits against the original continuous auc_simba fits.
# For binarized data, very steep k values are acceptable and do not invalidate a fit.

binary_fit_results = []
binary_fit_traces = []

for animal in animals:
    animal_data = (
        deplete_45_subset
        .query("id == @animal")
        .sort_values('trial')
        .copy()
    )

    y_raw = animal_data.auc_simba.values
    y_binary = (y_raw > 0).astype(float)
    x = np.arange(1, len(y_binary) + 1, dtype=float)

    if np.unique(y_binary).size < 2:
        print(f"  {animal}: SKIPPED (binary auc_simba has no variation)")
        continue

    p0 = [
        1.0,
        -1.0,
        np.median(x),
        0.0,
    ]
    bounds = (
        [-np.inf, -np.inf, np.min(x), -np.inf],
        [np.inf, np.inf, np.max(x), np.inf],
    )

    try:
        params, pcov = curve_fit(_sigmoid_model, x, y_binary, p0=p0, bounds=bounds, maxfev=30000)
        yhat = _sigmoid_model(x, *params)
        r, p_val = _safe_pearson(y_binary, yhat)
        rmse, aic, aicc, bic = _model_metrics(y_binary, yhat, n_params=4)
        _, _, checks = _sigmoid_quality_checks(x, params, pcov)

        binary_checks = dict(checks)
        binary_checks['k_plausible'] = True
        failed = [name for name, ok in binary_checks.items() if not ok]
        is_valid = len(failed) == 0
        reasons = 'ok' if is_valid else ';'.join(failed)

        L, k, x0, b = params
        binary_fit_results.append({
            'animal': animal,
            'n_trials': len(y_binary),
            'n_positive_trials': int(y_binary.sum()),
            'L': L,
            'k': k,
            'x0': x0,
            'b': b,
            'r': r,
            'p': p_val,
            'rmse': rmse,
            'aicc': aicc,
            'bic': bic,
            'is_valid': is_valid,
            'reasons': reasons,
            'x0_interior': binary_checks['x0_interior'],
            'k_plausible': binary_checks['k_plausible'],
            'ci_finite': binary_checks['ci_finite'],
            'asymptotes_covered': binary_checks['asymptotes_covered'],
        })

        binary_fit_traces.append({
            'animal': animal,
            'x': x,
            'y_binary': y_binary,
            'y_fit': yhat,
            'is_valid': is_valid,
        })

    except Exception as e:
        print(f"  {animal}: FIT FAILED ({e})")
        binary_fit_results.append({
            'animal': animal,
            'n_trials': len(y_binary),
            'n_positive_trials': int(y_binary.sum()),
            'L': np.nan,
            'k': np.nan,
            'x0': np.nan,
            'b': np.nan,
            'r': np.nan,
            'p': np.nan,
            'rmse': np.nan,
            'aicc': np.nan,
            'bic': np.nan,
            'is_valid': False,
            'reasons': 'fit_failed',
            'x0_interior': False,
            'k_plausible': True,
            'ci_finite': False,
            'asymptotes_covered': False,
        })

auc_simba_binary_sigmoid_fits = pd.DataFrame(binary_fit_results)

comparison_df = (
    time_moving_sigmoid_fits[['animal', 'rmse', 'aicc', 'bic', 'is_valid']]
    .rename(columns={
        'rmse': 'rmse_continuous',
        'aicc': 'aicc_continuous',
        'bic': 'bic_continuous',
        'is_valid': 'is_valid_continuous',
    })
    .merge(
        auc_simba_binary_sigmoid_fits[['animal', 'n_positive_trials', 'rmse', 'aicc', 'bic', 'is_valid', 'x0', 'k']],
        on='animal',
        how='outer'
    )
    .rename(columns={
        'rmse': 'rmse_binary',
        'aicc': 'aicc_binary',
        'bic': 'bic_binary',
        'is_valid': 'is_valid_binary',
        'x0': 'x0_binary',
        'k': 'k_binary',
    })
)

comparison_df['rmse_delta'] = comparison_df['rmse_binary'] - comparison_df['rmse_continuous']
comparison_df['aicc_delta'] = comparison_df['aicc_binary'] - comparison_df['aicc_continuous']
comparison_df['bic_delta'] = comparison_df['bic_binary'] - comparison_df['bic_continuous']
comparison_df['binary_better_rmse'] = comparison_df['rmse_delta'] < 0
comparison_df['binary_better_aicc'] = comparison_df['aicc_delta'] < 0
comparison_df['binary_better_bic'] = comparison_df['bic_delta'] < 0

binary_display_df = auc_simba_binary_sigmoid_fits.copy()
for col in ['L', 'k', 'x0', 'b', 'r', 'p', 'rmse', 'aicc', 'bic']:
    if col in binary_display_df.columns:
        binary_display_df[col] = binary_display_df[col].round(4)

comparison_display_df = comparison_df.copy()
for col in ['rmse_continuous', 'rmse_binary', 'rmse_delta', 'aicc_continuous', 'aicc_binary', 'aicc_delta', 'bic_continuous', 'bic_binary', 'bic_delta', 'x0_binary', 'k_binary']:
    if col in comparison_display_df.columns:
        comparison_display_df[col] = comparison_display_df[col].round(4)

print("\n" + "=" * 80)
print("SIGMOID FIT RESULTS FOR BINARIZED auc_simba (> 0 vs <= 0)")
print("=" * 80)
print(binary_display_df[['animal', 'n_trials', 'n_positive_trials', 'L', 'k', 'x0', 'b', 'r', 'p', 'rmse', 'aicc', 'is_valid', 'reasons']].to_string(index=False))

print("\n" + "=" * 80)
print("COMPARISON TO ORIGINAL CONTINUOUS auc_simba FITS")
print("=" * 80)
print(comparison_display_df[['animal', 'n_positive_trials', 'is_valid_continuous', 'is_valid_binary', 'rmse_continuous', 'rmse_binary', 'rmse_delta', 'aicc_continuous', 'aicc_binary', 'aicc_delta', 'bic_continuous', 'bic_binary', 'bic_delta']].to_string(index=False))

print("\nSummary of binary-vs-continuous fit quality:")
print(f"  Better RMSE: {int(comparison_df['binary_better_rmse'].fillna(False).sum())}/{len(comparison_df)} animals")
print(f"  Better AICc: {int(comparison_df['binary_better_aicc'].fillna(False).sum())}/{len(comparison_df)} animals")
print(f"  Better BIC:  {int(comparison_df['binary_better_bic'].fillna(False).sum())}/{len(comparison_df)} animals")

comparison_df

# %%
## Plot binarized auc_simba sigmoid fits and compare fit metrics
if 'binary_fit_traces' not in globals() or len(binary_fit_traces) == 0:
    print("No binarized fit traces available. Run the previous cell first.")
else:
    n_animals = len(binary_fit_traces)
    n_cols = 3
    n_rows = int(np.ceil(n_animals / n_cols))
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 3.5 * n_rows), squeeze=False)
    axes_flat = axes.flatten()
    
    for ax, trace in zip(axes_flat, binary_fit_traces):
        animal = trace['animal']
        row = auc_simba_binary_sigmoid_fits.loc[auc_simba_binary_sigmoid_fits['animal'] == animal].iloc[0]
        ax.scatter(trace['x'], trace['y_binary'], s=28, color='black', alpha=0.75, label='binarized auc_simba')
        ax.plot(trace['x'], trace['y_fit'], color=colors[2], linewidth=2, label='sigmoid fit')
        if np.isfinite(row['x0']):
            ax.axvline(row['x0'], color=colors[1], linestyle='--', linewidth=1.5, alpha=0.8)
        ax.set_title(
            f"{animal} | valid={bool(row['is_valid'])} | x0={row['x0']:.1f}",
            fontsize=10
        )
        ax.set_xlabel('Trial')
        ax.set_ylabel('auc_simba > 0')
        ax.set_ylim(-0.1, 1.1)
        ax.set_yticks([0, 1])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    
    for ax in axes_flat[n_animals:]:
        ax.set_visible(False)
    
    handles, labels = axes_flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', ncol=2, frameon=False)
    fig.suptitle('Binarized auc_simba sigmoid fits by animal', y=1.02, fontsize=14)
    fig.tight_layout()
    plt.show()
    
    if 'comparison_df' in globals() and not comparison_df.empty:
        metrics = [
            ('rmse_continuous', 'rmse_binary', 'RMSE'),
            ('aicc_continuous', 'aicc_binary', 'AICc'),
            ('bic_continuous', 'bic_binary', 'BIC'),
        ]
        
        fig, axes = plt.subplots(1, len(metrics), figsize=(4.5 * len(metrics), 4), squeeze=False)
        axes = axes[0]
        
        for ax, (cont_col, bin_col, label) in zip(axes, metrics):
            plot_df = comparison_df[[cont_col, bin_col, 'animal']].dropna().copy()
            if plot_df.empty:
                ax.set_visible(False)
                continue
            
            lower = min(plot_df[cont_col].min(), plot_df[bin_col].min())
            upper = max(plot_df[cont_col].max(), plot_df[bin_col].max())
            ax.scatter(plot_df[cont_col], plot_df[bin_col], color=colors[0], s=45, alpha=0.85)
            ax.plot([lower, upper], [lower, upper], color='0.5', linestyle='--', linewidth=1)
            for _, row in plot_df.iterrows():
                ax.text(row[cont_col], row[bin_col], row['animal'], fontsize=8, alpha=0.8)
            ax.set_xlabel(f'Continuous {label}')
            ax.set_ylabel(f'Binary {label}')
            ax.set_title(f'Binary vs continuous {label}')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
        
        fig.tight_layout()
        plt.show()
    else:
        print("comparison_df is not available, so metric comparison plots were skipped.")

# %%
## Replot dopamine and behavior aligned to fitted binarized auc_simba x0 values
if 'auc_simba_binary_sigmoid_fits' not in globals():
    print("Run the binarized auc_simba sigmoid fitting cell first.")
else:
    alignment_fits = (
        auc_simba_binary_sigmoid_fits
        .dropna(subset=['x0'])
        .query('is_valid == True')
        .copy()
    )

    if alignment_fits.empty:
        print("No valid binarized fits with finite x0 values are available for realignment.")
    else:
        fit_lookup = (
            alignment_fits[['animal', 'x0']]
            .rename(columns={'animal': 'id'})
            .assign(x0_trial=lambda df: df['x0'].round().astype(int))
        )

        realigned_auc_simba = (
            deplete_45_subset
            .merge(fit_lookup, on='id', how='inner')
            .assign(trial_aligned=lambda df: df['trial'] - df['x0_trial'])
            .sort_values(['id', 'trial'])
            .copy()
        )

        n_required = realigned_auc_simba['id'].nunique()
        realigned_complete = (
            realigned_auc_simba
            .groupby('trial_aligned', group_keys=False)
            .filter(lambda group: len(group) == n_required)
            .copy()
        )

        def _fit_group_sigmoid(df, group_col, value_col):
            grouped = df.groupby(group_col)[value_col]
            mean = grouped.mean().sort_index()
            sem = grouped.sem().sort_index().fillna(0)
            x_vals = mean.index.to_numpy(dtype=float)
            y_vals = mean.to_numpy(dtype=float)

            popt = None
            y_fit = None
            if len(x_vals) >= 4 and np.nanmax(y_vals) - np.nanmin(y_vals) > 1e-8:
                p0 = [
                    float(np.nanmax(y_vals) - np.nanmin(y_vals)),
                    -1.0,
                    float(np.median(x_vals)),
                    float(np.nanmin(y_vals)),
                ]
                bounds = (
                    [-np.inf, -np.inf, np.min(x_vals), -np.inf],
                    [np.inf, np.inf, np.max(x_vals), np.inf],
                )
                try:
                    popt, _ = curve_fit(_sigmoid_model, x_vals, y_vals, p0=p0, bounds=bounds, maxfev=30000)
                    y_fit = _sigmoid_model(x_vals, *popt)
                except Exception as exc:
                    print(f"Failed to fit {value_col} by {group_col}: {exc}")

            return {
                'x': x_vals,
                'mean': y_vals,
                'sem': sem.reindex(mean.index).to_numpy(dtype=float),
                'popt': popt,
                'y_fit': y_fit,
            }

        dopamine_original = _fit_group_sigmoid(realigned_auc_simba, 'trial', 'auc_snips')
        dopamine_realigned = _fit_group_sigmoid(realigned_complete, 'trial_aligned', 'auc_snips')
        behaviour_original = _fit_group_sigmoid(realigned_auc_simba, 'trial', 'auc_simba')
        behaviour_realigned = _fit_group_sigmoid(realigned_complete, 'trial_aligned', 'auc_simba')

        fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.6), sharey=False)

        panel_specs = [
            (
                axes[0],
                dopamine_original,
                behaviour_original,
                'trial',
                'Original trial order',
                None,
            ),
            (
                axes[1],
                dopamine_realigned,
                behaviour_realigned,
                'trial_aligned',
                'Aligned to binarized-fit x0',
                0,
            ),
        ]

        behaviour_axes = []
        for ax, dopamine_data, behaviour_data, x_label, title, vline in panel_specs:
            ax.plot(dopamine_data['x'], dopamine_data['mean'], linestyle='', marker='o', markersize=4.5, color=colors[3], markerfacecolor='white')
            ax.fill_between(
                dopamine_data['x'],
                dopamine_data['mean'] - dopamine_data['sem'],
                dopamine_data['mean'] + dopamine_data['sem'],
                color=colors[3],
                alpha=0.12,
            )
            if dopamine_data['y_fit'] is not None:
                ax.plot(dopamine_data['x'], dopamine_data['y_fit'], color=colors[3], linestyle='--', linewidth=2)

            ax2 = ax.twinx()
            behaviour_axes.append(ax2)
            ax2.plot(behaviour_data['x'], behaviour_data['mean'], linestyle='', marker='o', markersize=4.5, color=colors[1], markerfacecolor='white')
            ax2.fill_between(
                behaviour_data['x'],
                behaviour_data['mean'] - behaviour_data['sem'],
                behaviour_data['mean'] + behaviour_data['sem'],
                color=colors[1],
                alpha=0.12,
            )
            if behaviour_data['y_fit'] is not None:
                ax2.plot(behaviour_data['x'], behaviour_data['y_fit'], color=colors[1], linestyle='--', linewidth=2)

            if vline is not None:
                ax.axvline(vline, color='0.5', linestyle='--', linewidth=1, alpha=0.8)

            ax.set_title(title)
            ax.set_xlabel('Trial' if x_label == 'trial' else 'Trial relative to x0')
            ax.tick_params(axis='y', labelcolor=colors[3])
            ax2.tick_params(axis='y', labelcolor=colors[1])
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax2.spines['top'].set_visible(False)

        axes[0].set_ylabel('Dopamine AUC', color=colors[3])
        axes[1].set_ylabel('Dopamine AUC', color=colors[3])
        behaviour_axes[0].set_ylabel('Behavior AUC', color=colors[1], rotation=270, labelpad=14)
        behaviour_axes[1].set_ylabel('Behavior AUC', color=colors[1], rotation=270, labelpad=14)

        for idx, (dopamine_data, behaviour_data) in enumerate([
            (dopamine_original, behaviour_original),
            (dopamine_realigned, behaviour_realigned),
        ]):
            text_lines = []
            if dopamine_data['popt'] is not None:
                text_lines.append(f"DA k={dopamine_data['popt'][1]:.2f}")
            if behaviour_data['popt'] is not None:
                text_lines.append(f"Beh k={behaviour_data['popt'][1]:.2f}")
            if text_lines:
                axes[idx].text(0.98, 0.95, '\n'.join(text_lines), transform=axes[idx].transAxes, ha='right', va='top', fontsize=9)

        fig.suptitle(f"Dopamine and behavior realigned to binarized auc_simba x0 ({n_required} animals)", y=1.03)
        fig.tight_layout()
        plt.show()

        print(f"Included animals: {sorted(realigned_auc_simba['id'].unique())}")
        print(f"Original trials: {len(realigned_auc_simba)}")
        print(f"Complete realigned trials: {len(realigned_complete)}")
        print(f"Realigned trial range: {int(realigned_complete['trial_aligned'].min())} to {int(realigned_complete['trial_aligned'].max())}")

# %% [markdown]
# # Grid search to find binarization threshold (no median balance)

# %%
## Grid search binarization threshold and optimize k-differences after realignment
# Goal: find threshold giving the largest k-change for BOTH dopamine and behavior.

if 'deplete_45_subset' not in globals():
    print("Run the sigmoid-fitting setup cell first (needs deplete_45_subset).")
elif '_sigmoid_model' not in globals() or '_sigmoid_quality_checks' not in globals():
    print("Run the sigmoid helper-definition cell first (needs _sigmoid_model and _sigmoid_quality_checks).")
else:
    threshold_values = np.linspace(-50, 50, 20)

    def _fit_binary_transition_table(threshold):
        rows = []
        for animal in sorted(deplete_45_subset['id'].unique()):
            animal_data = (
                deplete_45_subset
                .query("id == @animal")
                .sort_values('trial')
                .copy()
            )

            y_raw = animal_data.auc_simba.to_numpy(dtype=float)
            y_binary = (y_raw > threshold).astype(float)
            x_vals = np.arange(1, len(y_binary) + 1, dtype=float)

            if np.unique(y_binary).size < 2:
                rows.append({
                    'animal': animal,
                    'x0': np.nan,
                    'k': np.nan,
                    'is_valid': False,
                    'reasons': 'no_binary_variation',
                })
                continue

            p0 = [1.0, -1.0, np.median(x_vals), 0.0]
            bounds = (
                [-np.inf, -np.inf, np.min(x_vals), -np.inf],
                [np.inf, np.inf, np.max(x_vals), np.inf],
            )

            try:
                params, pcov = curve_fit(_sigmoid_model, x_vals, y_binary, p0=p0, bounds=bounds, maxfev=30000)
                _, _, checks = _sigmoid_quality_checks(x_vals, params, pcov)

                # For binarized data, steep k is acceptable.
                checks_binary = dict(checks)
                checks_binary['k_plausible'] = True
                failed = [name for name, ok in checks_binary.items() if not ok]
                is_valid = len(failed) == 0
                reasons = 'ok' if is_valid else ';'.join(failed)

                rows.append({
                    'animal': animal,
                    'x0': float(params[2]),
                    'k': float(params[1]),
                    'is_valid': bool(is_valid),
                    'reasons': reasons,
                })
            except Exception:
                rows.append({
                    'animal': animal,
                    'x0': np.nan,
                    'k': np.nan,
                    'is_valid': False,
                    'reasons': 'fit_failed',
                })

        return pd.DataFrame(rows)

    def _fit_group_k(df, group_col, value_col):
        grouped = df.groupby(group_col)[value_col].mean().sort_index()
        x_vals = grouped.index.to_numpy(dtype=float)
        y_vals = grouped.to_numpy(dtype=float)

        if len(x_vals) < 4 or (np.nanmax(y_vals) - np.nanmin(y_vals)) < 1e-8:
            return np.nan

        p0 = [
            float(np.nanmax(y_vals) - np.nanmin(y_vals)),
            -1.0,
            float(np.median(x_vals)),
            float(np.nanmin(y_vals)),
        ]
        bounds = (
            [-np.inf, -np.inf, np.min(x_vals), -np.inf],
            [np.inf, np.inf, np.max(x_vals), np.inf],
        )

        try:
            params, _ = curve_fit(_sigmoid_model, x_vals, y_vals, p0=p0, bounds=bounds, maxfev=30000)
            return float(params[1])
        except Exception:
            return np.nan

    grid_rows = []

    for threshold in threshold_values:
        threshold_fits = _fit_binary_transition_table(threshold)
        valid_fits = threshold_fits.query("is_valid == True").dropna(subset=['x0']).copy()

        if valid_fits.empty:
            grid_rows.append({
                'threshold': float(threshold),
                'n_animals': 0,
                'n_complete_trials': 0,
                'k_da_orig': np.nan,
                'k_da_realigned': np.nan,
                'k_beh_orig': np.nan,
                'k_beh_realigned': np.nan,
                'delta_k_da': np.nan,
                'delta_k_beh': np.nan,
                'score_min': np.nan,
                'score_sum': np.nan,
            })
            continue

        fit_lookup = (
            valid_fits[['animal', 'x0']]
            .rename(columns={'animal': 'id'})
            .assign(x0_trial=lambda df: df['x0'].round().astype(int))
        )

        aligned_df = (
            deplete_45_subset
            .merge(fit_lookup, on='id', how='inner')
            .assign(trial_aligned=lambda df: df['trial'] - df['x0_trial'])
            .sort_values(['id', 'trial'])
            .copy()
        )

        n_animals = int(aligned_df['id'].nunique())
        complete_df = (
            aligned_df
            .groupby('trial_aligned', group_keys=False)
            .filter(lambda group: len(group) == n_animals)
            .copy()
        )

        k_da_orig = _fit_group_k(aligned_df, 'trial', 'auc_snips')
        k_da_realigned = _fit_group_k(complete_df, 'trial_aligned', 'auc_snips') if not complete_df.empty else np.nan
        k_beh_orig = _fit_group_k(aligned_df, 'trial', 'auc_simba')
        k_beh_realigned = _fit_group_k(complete_df, 'trial_aligned', 'auc_simba') if not complete_df.empty else np.nan

        delta_k_da = abs(k_da_realigned - k_da_orig) if np.isfinite(k_da_realigned) and np.isfinite(k_da_orig) else np.nan
        delta_k_beh = abs(k_beh_realigned - k_beh_orig) if np.isfinite(k_beh_realigned) and np.isfinite(k_beh_orig) else np.nan
        score_min = min(delta_k_da, delta_k_beh) if np.isfinite(delta_k_da) and np.isfinite(delta_k_beh) else np.nan
        score_sum = (delta_k_da + delta_k_beh) if np.isfinite(delta_k_da) and np.isfinite(delta_k_beh) else np.nan

        grid_rows.append({
            'threshold': float(threshold),
            'n_animals': n_animals,
            'n_complete_trials': int(len(complete_df)),
            'k_da_orig': k_da_orig,
            'k_da_realigned': k_da_realigned,
            'k_beh_orig': k_beh_orig,
            'k_beh_realigned': k_beh_realigned,
            'delta_k_da': delta_k_da,
            'delta_k_beh': delta_k_beh,
            'score_min': score_min,
            'score_sum': score_sum,
        })

    threshold_grid_results = pd.DataFrame(grid_rows)

    valid_scores = threshold_grid_results.dropna(subset=['score_min', 'score_sum']).copy()
    if valid_scores.empty:
        print("No valid thresholds produced both dopamine and behavior k-differences.")
    else:
        best_threshold_row = (
            valid_scores
            .sort_values(['score_min', 'score_sum'], ascending=False)
            .iloc[0]
        )

        print("Best threshold (maximize BOTH k differences):")
        print(
            f"  threshold={best_threshold_row['threshold']:.2f}, "
            f"Δk_dopamine={best_threshold_row['delta_k_da']:.4f}, "
            f"Δk_behavior={best_threshold_row['delta_k_beh']:.4f}, "
            f"score_min={best_threshold_row['score_min']:.4f}, "
            f"score_sum={best_threshold_row['score_sum']:.4f}, "
            f"animals={int(best_threshold_row['n_animals'])}"
        )

        display_cols = [
            'threshold', 'n_animals', 'n_complete_trials',
            'delta_k_da', 'delta_k_beh', 'score_min', 'score_sum',
            'k_da_orig', 'k_da_realigned', 'k_beh_orig', 'k_beh_realigned'
        ]
        display_df = (
            valid_scores[display_cols]
            .sort_values(['score_min', 'score_sum'], ascending=False)
            .reset_index(drop=True)
        )

        print("\nTop 10 thresholds:")
        print(display_df.head(10).round(4).to_string(index=False))

        fig, axes = plt.subplots(1, 2, figsize=(11.5, 4), squeeze=False)
        ax0, ax1 = axes[0]

        plot_df = threshold_grid_results.sort_values('threshold')
        ax0.plot(plot_df['threshold'], plot_df['delta_k_da'], marker='o', color=colors[3], label='Δk dopamine')
        ax0.plot(plot_df['threshold'], plot_df['delta_k_beh'], marker='o', color=colors[1], label='Δk behavior')
        ax0.axvline(best_threshold_row['threshold'], linestyle='--', color='0.4', linewidth=1)
        ax0.set_xlabel('Binarization threshold')
        ax0.set_ylabel('|k_realigned - k_original|')
        ax0.set_title('k-difference vs threshold')
        ax0.legend(frameon=False)
        ax0.spines['top'].set_visible(False)
        ax0.spines['right'].set_visible(False)

        ax1.plot(plot_df['threshold'], plot_df['score_min'], marker='o', color=colors[0], label='score_min = min(Δk_DA, Δk_Beh)')
        ax1.plot(plot_df['threshold'], plot_df['score_sum'], marker='o', color=colors[2], alpha=0.7, label='score_sum = Δk_DA + Δk_Beh')
        ax1.axvline(best_threshold_row['threshold'], linestyle='--', color='0.4', linewidth=1)
        ax1.set_xlabel('Binarization threshold')
        ax1.set_ylabel('Score')
        ax1.set_title('Grid-search objective scores')
        ax1.legend(frameon=False, fontsize=9)
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)

        fig.tight_layout()
        plt.show()

    threshold_grid_results

# %%
## Try custom binarization thresholds and inspect realignment effects
# Edit this list as needed
manual_thresholds = [-7.89]

if 'deplete_45_subset' not in globals():
    print("Run the sigmoid-fitting setup cell first (needs deplete_45_subset).")
elif '_sigmoid_model' not in globals() or '_sigmoid_quality_checks' not in globals():
    print("Run the sigmoid helper-definition cell first (needs _sigmoid_model and _sigmoid_quality_checks).")
else:
    def _fit_binary_transition_table_single_threshold(threshold):
        rows = []
        for animal in sorted(deplete_45_subset['id'].unique()):
            animal_data = (
                deplete_45_subset
                .query("id == @animal")
                .sort_values('trial')
                .copy()
            )

            y_raw = animal_data.auc_simba.to_numpy(dtype=float)
            y_binary = (y_raw > threshold).astype(float)
            x_vals = np.arange(1, len(y_binary) + 1, dtype=float)

            if np.unique(y_binary).size < 2:
                rows.append({'animal': animal, 'x0': np.nan, 'k': np.nan, 'is_valid': False, 'reasons': 'no_binary_variation'})
                continue

            p0 = [1.0, -1.0, np.median(x_vals), 0.0]
            bounds = (
                [-np.inf, -np.inf, np.min(x_vals), -np.inf],
                [np.inf, np.inf, np.max(x_vals), np.inf],
            )

            try:
                params, pcov = curve_fit(_sigmoid_model, x_vals, y_binary, p0=p0, bounds=bounds, maxfev=30000)
                _, _, checks = _sigmoid_quality_checks(x_vals, params, pcov)
                checks_binary = dict(checks)
                checks_binary['k_plausible'] = True
                failed = [name for name, ok in checks_binary.items() if not ok]
                is_valid = len(failed) == 0
                reasons = 'ok' if is_valid else ';'.join(failed)
                rows.append({'animal': animal, 'x0': float(params[2]), 'k': float(params[1]), 'is_valid': bool(is_valid), 'reasons': reasons})
            except Exception:
                rows.append({'animal': animal, 'x0': np.nan, 'k': np.nan, 'is_valid': False, 'reasons': 'fit_failed'})

        return pd.DataFrame(rows)

    def _fit_group_k(df, group_col, value_col):
        grouped = df.groupby(group_col)[value_col].mean().sort_index()
        x_vals = grouped.index.to_numpy(dtype=float)
        y_vals = grouped.to_numpy(dtype=float)

        if len(x_vals) < 4 or (np.nanmax(y_vals) - np.nanmin(y_vals)) < 1e-8:
            return np.nan

        p0 = [
            float(np.nanmax(y_vals) - np.nanmin(y_vals)),
            -1.0,
            float(np.median(x_vals)),
            float(np.nanmin(y_vals)),
        ]
        bounds = (
            [-np.inf, -np.inf, np.min(x_vals), -np.inf],
            [np.inf, np.inf, np.max(x_vals), np.inf],
        )

        try:
            params, _ = curve_fit(_sigmoid_model, x_vals, y_vals, p0=p0, bounds=bounds, maxfev=30000)
            return float(params[1])
        except Exception:
            return np.nan

    manual_rows = []

    for threshold in manual_thresholds:
        threshold_fits = _fit_binary_transition_table_single_threshold(float(threshold))
        valid_fits = threshold_fits.query("is_valid == True").dropna(subset=['x0']).copy()

        if valid_fits.empty:
            print(f"Threshold {threshold:.2f}: no valid fits")
            manual_rows.append({
                'threshold': float(threshold),
                'n_animals': 0,
                'n_complete_trials': 0,
                'k_da_orig': np.nan,
                'k_da_realigned': np.nan,
                'k_beh_orig': np.nan,
                'k_beh_realigned': np.nan,
                'delta_k_da': np.nan,
                'delta_k_beh': np.nan,
                'score_min': np.nan,
            })
            continue

        fit_lookup = (
            valid_fits[['animal', 'x0']]
            .rename(columns={'animal': 'id'})
            .assign(x0_trial=lambda df: df['x0'].round().astype(int))
        )

        aligned_df = (
            deplete_45_subset
            .merge(fit_lookup, on='id', how='inner')
            .assign(trial_aligned=lambda df: df['trial'] - df['x0_trial'])
            .sort_values(['id', 'trial'])
            .copy()
        )

        n_animals = int(aligned_df['id'].nunique())
        complete_df = (
            aligned_df
            .groupby('trial_aligned', group_keys=False)
            .filter(lambda group: len(group) == n_animals)
            .copy()
        )

        k_da_orig = _fit_group_k(aligned_df, 'trial', 'auc_snips')
        k_da_realigned = _fit_group_k(complete_df, 'trial_aligned', 'auc_snips') if not complete_df.empty else np.nan
        k_beh_orig = _fit_group_k(aligned_df, 'trial', 'auc_simba')
        k_beh_realigned = _fit_group_k(complete_df, 'trial_aligned', 'auc_simba') if not complete_df.empty else np.nan

        delta_k_da = abs(k_da_realigned - k_da_orig) if np.isfinite(k_da_realigned) and np.isfinite(k_da_orig) else np.nan
        delta_k_beh = abs(k_beh_realigned - k_beh_orig) if np.isfinite(k_beh_realigned) and np.isfinite(k_beh_orig) else np.nan
        score_min = min(delta_k_da, delta_k_beh) if np.isfinite(delta_k_da) and np.isfinite(delta_k_beh) else np.nan

        manual_rows.append({
            'threshold': float(threshold),
            'n_animals': n_animals,
            'n_complete_trials': int(len(complete_df)),
            'k_da_orig': k_da_orig,
            'k_da_realigned': k_da_realigned,
            'k_beh_orig': k_beh_orig,
            'k_beh_realigned': k_beh_realigned,
            'delta_k_da': delta_k_da,
            'delta_k_beh': delta_k_beh,
            'score_min': score_min,
        })

        def _group_mean_sem(df, group_col, value_col):
            grouped = df.groupby(group_col)[value_col]
            mean = grouped.mean().sort_index()
            sem = grouped.sem().sort_index().fillna(0)
            return mean.index.to_numpy(dtype=float), mean.to_numpy(dtype=float), sem.to_numpy(dtype=float)

        x_do, y_do, sem_do = _group_mean_sem(aligned_df, 'trial', 'auc_snips')
        x_dr, y_dr, sem_dr = _group_mean_sem(complete_df, 'trial_aligned', 'auc_snips')
        x_bo, y_bo, sem_bo = _group_mean_sem(aligned_df, 'trial', 'auc_simba')
        x_br, y_br, sem_br = _group_mean_sem(complete_df, 'trial_aligned', 'auc_simba')

        fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.6), sharey=False)

        axes[0].plot(x_do, y_do, linestyle='', marker='o', markersize=4.5, color=colors[3], markerfacecolor='white')
        axes[0].fill_between(x_do, y_do - sem_do, y_do + sem_do, color=colors[3], alpha=0.12)
        if np.isfinite(k_da_orig):
            p0_do = [float(np.nanmax(y_do) - np.nanmin(y_do)), -1.0, float(np.median(x_do)), float(np.nanmin(y_do))]
            bounds_do = ([-np.inf, -np.inf, np.min(x_do), -np.inf], [np.inf, np.inf, np.max(x_do), np.inf])
            try:
                popt_do, _ = curve_fit(_sigmoid_model, x_do, y_do, p0=p0_do, bounds=bounds_do, maxfev=30000)
                axes[0].plot(x_do, _sigmoid_model(x_do, *popt_do), color=colors[3], linestyle='--', linewidth=2)
            except Exception:
                pass

        ax0b = axes[0].twinx()
        ax0b.plot(x_bo, y_bo, linestyle='', marker='o', markersize=4.5, color=colors[1], markerfacecolor='white')
        ax0b.fill_between(x_bo, y_bo - sem_bo, y_bo + sem_bo, color=colors[1], alpha=0.12)
        if np.isfinite(k_beh_orig):
            p0_bo = [float(np.nanmax(y_bo) - np.nanmin(y_bo)), -1.0, float(np.median(x_bo)), float(np.nanmin(y_bo))]
            bounds_bo = ([-np.inf, -np.inf, np.min(x_bo), -np.inf], [np.inf, np.inf, np.max(x_bo), np.inf])
            try:
                popt_bo, _ = curve_fit(_sigmoid_model, x_bo, y_bo, p0=p0_bo, bounds=bounds_bo, maxfev=30000)
                ax0b.plot(x_bo, _sigmoid_model(x_bo, *popt_bo), color=colors[1], linestyle='--', linewidth=2)
            except Exception:
                pass

        axes[1].plot(x_dr, y_dr, linestyle='', marker='o', markersize=4.5, color=colors[3], markerfacecolor='white')
        axes[1].fill_between(x_dr, y_dr - sem_dr, y_dr + sem_dr, color=colors[3], alpha=0.12)
        if np.isfinite(k_da_realigned):
            p0_dr = [float(np.nanmax(y_dr) - np.nanmin(y_dr)), -1.0, float(np.median(x_dr)), float(np.nanmin(y_dr))]
            bounds_dr = ([-np.inf, -np.inf, np.min(x_dr), -np.inf], [np.inf, np.inf, np.max(x_dr), np.inf])
            try:
                popt_dr, _ = curve_fit(_sigmoid_model, x_dr, y_dr, p0=p0_dr, bounds=bounds_dr, maxfev=30000)
                axes[1].plot(x_dr, _sigmoid_model(x_dr, *popt_dr), color=colors[3], linestyle='--', linewidth=2)
            except Exception:
                pass

        ax1b = axes[1].twinx()
        ax1b.plot(x_br, y_br, linestyle='', marker='o', markersize=4.5, color=colors[1], markerfacecolor='white')
        ax1b.fill_between(x_br, y_br - sem_br, y_br + sem_br, color=colors[1], alpha=0.12)
        if np.isfinite(k_beh_realigned):
            p0_br = [float(np.nanmax(y_br) - np.nanmin(y_br)), -1.0, float(np.median(x_br)), float(np.nanmin(y_br))]
            bounds_br = ([-np.inf, -np.inf, np.min(x_br), -np.inf], [np.inf, np.inf, np.max(x_br), np.inf])
            try:
                popt_br, _ = curve_fit(_sigmoid_model, x_br, y_br, p0=p0_br, bounds=bounds_br, maxfev=30000)
                ax1b.plot(x_br, _sigmoid_model(x_br, *popt_br), color=colors[1], linestyle='--', linewidth=2)
            except Exception:
                pass

        axes[1].axvline(0, color='0.5', linestyle='--', linewidth=1)

        axes[0].set_title('Original trial order')
        axes[1].set_title('Aligned to threshold-specific x0')
        axes[0].set_xlabel('Trial')
        axes[1].set_xlabel('Trial relative to x0')
        axes[0].set_ylabel('Dopamine AUC', color=colors[3])
        axes[1].set_ylabel('Dopamine AUC', color=colors[3])
        ax0b.set_ylabel('Behavior AUC', color=colors[1], rotation=270, labelpad=14)
        ax1b.set_ylabel('Behavior AUC', color=colors[1], rotation=270, labelpad=14)

        axes[0].tick_params(axis='y', labelcolor=colors[3])
        axes[1].tick_params(axis='y', labelcolor=colors[3])
        ax0b.tick_params(axis='y', labelcolor=colors[1])
        ax1b.tick_params(axis='y', labelcolor=colors[1])

        axes[0].text(0.98, 0.95, f"DA k={k_da_orig:.2f}\nBeh k={k_beh_orig:.2f}", transform=axes[0].transAxes, ha='right', va='top', fontsize=9)
        axes[1].text(0.98, 0.95, f"DA k={k_da_realigned:.2f}\nBeh k={k_beh_realigned:.2f}", transform=axes[1].transAxes, ha='right', va='top', fontsize=9)

        for ax in axes:
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

        fig.suptitle(f"Threshold = {threshold:.2f} | animals={n_animals} | complete trials={len(complete_df)}", y=1.03)
        fig.tight_layout()
        plt.show()

        print(
            f"threshold={threshold:.2f} | animals={n_animals} | complete_trials={len(complete_df)} | "
            f"Δk_dopamine={delta_k_da:.4f} | Δk_behavior={delta_k_beh:.4f} | score_min={score_min:.4f}"
        )

    manual_threshold_effects = pd.DataFrame(manual_rows)
    print("\nSummary across tested thresholds:")
    print(manual_threshold_effects.sort_values('score_min', ascending=False).round(4).to_string(index=False))
    manual_threshold_effects

# %%
## Third plot: normalized sigmoid fits on realigned trials -5 to +5
# Uses manual_thresholds from the previous cell.

if 'manual_thresholds' not in globals():
    print("Run the custom-threshold cell first (needs manual_thresholds).")
elif '_fit_binary_transition_table_single_threshold' not in globals():
    print("Run the custom-threshold cell first (needs helper fit functions).")
else:
    window_min, window_max = -5, 5

    def _normalize_unit(y_vals):
        y_vals = np.asarray(y_vals, dtype=float)
        y_min = np.nanmin(y_vals)
        y_max = np.nanmax(y_vals)
        y_range = y_max - y_min
        if not np.isfinite(y_range) or y_range < 1e-8:
            return np.full_like(y_vals, np.nan)
        return (y_vals - y_min) / y_range

    def _fit_norm_sigmoid(x_vals, y_vals):
        x_vals = np.asarray(x_vals, dtype=float)
        y_vals = np.asarray(y_vals, dtype=float)
        ok = np.isfinite(x_vals) & np.isfinite(y_vals)
        x_vals = x_vals[ok]
        y_vals = y_vals[ok]
        if len(x_vals) < 4:
            return None, None

        p0 = [1.0, -1.0, float(np.median(x_vals)), 0.0]
        bounds = (
            [-np.inf, -np.inf, np.min(x_vals), -np.inf],
            [np.inf, np.inf, np.max(x_vals), np.inf],
        )
        try:
            params, _ = curve_fit(_sigmoid_model, x_vals, y_vals, p0=p0, bounds=bounds, maxfev=30000)
            yhat = _sigmoid_model(x_vals, *params)
            return params, yhat
        except Exception:
            return None, None

    n_plots = len(manual_thresholds)
    fig, axes = plt.subplots(1, n_plots, figsize=(5 * n_plots, 3.8), squeeze=False)
    axes = axes[0]

    for ax, threshold in zip(axes, manual_thresholds):
        threshold_fits = _fit_binary_transition_table_single_threshold(float(threshold))
        valid_fits = threshold_fits.query("is_valid == True").dropna(subset=['x0']).copy()

        if valid_fits.empty:
            ax.set_title(f"Threshold {threshold:.2f}\n(no valid fits)")
            ax.axis('off')
            continue

        fit_lookup = (
            valid_fits[['animal', 'x0']]
            .rename(columns={'animal': 'id'})
            .assign(x0_trial=lambda df: df['x0'].round().astype(int))
        )

        aligned_df = (
            deplete_45_subset
            .merge(fit_lookup, on='id', how='inner')
            .assign(trial_aligned=lambda df: df['trial'] - df['x0_trial'])
            .sort_values(['id', 'trial'])
            .copy()
        )

        n_animals = int(aligned_df['id'].nunique())
        complete_df = (
            aligned_df
            .groupby('trial_aligned', group_keys=False)
            .filter(lambda group: len(group) == n_animals)
            .copy()
        )

        window_df = complete_df.query("@window_min <= trial_aligned <= @window_max").copy()
        if window_df.empty:
            ax.set_title(f"Threshold {threshold:.2f}\n(no data in -5..+5)")
            ax.axis('off')
            continue

        da_mean = window_df.groupby('trial_aligned')['auc_snips'].mean().sort_index()
        beh_mean = window_df.groupby('trial_aligned')['auc_simba'].mean().sort_index()

        x_da = da_mean.index.to_numpy(dtype=float)
        y_da_norm = _normalize_unit(da_mean.to_numpy(dtype=float))
        x_beh = beh_mean.index.to_numpy(dtype=float)
        y_beh_norm = _normalize_unit(beh_mean.to_numpy(dtype=float))

        da_params, da_fit = _fit_norm_sigmoid(x_da, y_da_norm)
        beh_params, beh_fit = _fit_norm_sigmoid(x_beh, y_beh_norm)

        ax.plot(x_da, y_da_norm, linestyle='', marker='o', markersize=4.5, color=colors[3], markerfacecolor='white', label='Dopamine (norm)')
        ax.plot(x_beh, y_beh_norm, linestyle='', marker='o', markersize=4.5, color=colors[1], markerfacecolor='white', label='Behavior (norm)')

        if da_fit is not None:
            ax.plot(x_da, da_fit, color=colors[3], linestyle='--', linewidth=2)
        if beh_fit is not None:
            ax.plot(x_beh, beh_fit, color=colors[1], linestyle='--', linewidth=2)

        ax.axvline(0, color='0.5', linestyle='--', linewidth=1)
        ax.set_xlim(window_min - 0.2, window_max + 0.2)
        ax.set_ylim(-0.05, 1.05)
        ax.set_xlabel('Trial relative to x0')
        ax.set_ylabel('Normalized signal (0-1)')
        ax.set_title(f"Threshold {threshold:.2f} | realigned -5 to +5")
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        text_lines = []
        if da_params is not None:
            text_lines.append(f"DA k={da_params[1]:.2f}")
        if beh_params is not None:
            text_lines.append(f"Beh k={beh_params[1]:.2f}")
        if text_lines:
            ax.text(0.98, 0.95, '\n'.join(text_lines), transform=ax.transAxes, ha='right', va='top', fontsize=9)

    axes[0].legend(frameon=False, loc='lower right')
    fig.suptitle('Normalized sigmoid fits after realignment (-5 to +5 trials)', y=1.03)
    fig.tight_layout()
    plt.show()

# %%

# %% [markdown]
# ## Sigmoidal Fits for `simba_median_balance`
#
# Fit per-animal sigmoids to `simba_median_balance` (deplete + 45NaCl) using the same quality checks used above.

# %%
## Fit sigmoids to simba_median_balance (Deplete + 45NaCl)
if 'simba_median_balance' not in x_array.columns:
    print("Column 'simba_median_balance' not found in x_array.")
    print("Available simba-like columns:")
    simba_cols = [c for c in x_array.columns if 'simba' in str(c).lower()]
    print(simba_cols if simba_cols else "  none found")
else:
    simba_subset = (
        x_array
        .query("condition == 'deplete' & infusiontype == '45NaCl'")
        .copy()
    )

    simba_animals = sorted(simba_subset['id'].unique())
    print(f"Fitting sigmoids for {len(simba_animals)} animals using simba_median_balance")

    simba_balance_fit_results = []
    simba_balance_fit_traces = []

    for animal in simba_animals:
        animal_data = (
            simba_subset
            .query("id == @animal")
            .sort_values('trial')
            .copy()
        )

        y_raw = animal_data['simba_median_balance'].to_numpy(dtype=float)
        finite_mask = np.isfinite(y_raw)
        x = np.arange(1, len(y_raw) + 1, dtype=float)

        if finite_mask.sum() < 4:
            print(f"  {animal}: SKIPPED (fewer than 4 finite samples)")
            continue

        x_fit = x[finite_mask]
        y_fit_raw = y_raw[finite_mask]

        y_min, y_max = np.nanmin(y_fit_raw), np.nanmax(y_fit_raw)
        y_range = y_max - y_min
        if y_range < 1e-8:
            print(f"  {animal}: SKIPPED (no variation in simba_median_balance)")
            continue

        y_norm = (y_fit_raw - y_min) / y_range

        p0 = [
            1.0,
            -1.0,
            float(np.median(x_fit)),
            0.0,
        ]
        bounds = (
            [-np.inf, -2, np.min(x_fit), -np.inf],
            [np.inf, 2, np.max(x_fit), np.inf],
        )

        try:
            params, pcov = curve_fit(_sigmoid_model, x_fit, y_norm, p0=p0, bounds=bounds, maxfev=30000)
            yhat_norm = _sigmoid_model(x_fit, *params)
            yhat_raw = yhat_norm * y_range + y_min

            r, p_val = _safe_pearson(y_norm, yhat_norm)
            rmse, aic, aicc, bic = _model_metrics(y_norm, yhat_norm, n_params=4)
            _, _, base_checks = _sigmoid_quality_checks(x_fit, params, pcov)

            L, k, x0, b = params

            # For this analysis, steep slopes are allowed; only positive k should fail.
            checks = dict(base_checks)
            checks['k_plausible'] = True
            checks['k_negative'] = bool(np.isfinite(k) and (k < 0))

            failed = [name for name, ok in checks.items() if not ok]
            is_valid = len(failed) == 0
            reasons = 'ok' if is_valid else ';'.join(failed)

            simba_balance_fit_results.append({
                'animal': animal,
                'n_trials': int(len(y_raw)),
                'n_finite': int(finite_mask.sum()),
                'L': float(L),
                'k': float(k),
                'x0': float(x0),
                'b': float(b),
                'r': float(r) if np.isfinite(r) else np.nan,
                'p': float(p_val) if np.isfinite(p_val) else np.nan,
                'rmse': float(rmse),
                'aicc': float(aicc),
                'bic': float(bic),
                'is_valid': bool(is_valid),
                'reasons': reasons,
                'x0_interior': checks['x0_interior'],
                'k_plausible': checks['k_plausible'],
                'k_negative': checks['k_negative'],
                'ci_finite': checks['ci_finite'],
                'asymptotes_covered': checks['asymptotes_covered'],
            })

            simba_balance_fit_traces.append({
                'animal': animal,
                'x': x_fit,
                'y_raw': y_fit_raw,
                'y_fit': yhat_raw,
                'is_valid': bool(is_valid),
            })

        except Exception as e:
            print(f"  {animal}: FIT FAILED ({e})")
            simba_balance_fit_results.append({
                'animal': animal,
                'n_trials': int(len(y_raw)),
                'n_finite': int(finite_mask.sum()),
                'L': np.nan,
                'k': np.nan,
                'x0': np.nan,
                'b': np.nan,
                'r': np.nan,
                'p': np.nan,
                'rmse': np.nan,
                'aicc': np.nan,
                'bic': np.nan,
                'is_valid': False,
                'reasons': 'fit_failed',
                'x0_interior': False,
                'k_plausible': True,
                'k_negative': False,
                'ci_finite': False,
                'asymptotes_covered': False,
            })

    simba_median_balance_sigmoid_fits = pd.DataFrame(simba_balance_fit_results)

    if simba_median_balance_sigmoid_fits.empty:
        print("No simba_median_balance fits were produced.")
    else:
        display_df = simba_median_balance_sigmoid_fits.copy()
        for col in ['L', 'k', 'x0', 'b', 'r', 'p', 'rmse', 'aicc', 'bic']:
            if col in display_df.columns:
                display_df[col] = display_df[col].round(4)

        print("\n" + "=" * 80)
        print("SIGMOID FIT RESULTS FOR simba_median_balance (Deplete + 45NaCl)")
        print("=" * 80)
        print(display_df[['animal', 'n_trials', 'n_finite', 'L', 'k', 'x0', 'b', 'r', 'p', 'rmse', 'aicc', 'is_valid', 'reasons']].to_string(index=False))

        valid_fits = simba_median_balance_sigmoid_fits[simba_median_balance_sigmoid_fits['is_valid'] == True]
        print("\n" + "=" * 80)
        print(f"SUMMARY: {len(valid_fits)}/{len(simba_median_balance_sigmoid_fits)} fits passed all quality checks")
        print("=" * 80)

        if len(valid_fits) > 0:
            print(f"\nValid fits statistics:")
            print(f"  k (steepness):   mean = {valid_fits['k'].mean():.4f}, std = {valid_fits['k'].std():.4f}")
            print(f"  x0 (transition): mean = {valid_fits['x0'].mean():.2f}, std = {valid_fits['x0'].std():.2f}")
            print(f"  RMSE:            mean = {valid_fits['rmse'].mean():.4f}, std = {valid_fits['rmse'].std():.4f}")

# %%
## Plot simba_median_balance fits for all animals (including failed fits)
if 'simba_median_balance_sigmoid_fits' not in globals() or simba_median_balance_sigmoid_fits.empty:
    print("Run the simba_median_balance fitting cell first.")
elif 'x_array' not in globals():
    print("x_array not found. Run the data loading cell first.")
else:
    plot_subset = (
        x_array
        .query("condition == 'deplete' & infusiontype == '45NaCl'")
        .copy()
    )

    if 'simba_median_balance' not in plot_subset.columns:
        print("Column 'simba_median_balance' not found in x_array.")
    else:
        plot_animals = sorted(plot_subset['id'].unique())
        fit_table = simba_median_balance_sigmoid_fits.copy()
        fit_lookup = fit_table.set_index('animal') if not fit_table.empty else pd.DataFrame()

        n_animals = len(plot_animals)
        n_cols = 3
        n_rows = int(np.ceil(n_animals / n_cols))

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(5.2 * n_cols, 3.8 * n_rows), squeeze=False)
        axes_flat = axes.flatten()

        summary_counts = {
            'valid': 0,
            'quality_failed': 0,
            'fit_failed_or_missing': 0,
        }

        for ax, animal in zip(axes_flat, plot_animals):
            animal_df = plot_subset.loc[plot_subset['id'] == animal].sort_values('trial').copy()
            y_raw = animal_df['simba_median_balance'].to_numpy(dtype=float)
            x_all = np.arange(1, len(y_raw) + 1, dtype=float)
            finite_mask = np.isfinite(y_raw)

            ax.scatter(x_all[finite_mask], y_raw[finite_mask], s=26, color='black', alpha=0.75, label='data')
            if (~finite_mask).any():
                ax.scatter(x_all[~finite_mask], np.zeros((~finite_mask).sum()), s=26, marker='x', color='0.6', alpha=0.8, label='NaN')

            status = 'fit missing'
            fit_color = '#999999'

            if animal in fit_lookup.index:
                row = fit_lookup.loc[animal]
                has_params = np.all(np.isfinite([row.get('L', np.nan), row.get('k', np.nan), row.get('x0', np.nan), row.get('b', np.nan)]))

                if has_params and finite_mask.sum() >= 4:
                    x_fit = x_all[finite_mask]
                    y_fit = _sigmoid_model(x_fit, row['L'], row['k'], row['x0'], row['b'])
                    is_valid = bool(row.get('is_valid', False))

                    if is_valid:
                        status = 'valid'
                        fit_color = '#2E8B57'
                        summary_counts['valid'] += 1
                    else:
                        status = f"invalid: {row.get('reasons', 'quality_failed')}"
                        fit_color = '#C74632'
                        summary_counts['quality_failed'] += 1

                    ax.plot(x_fit, y_fit, color=fit_color, linewidth=2.0, linestyle='--', label='fit')

                    if np.isfinite(row.get('x0', np.nan)):
                        ax.axvline(row['x0'], color=fit_color, linestyle=':', linewidth=1.2, alpha=0.9)
                else:
                    status = f"fit failed: {row.get('reasons', 'fit_failed')}"
                    summary_counts['fit_failed_or_missing'] += 1
            else:
                summary_counts['fit_failed_or_missing'] += 1

            ax.set_title(f"{animal} | {status}", fontsize=10)
            ax.set_xlabel('Trial')
            ax.set_ylabel('simba_median_balance')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

        for ax in axes_flat[n_animals:]:
            ax.set_visible(False)

        from matplotlib.lines import Line2D
        legend_handles = [
            Line2D([0], [0], marker='o', color='none', markerfacecolor='black', markersize=6, label='Data'),
            Line2D([0], [0], color='#2E8B57', linestyle='--', linewidth=2, label='Valid fit'),
            Line2D([0], [0], color='#C74632', linestyle='--', linewidth=2, label='Invalid fit (quality check failed)'),
            Line2D([0], [0], color='#999999', linestyle='--', linewidth=2, label='Fit failed/missing'),
        ]
        fig.legend(handles=legend_handles, loc='upper center', ncol=2, frameon=False)

        fig.suptitle(
            "simba_median_balance sigmoid fits by animal (deplete + 45NaCl)\n"
            f"valid={summary_counts['valid']} | quality_failed={summary_counts['quality_failed']} | fit_failed_or_missing={summary_counts['fit_failed_or_missing']}",
            y=1.03,
            fontsize=13,
        )
        fig.tight_layout()
        plt.show()

# %%
## Plot simba_median_balance fits for all animals (including failed fits)
if 'simba_median_balance_sigmoid_fits' not in globals() or simba_median_balance_sigmoid_fits.empty:
    print("Run the simba_median_balance fitting cell first.")
elif 'x_array' not in globals():
    print("x_array not found. Run the data loading cell first.")
else:
    plot_subset = (
        x_array
        .query("condition == 'deplete' & infusiontype == '45NaCl'")
        .copy()
    )

    if 'simba_median_balance' not in plot_subset.columns:
        print("Column 'simba_median_balance' not found in x_array.")
    else:
        plot_animals = sorted(plot_subset['id'].unique())
        fit_table = simba_median_balance_sigmoid_fits.copy()
        fit_lookup = fit_table.set_index('animal') if not fit_table.empty else pd.DataFrame()

        n_animals = len(plot_animals)
        n_cols = 3
        n_rows = int(np.ceil(n_animals / n_cols))

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(5.2 * n_cols, 3.8 * n_rows), squeeze=False)
        axes_flat = axes.flatten()

        summary_counts = {
            'valid': 0,
            'quality_failed': 0,
            'fit_failed_or_missing': 0,
        }

        for ax, animal in zip(axes_flat, plot_animals):
            animal_df = plot_subset.loc[plot_subset['id'] == animal].sort_values('trial').copy()
            y_raw = animal_df['simba_median_balance'].to_numpy(dtype=float)
            x_all = np.arange(1, len(y_raw) + 1, dtype=float)
            finite_mask = np.isfinite(y_raw)

            ax.scatter(x_all[finite_mask], y_raw[finite_mask], s=26, color='black', alpha=0.75)

            status = 'fit missing'
            fit_color = '#999999'

            if animal in fit_lookup.index:
                row = fit_lookup.loc[animal]
                has_params = np.all(np.isfinite([row.get('L', np.nan), row.get('k', np.nan), row.get('x0', np.nan), row.get('b', np.nan)]))

                if has_params and finite_mask.sum() >= 4:
                    x_fit = x_all[finite_mask]
                    y_fit = _sigmoid_model(x_fit, row['L'], row['k'], row['x0'], row['b'])
                    is_valid = bool(row.get('is_valid', False))

                    if is_valid:
                        status = 'valid'
                        fit_color = '#2E8B57'
                        summary_counts['valid'] += 1
                    else:
                        status = f"invalid: {row.get('reasons', 'quality_failed')}"
                        fit_color = '#C74632'
                        summary_counts['quality_failed'] += 1

                    ax.plot(x_fit, y_fit, color=fit_color, linewidth=2.0, linestyle='--')

                    if np.isfinite(row.get('x0', np.nan)):
                        ax.axvline(row['x0'], color=fit_color, linestyle=':', linewidth=1.2, alpha=0.9)
                else:
                    status = f"fit failed: {row.get('reasons', 'fit_failed')}"
                    summary_counts['fit_failed_or_missing'] += 1
            else:
                summary_counts['fit_failed_or_missing'] += 1

            ax.set_title(f"{animal} | {status}", fontsize=10)
            ax.set_xlabel('Trial')
            ax.set_ylabel('simba_median_balance')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

        for ax in axes_flat[n_animals:]:
            ax.set_visible(False)

        from matplotlib.lines import Line2D
        legend_handles = [
            Line2D([0], [0], marker='o', color='none', markerfacecolor='black', markersize=6, label='Data'),
            Line2D([0], [0], color='#2E8B57', linestyle='--', linewidth=2, label='Valid fit'),
            Line2D([0], [0], color='#C74632', linestyle='--', linewidth=2, label='Invalid fit (quality check failed)'),
            Line2D([0], [0], color='#999999', linestyle='--', linewidth=2, label='Fit failed/missing'),
        ]
        fig.legend(handles=legend_handles, loc='upper center', ncol=2, frameon=False)

        fig.suptitle(
            "simba_median_balance sigmoid fits by animal (deplete + 45NaCl)\n"
            f"valid={summary_counts['valid']} | quality_failed={summary_counts['quality_failed']} | fit_failed_or_missing={summary_counts['fit_failed_or_missing']}",
            y=1.03,
            fontsize=13,
        )
        fig.tight_layout()
        plt.show()

# %%
## Figure 3-style plot: dopamine + behavior realigned to simba_median_balance transition points
if 'simba_median_balance_sigmoid_fits' not in globals() or simba_median_balance_sigmoid_fits.empty:
    print("Run the simba_median_balance fitting cell first.")
elif 'deplete_45_subset' not in globals():
    deplete_45_subset = x_array.query("condition == 'deplete' & infusiontype == '45NaCl'").copy()

if 'simba_median_balance_sigmoid_fits' in globals() and not simba_median_balance_sigmoid_fits.empty:
    transition_fits = (
        simba_median_balance_sigmoid_fits
        .dropna(subset=['x0'])
        .query('is_valid == True')
        .copy()
    )

    if transition_fits.empty:
        print("No valid simba_median_balance transition points available for realignment.")
    else:
        fit_lookup = (
            transition_fits[['animal', 'x0', 'k']]
            .rename(columns={'animal': 'id', 'x0': 'x0_fit', 'k': 'k_fit'})
            .assign(x0_trial=lambda df: df['x0_fit'].round().astype(int))
        )

        realigned_df = (
            deplete_45_subset
            .merge(fit_lookup, on='id', how='inner')
            .assign(trial_aligned=lambda df: df['trial'] - df['x0_trial'])
            .sort_values(['id', 'trial'])
            .copy()
        )

        n_required = int(realigned_df['id'].nunique())
        realigned_complete = (
            realigned_df
            .groupby('trial_aligned', group_keys=False)
            .filter(lambda group: len(group) == n_required)
            .copy()
        )

        def _fit_group_sigmoid(df, group_col, value_col):
            grouped = df.groupby(group_col)[value_col]
            mean = grouped.mean().sort_index()
            sem = grouped.sem().sort_index().fillna(0)
            x_vals = mean.index.to_numpy(dtype=float)
            y_vals = mean.to_numpy(dtype=float)

            popt = None
            y_fit = None
            if len(x_vals) >= 4 and np.nanmax(y_vals) - np.nanmin(y_vals) > 1e-8:
                p0 = [
                    float(np.nanmax(y_vals) - np.nanmin(y_vals)),
                    -1.0,
                    float(np.median(x_vals)),
                    float(np.nanmin(y_vals)),
                ]
                bounds = (
                    [-np.inf, -np.inf, np.min(x_vals), -np.inf],
                    [np.inf, np.inf, np.max(x_vals), np.inf],
                )
                try:
                    popt, _ = curve_fit(_sigmoid_model, x_vals, y_vals, p0=p0, bounds=bounds, maxfev=30000)
                    y_fit = _sigmoid_model(x_vals, *popt)
                except Exception:
                    popt, y_fit = None, None

            return {
                'x': x_vals,
                'mean': y_vals,
                'sem': sem.reindex(mean.index).to_numpy(dtype=float),
                'popt': popt,
                'y_fit': y_fit,
            }

        dopamine_original = _fit_group_sigmoid(realigned_df, 'trial', 'auc_snips')
        dopamine_realigned = _fit_group_sigmoid(realigned_complete, 'trial_aligned', 'auc_snips')
        behavior_original = _fit_group_sigmoid(realigned_df, 'trial', 'simba_median_balance')
        behavior_realigned = _fit_group_sigmoid(realigned_complete, 'trial_aligned', 'simba_median_balance')

        fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.6), sharey=False)

        panel_specs = [
            (axes[0], dopamine_original, behavior_original, 'Original trial order', False),
            (axes[1], dopamine_realigned, behavior_realigned, 'Aligned to simba_median_balance x0', True),
        ]

        behavior_axes = []
        for ax, dopamine_data, behavior_data, title, add_vline in panel_specs:
            ax.plot(dopamine_data['x'], dopamine_data['mean'], linestyle='', marker='o', markersize=4.5, color=colors[3], markerfacecolor='white')
            ax.fill_between(
                dopamine_data['x'],
                dopamine_data['mean'] - dopamine_data['sem'],
                dopamine_data['mean'] + dopamine_data['sem'],
                color=colors[3],
                alpha=0.12,
            )
            if dopamine_data['y_fit'] is not None:
                ax.plot(dopamine_data['x'], dopamine_data['y_fit'], color=colors[3], linestyle='--', linewidth=2)

            ax2 = ax.twinx()
            behavior_axes.append(ax2)
            ax2.plot(behavior_data['x'], behavior_data['mean'], linestyle='', marker='o', markersize=4.5, color=colors[1], markerfacecolor='white')
            ax2.fill_between(
                behavior_data['x'],
                behavior_data['mean'] - behavior_data['sem'],
                behavior_data['mean'] + behavior_data['sem'],
                color=colors[1],
                alpha=0.12,
            )
            if behavior_data['y_fit'] is not None:
                ax2.plot(behavior_data['x'], behavior_data['y_fit'], color=colors[1], linestyle='--', linewidth=2)

            if add_vline:
                ax.axvline(0, color='0.5', linestyle='--', linewidth=1, alpha=0.8)

            ax.set_title(title)
            ax.set_xlabel('Trial' if not add_vline else 'Trial relative to x0')
            ax.tick_params(axis='y', labelcolor=colors[3])
            ax2.tick_params(axis='y', labelcolor=colors[1])
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax2.spines['top'].set_visible(False)

        axes[0].set_ylabel('Dopamine AUC', color=colors[3])
        axes[1].set_ylabel('Dopamine AUC', color=colors[3])
        behavior_axes[0].set_ylabel('simba_median_balance', color=colors[1], rotation=270, labelpad=14)
        behavior_axes[1].set_ylabel('simba_median_balance', color=colors[1], rotation=270, labelpad=14)

        for idx, (dopamine_data, behavior_data) in enumerate([
            (dopamine_original, behavior_original),
            (dopamine_realigned, behavior_realigned),
        ]):
            text_lines = []
            if dopamine_data['popt'] is not None:
                text_lines.append(f"DA k={dopamine_data['popt'][1]:.2f}")
            if behavior_data['popt'] is not None:
                text_lines.append(f"Beh k={behavior_data['popt'][1]:.2f}")
            if text_lines:
                axes[idx].text(0.98, 0.95, '\n'.join(text_lines), transform=axes[idx].transAxes, ha='right', va='top', fontsize=9)

        fig.suptitle(f"Dopamine and behavior realigned to simba_median_balance transition points ({n_required} animals)", y=1.03)
        fig.tight_layout()
        plt.show()

        print(f"Included animals: {sorted(realigned_df['id'].unique())}")
        print(f"Original rows: {len(realigned_df)}")
        print(f"Complete realigned rows: {len(realigned_complete)}")
        if not realigned_complete.empty:
            print(f"Realigned trial range: {int(realigned_complete['trial_aligned'].min())} to {int(realigned_complete['trial_aligned'].max())}")

# %% [markdown]
# # Using signtest with median balance to "binarize" behavioral data

# %%
## Trying with "sign-tested / binarized" data
# Use number of bins to binarize, threshold

nbins = 70/100 # based on sig level after a sign test

def threshold_behav(col, nbins=70/100):
    
    newcol = []
    for val in col:
        # print(val)
        if val > nbins: newcol.append(1)
        elif val < -nbins: newcol.append(-1)
        else: newcol.append(0)
        
    return newcol

def add_signtest_behav(x_array):
    
    return (
        x_array.assign(signtest_behav=lambda df_: threshold_behav(df_.simba_median_balance)
    ))
    
# x_array = add_signtest_behav(x_array)
deplete_45_subset = add_signtest_behav(deplete_45_subset)


# %%
a = np.array([1, 1, 1, 0, 1, -1, -1, 0, 0, 0, 1])
b = np.arange(0, len(a))

def binarize_signtest(data_in):
    
    x = np.arange(0, len(data_in))
    
    a_pos = np.ones(sum((data_in > 0)))
    a_pos_x = x[data_in > 0]

    a_neg = np.zeros(sum((data_in < 1)))
    a_neg_x = x[data_in < 0]

    data_out = np.concatenate([a_pos, a_neg])
    data_out_x = np.concatenate([a_pos_x, a_neg_x])
    
    s = np.argsort(data_out_x)
    
    return data_out[s], data_out_x[s]
    
x, y = binarize_signtest(a)


# %%
## Re-fit sigmoids after binarizing auc_simba (above vs below zero)
# Compare binary fits against the original continuous auc_simba fits.
# For binarized data, very steep k values are acceptable and do not invalidate a fit.

binary_fit_results = []
binary_fit_traces = []

for animal in animals:
    animal_data = (
        deplete_45_subset
        .query("id == @animal")
        .sort_values('trial')
        .copy()
    )

    y_raw = animal_data.signtest_behav
    
    y_binary = y_raw
    
    # trying with binarized signtest (just positive)
    y_binary = (y_raw > 0).astype(float)
    x = np.arange(1, len(y_binary) + 1, dtype=float)
    
    # trying with binarized signtest
    y_binary, x = binarize_signtest(y_raw)
    
    # trying with quantized (3 levels)
    # y_binary = y_raw
    # x = np.arange(1, len(y_binary) + 1, dtype=float)

    if np.unique(y_binary).size < 2:
        print(f"  {animal}: SKIPPED (binary auc_simba has no variation)")
        continue

    p0 = [
        1.0,
        -1.0,
        np.median(x),
        0.0,
    ]
    bounds = (
        [-np.inf, -np.inf, np.min(x), -np.inf],
        [np.inf, np.inf, np.max(x), np.inf],
    )

    try:
        params, pcov = curve_fit(_sigmoid_model, x, y_binary, p0=p0, bounds=bounds, maxfev=30000)
        yhat = _sigmoid_model(x, *params)
        r, p_val = _safe_pearson(y_binary, yhat)
        rmse, aic, aicc, bic = _model_metrics(y_binary, yhat, n_params=4)
        _, _, checks = _sigmoid_quality_checks(x, params, pcov)

        binary_checks = dict(checks)
        binary_checks['k_plausible'] = True
        failed = [name for name, ok in binary_checks.items() if not ok]
        is_valid = len(failed) == 0
        reasons = 'ok' if is_valid else ';'.join(failed)

        L, k, x0, b = params
        binary_fit_results.append({
            'animal': animal,
            'n_trials': len(y_binary),
            'n_positive_trials': int(y_binary.sum()),
            'L': L,
            'k': k,
            'x0': x0,
            'b': b,
            'r': r,
            'p': p_val,
            'rmse': rmse,
            'aicc': aicc,
            'bic': bic,
            'is_valid': is_valid,
            'reasons': reasons,
            'x0_interior': binary_checks['x0_interior'],
            'k_plausible': binary_checks['k_plausible'],
            'ci_finite': binary_checks['ci_finite'],
            'asymptotes_covered': binary_checks['asymptotes_covered'],
        })

        binary_fit_traces.append({
            'animal': animal,
            'x': x,
            'y_binary': y_to_use,
            'y_fit': yhat,
            'is_valid': is_valid,
        })

    except Exception as e:
        print(f"  {animal}: FIT FAILED ({e})")
        binary_fit_results.append({
            'animal': animal,
            'n_trials': len(y_to_use),
            'n_positive_trials': int(y_binary.sum()),
            'L': np.nan,
            'k': np.nan,
            'x0': np.nan,
            'b': np.nan,
            'r': np.nan,
            'p': np.nan,
            'rmse': np.nan,
            'aicc': np.nan,
            'bic': np.nan,
            'is_valid': False,
            'reasons': 'fit_failed',
            'x0_interior': False,
            'k_plausible': True,
            'ci_finite': False,
            'asymptotes_covered': False,
        })

auc_simba_binary_sigmoid_fits = pd.DataFrame(binary_fit_results)

# comparison_df = (
#     time_moving_sigmoid_fits[['animal', 'rmse', 'aicc', 'bic', 'is_valid']]
#     .rename(columns={
#         'rmse': 'rmse_continuous',
#         'aicc': 'aicc_continuous',
#         'bic': 'bic_continuous',
#         'is_valid': 'is_valid_continuous',
#     })
#     .merge(
#         auc_simba_binary_sigmoid_fits[['animal', 'n_positive_trials', 'rmse', 'aicc', 'bic', 'is_valid', 'x0', 'k']],
#         on='animal',
#         how='outer'
#     )
#     .rename(columns={
#         'rmse': 'rmse_binary',
#         'aicc': 'aicc_binary',
#         'bic': 'bic_binary',
#         'is_valid': 'is_valid_binary',
#         'x0': 'x0_binary',
#         'k': 'k_binary',
#     })
# )

# comparison_df['rmse_delta'] = comparison_df['rmse_binary'] - comparison_df['rmse_continuous']
# comparison_df['aicc_delta'] = comparison_df['aicc_binary'] - comparison_df['aicc_continuous']
# comparison_df['bic_delta'] = comparison_df['bic_binary'] - comparison_df['bic_continuous']
# comparison_df['binary_better_rmse'] = comparison_df['rmse_delta'] < 0
# comparison_df['binary_better_aicc'] = comparison_df['aicc_delta'] < 0
# comparison_df['binary_better_bic'] = comparison_df['bic_delta'] < 0

binary_display_df = auc_simba_binary_sigmoid_fits.copy()
for col in ['L', 'k', 'x0', 'b', 'r', 'p', 'rmse', 'aicc', 'bic']:
    if col in binary_display_df.columns:
        binary_display_df[col] = binary_display_df[col].round(4)

# comparison_display_df = comparison_df.copy()
# for col in ['rmse_continuous', 'rmse_binary', 'rmse_delta', 'aicc_continuous', 'aicc_binary', 'aicc_delta', 'bic_continuous', 'bic_binary', 'bic_delta', 'x0_binary', 'k_binary']:
#     if col in comparison_display_df.columns:
#         comparison_display_df[col] = comparison_display_df[col].round(4)

print("\n" + "=" * 80)
print("SIGMOID FIT RESULTS FOR BINARIZED auc_simba (> 0 vs <= 0)")
print("=" * 80)
print(binary_display_df[['animal', 'n_trials', 'n_positive_trials', 'L', 'k', 'x0', 'b', 'r', 'p', 'rmse', 'aicc', 'is_valid', 'reasons']].to_string(index=False))

# print("\n" + "=" * 80)
# print("COMPARISON TO ORIGINAL CONTINUOUS auc_simba FITS")
# print("=" * 80)
# print(comparison_display_df[['animal', 'n_positive_trials', 'is_valid_continuous', 'is_valid_binary', 'rmse_continuous', 'rmse_binary', 'rmse_delta', 'aicc_continuous', 'aicc_binary', 'aicc_delta', 'bic_continuous', 'bic_binary', 'bic_delta']].to_string(index=False))

# print("\nSummary of binary-vs-continuous fit quality:")
# print(f"  Better RMSE: {int(comparison_df['binary_better_rmse'].fillna(False).sum())}/{len(comparison_df)} animals")
# print(f"  Better AICc: {int(comparison_df['binary_better_aicc'].fillna(False).sum())}/{len(comparison_df)} animals")
# print(f"  Better BIC:  {int(comparison_df['binary_better_bic'].fillna(False).sum())}/{len(comparison_df)} animals")

# comparison_df

# %%
auc_simba_binary_sigmoid_fits

# %%
deplete_45_subset.columns

# %%
## Plot binarized auc_simba sigmoid fits and compare fit metrics
if 'binary_fit_traces' not in globals() or len(binary_fit_traces) == 0:
    print("No binarized fit traces available. Run the previous cell first.")
else:
    n_animals = len(binary_fit_traces)
    n_cols = 3
    n_rows = int(np.ceil(n_animals / n_cols))
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 3.5 * n_rows), squeeze=False)
    axes_flat = axes.flatten()
    
    for ax, trace in zip(axes_flat, binary_fit_traces):
        animal = trace['animal']
        row = auc_simba_binary_sigmoid_fits.loc[auc_simba_binary_sigmoid_fits['animal'] == animal].iloc[0]
        ax.scatter(trace['x'], trace['y_binary'], s=28, color='black', alpha=0.75, label='binarized auc_simba')
        ax.plot(trace['x'], trace['y_fit'], color=colors[2], linewidth=2, label='sigmoid fit')
        if np.isfinite(row['x0']):
            ax.axvline(row['x0'], color=colors[1], linestyle='--', linewidth=1.5, alpha=0.8)
        ax.set_title(
            f"{animal} | valid={bool(row['is_valid'])} | x0={row['x0']:.1f}",
            fontsize=10
        )
        ax.set_xlabel('Trial')
        ax.set_ylabel('auc_simba > 0')
        ax.set_ylim(-1.1, 1.1)
        ax.set_yticks([0, 1])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    
    for ax in axes_flat[n_animals:]:
        ax.set_visible(False)
    
    handles, labels = axes_flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', ncol=2, frameon=False)
    fig.suptitle('Binarized auc_simba sigmoid fits by animal', y=1.02, fontsize=14)
    fig.tight_layout()
    plt.show()
    
    if 'comparison_df' in globals() and not comparison_df.empty:
        metrics = [
            ('rmse_continuous', 'rmse_binary', 'RMSE'),
            ('aicc_continuous', 'aicc_binary', 'AICc'),
            ('bic_continuous', 'bic_binary', 'BIC'),
        ]
        
        fig, axes = plt.subplots(1, len(metrics), figsize=(4.5 * len(metrics), 4), squeeze=False)
        axes = axes[0]
        
        for ax, (cont_col, bin_col, label) in zip(axes, metrics):
            plot_df = comparison_df[[cont_col, bin_col, 'animal']].dropna().copy()
            if plot_df.empty:
                ax.set_visible(False)
                continue
            
            lower = min(plot_df[cont_col].min(), plot_df[bin_col].min())
            upper = max(plot_df[cont_col].max(), plot_df[bin_col].max())
            ax.scatter(plot_df[cont_col], plot_df[bin_col], color=colors[0], s=45, alpha=0.85)
            ax.plot([lower, upper], [lower, upper], color='0.5', linestyle='--', linewidth=1)
            for _, row in plot_df.iterrows():
                ax.text(row[cont_col], row[bin_col], row['animal'], fontsize=8, alpha=0.8)
            ax.set_xlabel(f'Continuous {label}')
            ax.set_ylabel(f'Binary {label}')
            ax.set_title(f'Binary vs continuous {label}')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
        
        fig.tight_layout()
        plt.show()
    else:
        print("comparison_df is not available, so metric comparison plots were skipped.")

# %%
## Replot dopamine and behavior aligned to fitted binarized auc_simba x0 values
if 'auc_simba_binary_sigmoid_fits' not in globals():
    print("Run the binarized auc_simba sigmoid fitting cell first.")
else:
    alignment_fits = (
        auc_simba_binary_sigmoid_fits
        .dropna(subset=['x0'])
        .query('is_valid == True')
        .copy()
    )

    if alignment_fits.empty:
        print("No valid binarized fits with finite x0 values are available for realignment.")
    else:
        fit_lookup = (
            alignment_fits[['animal', 'x0']]
            .rename(columns={'animal': 'id'})
            .assign(x0_trial=lambda df: df['x0'].round().astype(int))
        )

        rats_to_exclude = ["PB27", "PB71"]
        realigned_auc_simba = (
            deplete_45_subset
            .query("id not in @rats_to_exclude")
            .merge(fit_lookup, on='id', how='inner')
            .assign(trial_aligned=lambda df: df['trial'] - df['x0_trial'])
            .sort_values(['id', 'trial'])
            .copy()
        )

        n_required = realigned_auc_simba['id'].nunique()
        realigned_complete = (
            realigned_auc_simba
            .groupby('trial_aligned', group_keys=False)
            .filter(lambda group: len(group) == n_required)
            .copy()
        )

        def _fit_group_sigmoid(df, group_col, value_col):
            grouped = df.groupby(group_col)[value_col]
            mean = grouped.mean().sort_index()
            sem = grouped.sem().sort_index().fillna(0)
            x_vals = mean.index.to_numpy(dtype=float)
            y_vals = mean.to_numpy(dtype=float)

            popt = None
            y_fit = None
            if len(x_vals) >= 4 and np.nanmax(y_vals) - np.nanmin(y_vals) > 1e-8:
                p0 = [
                    float(np.nanmax(y_vals) - np.nanmin(y_vals)),
                    -1.0,
                    float(np.median(x_vals)),
                    float(np.nanmin(y_vals)),
                ]
                bounds = (
                    [-np.inf, -np.inf, np.min(x_vals), -np.inf],
                    [np.inf, np.inf, np.max(x_vals), np.inf],
                )
                try:
                    popt, _ = curve_fit(_sigmoid_model, x_vals, y_vals, p0=p0, bounds=bounds, maxfev=30000)
                    y_fit = _sigmoid_model(x_vals, *popt)
                except Exception as exc:
                    print(f"Failed to fit {value_col} by {group_col}: {exc}")

            return {
                'x': x_vals,
                'mean': y_vals,
                'sem': sem.reindex(mean.index).to_numpy(dtype=float),
                'popt': popt,
                'y_fit': y_fit,
            }

        dopamine_original = _fit_group_sigmoid(realigned_auc_simba, 'trial', 'auc_snips')
        dopamine_realigned = _fit_group_sigmoid(realigned_complete, 'trial_aligned', 'auc_snips')
        behaviour_original = _fit_group_sigmoid(realigned_auc_simba, 'trial', 'simba_median_balance')
        behaviour_realigned = _fit_group_sigmoid(realigned_complete, 'trial_aligned', 'simba_median_balance')

        fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.6), sharey=False)

        panel_specs = [
            (
                axes[0],
                dopamine_original,
                behaviour_original,
                'trial',
                'Original trial order',
                None,
            ),
            (
                axes[1],
                dopamine_realigned,
                behaviour_realigned,
                'trial_aligned',
                'Aligned to binarized-fit x0',
                0,
            ),
        ]

        behaviour_axes = []
        for ax, dopamine_data, behaviour_data, x_label, title, vline in panel_specs:
            ax.plot(dopamine_data['x'], dopamine_data['mean'], linestyle='', marker='o', markersize=4.5, color=colors[3], markerfacecolor='white')
            ax.fill_between(
                dopamine_data['x'],
                dopamine_data['mean'] - dopamine_data['sem'],
                dopamine_data['mean'] + dopamine_data['sem'],
                color=colors[3],
                alpha=0.12,
            )
            if dopamine_data['y_fit'] is not None:
                ax.plot(dopamine_data['x'], dopamine_data['y_fit'], color=colors[3], linestyle='--', linewidth=2)

            ax2 = ax.twinx()
            behaviour_axes.append(ax2)
            ax2.plot(behaviour_data['x'], behaviour_data['mean'], linestyle='', marker='o', markersize=4.5, color=colors[1], markerfacecolor='white')
            ax2.fill_between(
                behaviour_data['x'],
                behaviour_data['mean'] - behaviour_data['sem'],
                behaviour_data['mean'] + behaviour_data['sem'],
                color=colors[1],
                alpha=0.12,
            )
            if behaviour_data['y_fit'] is not None:
                ax2.plot(behaviour_data['x'], behaviour_data['y_fit'], color=colors[1], linestyle='--', linewidth=2)

            if vline is not None:
                ax.axvline(vline, color='0.5', linestyle='--', linewidth=1, alpha=0.8)

            ax.set_title(title)
            ax.set_xlabel('Trial' if x_label == 'trial' else 'Trial relative to x0')
            ax.tick_params(axis='y', labelcolor=colors[3])
            ax2.tick_params(axis='y', labelcolor=colors[1])
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax2.spines['top'].set_visible(False)

        axes[0].set_ylabel('Dopamine AUC', color=colors[3])
        axes[1].set_ylabel('Dopamine AUC', color=colors[3])
        behaviour_axes[0].set_ylabel('Behavior AUC', color=colors[1], rotation=270, labelpad=14)
        behaviour_axes[1].set_ylabel('Behavior AUC', color=colors[1], rotation=270, labelpad=14)

        for idx, (dopamine_data, behaviour_data) in enumerate([
            (dopamine_original, behaviour_original),
            (dopamine_realigned, behaviour_realigned),
        ]):
            text_lines = []
            if dopamine_data['popt'] is not None:
                text_lines.append(f"DA k={dopamine_data['popt'][1]:.2f}")
            if behaviour_data['popt'] is not None:
                text_lines.append(f"Beh k={behaviour_data['popt'][1]:.2f}")
            if text_lines:
                axes[idx].text(0.98, 0.95, '\n'.join(text_lines), transform=axes[idx].transAxes, ha='right', va='top', fontsize=9)

        fig.suptitle(f"Dopamine and behavior realigned to binarized auc_simba x0 ({n_required} animals)", y=1.03)
        fig.tight_layout()
        plt.show()

        print(f"Included animals: {sorted(realigned_auc_simba['id'].unique())}")
        print(f"Original trials: {len(realigned_auc_simba)}")
        print(f"Complete realigned trials: {len(realigned_complete)}")
        print(f"Realigned trial range: {int(realigned_complete['trial_aligned'].min())} to {int(realigned_complete['trial_aligned'].max())}")

# %%
n_required

# %%
### Questions and things to try

# add in actual "sign test" binarization
# check whether it is actually doing what we think - should go from 1 to -1 - can I run binarized model on that?
# plot what the results look like after binarization for all rats
# move sigmoidal functions and fits into src file
# how can you get ks to become negative (or positive) even when starting with the opposite sign



# plot realignment in some different conditions
# in all cases exclude any fits that go positive (basically means PB26 and PB71 alwasy excluded)
# 1) raw median balance fits
# 2) binarized using signtest (70%) vals
# 3) quantized (using 1, 0, -1) based on signtest


# %% [markdown]
# ## Analyze Velocity of Change (Derivative) of Behavior
#
# Calculate the derivative of behavior over trials, smooth by a few time bins, and find the point of steepest change for each rat.

# %%
from scipy.ndimage import uniform_filter1d

# Set parameters for smoothing and analysis
smoothing_window = 5  # Number of trials to smooth over (should be odd)

def compute_smoothed_derivative(signal, window_size=3):
    """
    Compute the derivative of a signal and smooth it.
    
    Parameters:
    -----------
    signal : array-like
        The behavior signal (e.g., simba_median_balance values)
    window_size : int
        Size of smoothing window (should be odd for symmetry)
        
    Returns:
    --------
    derivative : array
        Smoothed derivative of the signal
    max_velocity_idx : int
        Index of maximum negative velocity (most negative value)
    max_velocity : float
        Maximum negative velocity value (most negative)
    """
    # Compute raw derivative using differences
    derivative = np.diff(signal, prepend=signal[0])
    
    # Smooth the derivative
    smoothed_deriv = uniform_filter1d(derivative, size=window_size, mode='nearest')
    
    # Find point of most negative velocity (steepest negative change)
    max_velocity_idx = np.argmin(smoothed_deriv)
    max_velocity = smoothed_deriv[max_velocity_idx]
    
    return smoothed_deriv, max_velocity_idx, max_velocity

# Calculate derivatives for each rat and behavior metric
df_dep_45 = x_array.query("condition == 'deplete' & infusiontype == '45NaCl'").copy()

velocity_results = []

for rat in df_dep_45.id.unique():
    rat_data = df_dep_45.loc[df_dep_45.id == rat].copy()
    rat_data = rat_data.sort_values('trial')
    
    # Analyze simba_median_balance
    simba_signal = rat_data['simba_median_balance'].to_numpy()
    deriv_simba, max_idx_simba, max_vel_simba = compute_smoothed_derivative(
        simba_signal, window_size=smoothing_window
    )
    
    velocity_results.append({
        'id': rat,
        'n_trials': len(rat_data),
        'simba_max_velocity_trial': rat_data.iloc[max_idx_simba]['trial'],
        'simba_max_velocity_value': max_vel_simba,
    })

velocity_df = pd.DataFrame(velocity_results)
print("Summary of maximum velocity points by rat:")
print(velocity_df.to_string())
print(f"\nMean trial of max simba_median_balance velocity: {velocity_df['simba_max_velocity_trial'].mean():.1f}")

# %%
# Visualize behavior and velocity of change for each rat
fig, axes = plt.subplots(len(df_dep_45.id.unique()), 2, figsize=(14, 2.5*len(df_dep_45.id.unique())))

if len(df_dep_45.id.unique()) == 1:
    axes = axes.reshape(1, -1)

for row, rat in enumerate(sorted(df_dep_45.id.unique())):
    rat_data = df_dep_45.loc[df_dep_45.id == rat].copy()
    rat_data = rat_data.sort_values('trial')
    
    trial_nums = rat_data['trial'].to_numpy()
    
    # Left plot: simba_median_balance with velocity overlay
    simba_signal = rat_data['simba_median_balance'].to_numpy()
    deriv_simba, max_idx_simba, max_vel_simba = compute_smoothed_derivative(
        simba_signal, window_size=smoothing_window
    )
    
    ax_simba = axes[row, 0]
    ax_simba.plot(trial_nums, simba_signal, marker='o', linestyle='-', 
                   color=colors[0], alpha=0.7, label='SIMBA Median Balance', markersize=4)
    ax_simba.axvline(trial_nums[max_idx_simba], color=colors[2], linestyle='--', 
                      linewidth=2, alpha=0.7, label=f'Max velocity (trial {int(trial_nums[max_idx_simba])})')
    ax_simba.set_ylabel('SIMBA Median Balance')
    ax_simba.set_xlabel('Trial Number')
    ax_simba.set_title(f'Rat {rat}: SIMBA Median Balance')
    ax_simba.legend(fontsize=9)
    sns.despine(ax=ax_simba, offset=5)
    
    # Right plot: velocity of change
    ax_velocity = axes[row, 1]
    ax_velocity.plot(trial_nums, deriv_simba, marker='o', linestyle='-', 
                    color=colors[1], alpha=0.7, label='Smoothed Derivative', markersize=4)
    ax_velocity.axhline(0, color='gray', linestyle='-', linewidth=0.5, alpha=0.5)
    ax_velocity.axvline(trial_nums[max_idx_simba], color=colors[2], linestyle='--', 
                       linewidth=2, alpha=0.7, label=f'Max velocity')
    ax_velocity.scatter(trial_nums[max_idx_simba], max_vel_simba, 
                       color=colors[2], s=100, marker='*', zorder=5, edgecolors='black', linewidth=1)
    ax_velocity.set_ylabel('Rate of Change')
    ax_velocity.set_xlabel('Trial Number')
    ax_velocity.set_title(f'Velocity of Change (smoothing window={smoothing_window})')
    ax_velocity.legend(fontsize=9)
    sns.despine(ax=ax_velocity, offset=5)

fig.suptitle("Behavior Dynamics: Signal and Velocity of Change", y=1.00)
fig.tight_layout()
plt.show()

# %%
# Visualize SIMBA median balance transition points
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

# Plot 1: Transition points across rats
rats = velocity_df['id'].values
x_pos = np.arange(len(rats))

ax1.bar(x_pos, velocity_df['simba_max_velocity_trial'], 
        label='Max Velocity Trial', color=colors[1], alpha=0.8)

ax1.set_xlabel('Rat ID')
ax1.set_ylabel('Trial Number of Maximum Velocity')
ax1.set_title('SIMBA Median Balance: Max Velocity Transition Points')
ax1.set_xticks(x_pos)
ax1.set_xticklabels(rats, rotation=45)
ax1.legend()
sns.despine(ax=ax1, offset=5)

# Plot 2: Distribution of velocity values at transition points
ax2.scatter(velocity_df['simba_max_velocity_trial'], 
           velocity_df['simba_max_velocity_value'],
           s=150, alpha=0.7, color=colors[1], edgecolors='black', linewidth=1)
for idx, row in velocity_df.iterrows():
    ax2.annotate(f"PB{row['id']}", 
                (row['simba_max_velocity_trial'], row['simba_max_velocity_value']),
                fontsize=9, ha='center', va='bottom')

ax2.set_xlabel('Trial Number of Max Velocity')
ax2.set_ylabel('Velocity Magnitude')
ax2.set_title('SIMBA Median Balance: Velocity Magnitude at Transition')
sns.despine(ax=ax2, offset=5)

fig.tight_layout()
plt.show()

# %%
# Explore sensitivity to smoothing window size
window_sizes = [1, 3, 5, 7, 9]
sensitivity_results = []

for window_size in window_sizes:
    for rat in df_dep_45.id.unique():
        rat_data = df_dep_45.loc[df_dep_45.id == rat].copy()
        rat_data = rat_data.sort_values('trial')
        
        simba_signal = rat_data['simba_median_balance'].to_numpy()
        deriv_simba, max_idx_simba, max_vel_simba = compute_smoothed_derivative(
            simba_signal, window_size=window_size
        )
        
        sensitivity_results.append({
            'id': rat,
            'window_size': window_size,
            'max_velocity_trial': rat_data.iloc[max_idx_simba]['trial'],
            'max_velocity_value': max_vel_simba,
        })

sensitivity_df = pd.DataFrame(sensitivity_results)

# Visualize sensitivity to window size
fig, axes = plt.subplots(1, 2, figsize=(13, 4))

# Plot 1: How transition trial varies with window size
for rat in df_dep_45.id.unique():
    rat_sens = sensitivity_df[sensitivity_df['id'] == rat]
    axes[0].plot(rat_sens['window_size'], rat_sens['max_velocity_trial'], 
                marker='o', label=f'PB{rat}', alpha=0.7)

axes[0].set_xlabel('Smoothing Window Size (trials)')
axes[0].set_ylabel('Trial Number of Max Velocity')
axes[0].set_title('Stability of Transition Points vs. Smoothing Window')
axes[0].set_xticks(window_sizes)
axes[0].legend(fontsize=9)
axes[0].grid(True, alpha=0.3)
sns.despine(ax=axes[0], offset=5)

# Plot 2: How velocity magnitude varies with window size
for rat in df_dep_45.id.unique():
    rat_sens = sensitivity_df[sensitivity_df['id'] == rat]
    axes[1].plot(rat_sens['window_size'], rat_sens['max_velocity_value'], 
                marker='s', label=f'PB{rat}', alpha=0.7)

axes[1].set_xlabel('Smoothing Window Size (trials)')
axes[1].set_ylabel('Max Velocity Magnitude')
axes[1].set_title('Velocity Magnitude vs. Smoothing Window')
axes[1].set_xticks(window_sizes)
axes[1].legend(fontsize=9)
axes[1].grid(True, alpha=0.3)
sns.despine(ax=axes[1], offset=5)

fig.tight_layout()
plt.show()

print("\nSensitivity Analysis Summary:")
print("Window size vs. stability of transition point detection")
for rat in df_dep_45.id.unique():
    rat_sens = sensitivity_df[sensitivity_df['id'] == rat]
    trial_range = rat_sens['max_velocity_trial'].max() - rat_sens['max_velocity_trial'].min()
    print(f"  Rat {rat}: Trial range = {trial_range} (min={rat_sens['max_velocity_trial'].min()}, max={rat_sens['max_velocity_trial'].max()})")

# %%
