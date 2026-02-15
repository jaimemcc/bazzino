"""
Assemble all data for Bazzino & Roitman sodium appetite project.

This script runs the full data pipeline:
  1. Assemble behavioural (DLC) data — movement metric snips + x_array
  2. Assemble photometry data — photometry snips + x_array
  3. Equalize photometry and DLC data so they match in length
  4. Spectral clustering of photometry data (labels 0/1 added to x_array)
  5. Compute cluster distances (clusterness, euclidean diff)
  6. Sigmoidal transition fitting (deplete + 45NaCl only)
  7. Combine behaviour and photometry x_arrays, realign trials

Outputs a single pickle: data/assembled_data.pickle
containing: x_array, snips_photo, snips_movement, fits_df, z_dep45

Usage:
    python src/assemble_all_data.py

Configurable parameters are collected at the top of the script in PARAMS dict.
"""

from pathlib import Path
import sys
import os
import warnings

os.environ["OMP_NUM_THREADS"] = "2"  # Avoid sklearn threading warning
warnings.filterwarnings("ignore", category=FutureWarning)

import numpy as np
import pandas as pd
import pickle
import dill
import tdt
import trompy as tp

from scipy.optimize import curve_fit
from scipy.spatial.distance import cdist
from sklearn.decomposition import PCA
from sklearn.cluster import SpectralClustering
from sklearn.metrics import silhouette_score

# Add src to path so we can import local modules
SRC_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SRC_DIR))
from extract_behav_parameters import (
    get_ttls, read_DLC_csv, interpolate_low_likehood,
    calc_bodypart_movement, calc_angular_velocity,
    get_behav_snips, smooth_array,
)

# ──────────────────────────────────────────────────────────────────────
# CONFIGURABLE PARAMETERS — edit these to change defaults
# ──────────────────────────────────────────────────────────────────────
PARAMS = {
    # ── Paths ──
    "data_folder": Path("data"),
    "results_folder": Path("results"),
    "tank_folder": Path("D:/TestData/bazzino/from_paula"),
    # "dlc_folder": Path("D:/TestData/bazzino/output_csv_shuffle4"), #office
    "dlc_folder": Path("C:/Users/jmc010/Data/bazzino/Output DLC shuffle 4 csv files"), #laptop

    # ── Behavioural metric ──
    # Options: "movement" or "angular_velocity"
    "behav_metric": "movement",

    # ── DLC parameters ──
    "dlc_likelihood_threshold": 0.6,
    "dlc_bodyparts": ["r_ear", "l_ear", "head_base"],  # Ears + head_base (matches old pipeline exactly)
    "dlc_smooth_method": "gaussian",  # "gaussian", "moving_avg", "savgol", or None
    "dlc_smooth_window": 10,
    "dlc_zscore_to_baseline": False,

    # ── Photometry parameters ──
    "photo_pre_seconds": 5,
    "photo_post_seconds": 15,
    "photo_bins": 200,

    # ── Conditions to exclude ──
    "conditions_to_exclude": ["thirsty", "replete_exp"],

    # ── Snip parameters for AUC calculation ──
    "auc_start_bin": 50,   # bin index for start of infusion window
    "auc_end_bin": 150,    # bin index for end of infusion window

    # ── Velocity smoothing ──
    "vel_smooth_window": 5,

    # ── Clustering ──
    "num_retained_pcs": 3,
    "n_clusters": None,           # None = auto-select via silhouette sweep; or set to int to force
    "max_n_clusters": 9,          # Maximum number of clusters to test in sweep
    "clustering_affinity": "sigmoid",
    "clustering_assign_labels": "discretize",

    # ── Movement analysis parameters ──
    "normalize_movement": True,              # Normalize to [0,1] by default
    "movement_threshold": 0.02,              # Threshold for normalized movement
    "calculate_raw_movement": True,         # Optional: also analyze raw pixels
    "movement_threshold_raw": 0.5,           # Threshold for raw pixel movement (pixels/frame)

    # ── Sigmoidal transition fit ──
    "transition_condition": "deplete",
    "transition_infusion": "45NaCl",
    "logistic_maxfev": 60000,

    # ── Output ──
    "output_filename": "assembled_data.pickle",

    # ── Caching ──
    # Set these to True to load from cached pickle instead of re-extracting.
    # Cache files are saved automatically after each step.
    # If the cache file doesn't exist, the step runs from scratch regardless.
    "cache_behav": False,          # Skip DLC extraction, load from cache
    "cache_photo": True,          # Skip TDT extraction, load from cache
    "cache_clustering": True,     # Skip PCA + spectral clustering, load from cache
    "cache_transitions": False,    # Skip sigmoidal fitting, load from cache

    # Cache filenames (in data_folder)
    "cache_behav_file": "_cache_behav.pickle",
    "cache_photo_file": "_cache_photo.pickle",
    "cache_clustering_file": "_cache_clustering.pickle",
    "cache_transitions_file": "_cache_transitions.pickle",
}


# ──────────────────────────────────────────────────────────────────────
# HELPER FUNCTIONS
# ──────────────────────────────────────────────────────────────────────

def _condition_map():
    return {
        "Sodium Depleted": "deplete",
        "Sodium Replete": "replete",
        "Sodium Replete Experienced": "replete_exp",
        "Thirsty": "thirsty",
    }

def _infusion_map():
    return {45: "45NaCl", 1: "10NaCl"}


def _load_cache(cache_path, label):
    """Try to load a cached pickle. Returns data or None."""
    if cache_path.exists():
        print(f"  Loading cached {label} from {cache_path}")
        with open(cache_path, "rb") as f:
            cached = dill.load(f)
        if "_cached_at" in cached:
            print(f"    cached at: {cached['_cached_at']}")
        if "_cached_params" in cached:
            p = cached["_cached_params"]
            # Print a few key params so you can verify they match
            summary = {k: v for k, v in p.items()
                       if not k.startswith("cache") and k not in ("data_folder", "results_folder")}
            print(f"    cached with params: {summary}")
        return cached
    else:
        print(f"  No cache found at {cache_path}, running {label} from scratch.")
        return None


def _save_cache(data, cache_path, label, params=None):
    """Save intermediate result to a cache pickle, including params and timestamp."""
    from datetime import datetime
    data["_cached_at"] = datetime.now().isoformat()
    if params is not None:
        # Convert Path objects to strings for cleaner serialization
        data["_cached_params"] = {k: str(v) if isinstance(v, Path) else v
                                   for k, v in params.items()}
    with open(cache_path, "wb") as f:
        dill.dump(data, f)
    print(f"  Cached {label} saved to {cache_path}")


def get_movement_vector(stub, dlc_folder, params):
    """Get movement metric from DLC data for a single session."""
    df = read_DLC_csv(stub, dlc_folder)
    df = interpolate_low_likehood(df, threshold=params["dlc_likelihood_threshold"])
    
    # Get normalized movement (default analysis)
    movement_norm = calc_bodypart_movement(
        df,
        include_bodyparts=params["dlc_bodyparts"],
        smooth_method=params["dlc_smooth_method"],
        smooth_window=params["dlc_smooth_window"],
        normalize=True,
    )
    
    # Optionally also get raw pixel movement
    if params.get("calculate_raw_movement", False):
        movement_raw = calc_bodypart_movement(
            df,
            include_bodyparts=params["dlc_bodyparts"],
            smooth_method=params["dlc_smooth_method"],
            smooth_window=params["dlc_smooth_window"],
            normalize=False,
        )
        return movement_norm, movement_raw
    else:
        return movement_norm


def get_angular_velocity_vector(stub, dlc_folder, params):
    """Get angular velocity from DLC data for a single session."""
    df = read_DLC_csv(stub, dlc_folder)
    df = interpolate_low_likehood(df, threshold=params["dlc_likelihood_threshold"])
    df = calc_angular_velocity(df, rightear="r_ear", leftear="l_ear")
    return df.d_angle_deg


def get_behav_vector(stub, dlc_folder, params):
    """Dispatch to the correct behavioural metric function."""
    if params["behav_metric"] == "movement":
        return get_movement_vector(stub, dlc_folder, params)
    elif params["behav_metric"] == "angular_velocity":
        return get_angular_velocity_vector(stub, dlc_folder, params)
    else:
        raise ValueError(f"Unknown behav_metric: {params['behav_metric']}")


# ──────────────────────────────────────────────────────────────────────
# STEP 1 & 2: Assemble photometry and behavioural snips
# ──────────────────────────────────────────────────────────────────────

def get_photometry_snips(tank, params):
    """Extract photometry snips from a TDT tank."""
    data = tdt.read_block(tank)
    blue = data.streams["x65A"].data
    uv = data.streams["x05A"].data
    fs = data.streams["x05A"].fs

    filtered_sig = tp.processdata(blue, uv, fs=fs)
    sol = data.epocs.sol_.onset
    snips = tp.snipper(
        filtered_sig, sol, fs=fs,
        pre=params["photo_pre_seconds"],
        post=params["photo_post_seconds"],
        bins=params["photo_bins"],
    )[0]
    return snips


def assemble_photometry(params):
    """Step 2: Assemble photometry snips and x_array from TDT tanks."""
    print("\n" + "=" * 60)
    print("STEP 2: Assembling photometry data from TDT tanks")
    print("=" * 60)

    data_folder = params["data_folder"]
    tank_folder = params["tank_folder"]

    meta_10 = pd.read_csv(data_folder / "10NaCl_FileKey.csv")
    meta_45 = pd.read_csv(data_folder / "45NaCl_FileKey.csv")

    def _assemble(metadata, infusion_label):
        snips_list, x_list = [], []
        for _, row in metadata.iterrows():
            tank = tank_folder / row["Folder"]
            print(f"  Photometry: {row['Folder']}", end=" ... ")
            try:
                s = get_photometry_snips(tank, params)
                n = len(s)
                print(f"{n} trials")
                snips_list.append(s)
                x_list.append(pd.DataFrame({
                    "trial": np.arange(n),
                    "id": row["Subject"],
                    "condition": row["Physiological state"],
                    "infusiontype": infusion_label,
                }))
            except Exception as e:
                print(f"ERROR: {e}")
        return snips_list, x_list

    snips_10, x_10 = _assemble(meta_10, "10NaCl")
    snips_45, x_45 = _assemble(meta_45, "45NaCl")

    snips_all = np.vstack([np.concatenate(snips_10), np.concatenate(snips_45)])
    x_all = (
        pd.concat(x_10 + x_45, ignore_index=True)
        .replace({"condition": _condition_map()})
    )

    # Add sex info
    subjects_df = (
        pd.concat([
            pd.read_csv(data_folder / "10NaCl_SubjectKey.csv").iloc[:, :2],
            pd.read_csv(data_folder / "45NaCl_SubjectKey.csv").iloc[:, :2],
        ])
        .drop_duplicates()
        .rename(columns={"Subject": "id", "Sex": "sex"})
    )
    x_all = pd.merge(x_all, subjects_df, on="id", how="left")

    # Filter conditions
    mask = ~x_all.condition.isin(params["conditions_to_exclude"])
    snips_all = snips_all[mask.values]
    x_all = x_all[mask].reset_index(drop=True)

    print(f"  Photometry assembled: {snips_all.shape[0]} trials, {snips_all.shape[1]} bins")
    return snips_all, x_all


def assemble_behaviour(params):
    """Step 1: Assemble behavioural (DLC) snips and x_array."""
    print("\n" + "=" * 60)
    print(f"STEP 1: Assembling behavioural data ({params['behav_metric']})")
    print("=" * 60)

    data_folder = params["data_folder"]
    dlc_folder = params["dlc_folder"]

    meta_df = pd.concat([
        pd.read_csv(data_folder / "10NaCl_FileKey.csv"),
        pd.read_csv(data_folder / "45NaCl_FileKey.csv"),
    ])

    snips_list, snips_list_raw, x_list = [], [], []
    for _, row in meta_df.iterrows():
        stub = row["Folder"]
        print(f"  Behaviour: {stub}", end=" ... ")
        try:
            behav_vec_result = get_behav_vector(stub, dlc_folder, params)
            
            # Handle case where get_behav_vector returns tuple (norm, raw) or single vector
            if isinstance(behav_vec_result, tuple):
                behav_vec, behav_vec_raw = behav_vec_result
            else:
                behav_vec = behav_vec_result
                behav_vec_raw = None
            
            solenoid_ts = get_ttls(stub, data_folder)
            s = get_behav_snips(
                solenoid_ts=solenoid_ts,
                behav_vector=behav_vec,
                zscore_to_baseline=params["dlc_zscore_to_baseline"],
            )
            n = len(s)
            print(f"{n} trials")
            snips_list.append(s)
            
            # Also create snips for raw movement if available
            if behav_vec_raw is not None:
                s_raw = get_behav_snips(
                    solenoid_ts=solenoid_ts,
                    behav_vector=behav_vec_raw,
                    zscore_to_baseline=False,  # Don't zscore raw pixels
                )
                snips_list_raw.append(s_raw)

            infusion_label = "45NaCl" if row["TreatNum"] == 45 else "10NaCl"
            x_list.append(pd.DataFrame({
                "trial": np.arange(n),
                "id": row["Subject"],
                "condition": row["Physiological state"],
                "infusiontype": infusion_label,
            }))
        except Exception as e:
            print(f"ERROR: {e}")

    snips_all = np.concatenate(snips_list)
    snips_all_raw = np.concatenate(snips_list_raw) if snips_list_raw else None
    x_all = (
        pd.concat(x_list, ignore_index=True)
        .replace({"condition": _condition_map()})
    )

    # Add sex info
    subjects_df = (
        pd.concat([
            pd.read_csv(data_folder / "10NaCl_SubjectKey.csv").iloc[:, :2],
            pd.read_csv(data_folder / "45NaCl_SubjectKey.csv").iloc[:, :2],
        ])
        .drop_duplicates()
        .rename(columns={"Subject": "id", "Sex": "sex"})
    )
    x_all = pd.merge(x_all, subjects_df, on="id", how="left")

    # Filter conditions
    mask = ~x_all.condition.isin(params["conditions_to_exclude"])
    snips_all = snips_all[mask.values]
    if snips_all_raw is not None:
        snips_all_raw = snips_all_raw[mask.values]
    x_all = x_all[mask].reset_index(drop=True)

    print(f"  Behaviour assembled: {snips_all.shape[0]} trials, {snips_all.shape[1]} bins")
    if snips_all_raw is not None:
        print(f"  Raw movement snips also assembled")
    return snips_all, snips_all_raw, x_all


# ──────────────────────────────────────────────────────────────────────
# STEP 3: Equalize and combine
# ──────────────────────────────────────────────────────────────────────

def equalize_datasets(x_photo, snips_photo, x_behav, snips_behav, snips_behav_raw=None):
    """
    Step 3: Make sure photometry and behaviour datasets match row-for-row.
    Finds common rows based on (trial, id, condition, infusiontype) and
    removes extras from whichever is longer.
    """
    print("\n" + "=" * 60)
    print("STEP 3: Equalizing photometry and behaviour datasets")
    print("=" * 60)

    merge_cols = ["trial", "id", "condition", "infusiontype"]
    df_p = x_photo[merge_cols].reset_index(drop=True)
    df_b = x_behav[merge_cols].reset_index(drop=True)

    if df_p.equals(df_b):
        print("  Datasets already aligned.")
        return x_photo, snips_photo, x_behav, snips_behav, snips_behav_raw

    # Find common rows via inner merge
    merged = pd.merge(df_p.assign(_idx_p=df_p.index),
                       df_b.assign(_idx_b=df_b.index),
                       on=merge_cols, how="inner")

    idx_p = merged["_idx_p"].values
    idx_b = merged["_idx_b"].values

    x_photo = x_photo.iloc[idx_p].reset_index(drop=True)
    snips_photo = snips_photo[idx_p]
    x_behav = x_behav.iloc[idx_b].reset_index(drop=True)
    snips_behav = snips_behav[idx_b]

    if snips_behav_raw is not None:
        snips_behav_raw = snips_behav_raw[idx_b]

    print(f"  After equalization: {len(x_photo)} trials (photo), {len(x_behav)} trials (behav)")
    assert len(x_photo) == len(x_behav), "Datasets still not aligned!"
    return x_photo, snips_photo, x_behav, snips_behav, snips_behav_raw


# ──────────────────────────────────────────────────────────────────────
# STEP 4: Clustering of photometry data
# ──────────────────────────────────────────────────────────────────────

def cluster_photometry(snips_photo, x_array, params):
    """
    Step 4: PCA + spectral clustering on photometry snips.
    Adds cluster_photo column (0 or 1) to x_array.
    Returns x_array, pca_transformed data.
    """
    print("\n" + "=" * 60)
    print("STEP 4: Spectral clustering of photometry data")
    print("=" * 60)

    num_pcs = params["num_retained_pcs"]

    # PCA
    pca = PCA(n_components=snips_photo.shape[1], whiten=True)
    pca.fit(snips_photo)
    transformed = pca.transform(snips_photo)

    x = 100 * pca.explained_variance_ratio_
    xprime = x - (x[0] + (x[-1] - x[0]) / (x.size - 1) * np.arange(x.size))
    auto_pcs = np.argmin(xprime)
    print(f"  Auto-detected PCs to keep: {auto_pcs}, using: {num_pcs}")

    # Determine number of clusters
    n_clusters = params["n_clusters"]
    if n_clusters is None:
        # Sweep through cluster numbers and pick best silhouette score
        max_k = params["max_n_clusters"]
        possible_n = np.arange(2, max_k + 1)
        sil_scores = np.nan * np.ones(possible_n.size)

        for idx, k in enumerate(possible_n):
            m = SpectralClustering(
                n_clusters=k,
                affinity=params["clustering_affinity"],
                assign_labels=params["clustering_assign_labels"],
                random_state=0,
            )
            m.fit(transformed[:, :num_pcs])
            sil_scores[idx] = silhouette_score(
                transformed[:, :num_pcs], m.labels_, metric="cosine"
            )
            print(f"    k={k}: silhouette={sil_scores[idx]:.3f}")

        best_idx = np.nanargmax(sil_scores)
        n_clusters = int(possible_n[best_idx])
        print(f"  Best silhouette at k={n_clusters} (score={sil_scores[best_idx]:.3f})")
    else:
        print(f"  Using fixed n_clusters={n_clusters}")

    # Final clustering with chosen n_clusters
    model = SpectralClustering(
        n_clusters=n_clusters,
        affinity=params["clustering_affinity"],
        assign_labels=params["clustering_assign_labels"],
        random_state=0,
    )
    model.fit(transformed[:, :num_pcs])
    sil = silhouette_score(transformed[:, :num_pcs], model.labels_, metric="cosine")
    print(f"  Final clustering: k={n_clusters}, silhouette={sil:.3f}")

    # Reorder clusters so cluster 0 = most positive response during infusion
    # This matches the reorder_clusters function from spectral_clustering_all_trials.ipynb
    pre_window = int(params["photo_pre_seconds"] * 10)  # 50 bins
    uniquelabels = list(set(model.labels_))
    responses = np.nan * np.ones((len(uniquelabels),))
    for l, label in enumerate(uniquelabels):
        responses[l] = np.mean(snips_photo[model.labels_ == label, pre_window:2 * pre_window])
    temp = np.argsort(responses).astype(int)[::-1]
    temp = np.array([np.where(temp == a)[0][0] for a in uniquelabels])
    newlabels = np.array([temp[a] for a in list(np.digitize(model.labels_, uniquelabels) - 1)])

    x_array = x_array.assign(cluster_photo=newlabels)
    print(f"  Cluster counts: {dict(zip(*np.unique(newlabels, return_counts=True)))}")

    return x_array, transformed


# ──────────────────────────────────────────────────────────────────────
# STEP 5: Cluster distances (clusterness + euclidean)
# ──────────────────────────────────────────────────────────────────────

def compute_cluster_distances(x_array, pca_data, params):
    """
    Step 5: Compute clusterness (projection onto cluster vector) and
    euclidean distance difference for each trial.
    """
    print("\n" + "=" * 60)
    print("STEP 5: Computing cluster distances")
    print("=" * 60)

    num_pcs = params["num_retained_pcs"]
    pca_subset = pca_data[:, :num_pcs]

    centroid_0 = pca_subset[x_array.cluster_photo == 0].mean(axis=0)
    centroid_1 = pca_subset[x_array.cluster_photo == 1].mean(axis=0)

    # Clusterness (projection)
    cluster_vector = centroid_0 - centroid_1
    cluster_vector_norm = cluster_vector / np.linalg.norm(cluster_vector)
    projections = np.dot(pca_subset - centroid_1, cluster_vector_norm)
    normalized = (projections - projections.min()) / (projections.max() - projections.min())
    x_array = x_array.assign(clusterness_photo=normalized)

    # Euclidean distance difference
    centroids = np.vstack([centroid_0, centroid_1])
    distances = cdist(pca_subset, centroids, metric="euclidean")
    x_array = x_array.assign(euclidean_diff=distances[:, 1] - distances[:, 0])

    print("  Added clusterness_photo and euclidean_diff columns.")
    return x_array


# ──────────────────────────────────────────────────────────────────────
# STEP 6: Sigmoidal transition fitting
# ──────────────────────────────────────────────────────────────────────

def _logistic4(x, A, L, x0, k):
    return A + (L - A) / (1 + np.exp(-k * (x - x0)))

def _logistic3(x, L, x0, k):
    return L / (1 + np.exp(-k * (x - x0)))


def fit_logistic_per_series(y, x=None, prefer_4p=True, direction=None, maxfev=60000):
    """Fit a logistic curve to binary/near-binary data with robust inits."""
    if x is None:
        x = np.arange(len(y), dtype=float)
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    # Normalize x
    x_mean, x_std = float(np.mean(x)), float(np.std(x))
    if not np.isfinite(x_std) or x_std == 0:
        x_std = 1.0
    x_norm = (x - x_mean) / x_std
    y_clip = np.clip(y, 1e-4, 1 - 1e-4)

    y_min, y_max = float(np.min(y_clip)), float(np.max(y_clip))
    A_init, L_init, x0_init = y_min, y_max, 0.0

    if direction is None:
        try:
            c = float(np.corrcoef(x, y_clip)[0, 1])
        except Exception:
            c = 0.0
        if not np.isfinite(c):
            c = 0.0
        sign = 1.0 if c >= 0 else -1.0
    else:
        sign = 1.0 if direction == "increasing" else -1.0

    k_mags = [0.5, 1.0, 2.0]

    def try_fit(func, p0_list, bounds):
        best, best_rss = None, np.inf
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

    res4 = None
    if prefer_4p:
        p0s_4 = [[A_init, L_init, x0_init, sign * km] for km in k_mags]
        bnds_4 = ([-0.1, 0.4, -3.0, -10.0], [0.6, 1.6, 3.0, 10.0])
        res4 = try_fit(_logistic4, p0s_4, bnds_4)

    p0s_3 = [[L_init, x0_init, sign * km] for km in k_mags]
    bnds_3 = ([0.4, -3.0, -10.0], [1.6, 3.0, 10.0])
    res3 = try_fit(_logistic3, p0s_3, bnds_3)

    if res4 is not None:
        popt, y_hat = res4
        A, L, x0n, k = map(float, popt)
        x0_orig = x0n * x_std + x_mean
        # Calculate R-squared
        ss_res = float(np.sum((y_clip - y_hat) ** 2))
        ss_tot = float(np.sum((y_clip - np.mean(y_clip)) ** 2))
        r_squared = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 0.0
        return {"model": "logistic4", "x0_orig": x0_orig, "k": k, "r_squared": r_squared,
                "params": {"A": A, "L": L, "x0_norm": x0n, "x0_orig": x0_orig, "k": k},
                "y_hat": y_hat, "success": True, "note": ""}
    elif res3 is not None:
        popt, y_hat = res3
        L, x0n, k = map(float, popt)
        x0_orig = x0n * x_std + x_mean
        # Calculate R-squared
        ss_res = float(np.sum((y_clip - y_hat) ** 2))
        ss_tot = float(np.sum((y_clip - np.mean(y_clip)) ** 2))
        r_squared = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 0.0
        # For 3p model, A is implicitly 0 (the function is A + (L-A)/(1+exp(...)) with A=0)
        A = 0.0
        return {"model": "logistic3", "x0_orig": x0_orig, "k": k, "r_squared": r_squared,
                "params": {"A": A, "L": L, "x0_norm": x0n, "x0_orig": x0_orig, "k": k},
                "y_hat": y_hat, "success": True, "note": "4p failed; used 3p"}
    else:
        return {"model": None, "x0_orig": np.nan, "k": np.nan, "r_squared": np.nan,
                "params": {}, "y_hat": None, "success": False, "note": "fit failed"}


def find_sigmoidal_transitions(x_array, params):
    """
    Step 6: Fit sigmoidal transitions per rat for deplete + 45NaCl.
    Uses raw cluster assignments (binary).
    Returns fits_df with transition points.
    """
    print("\n" + "=" * 60)
    print("STEP 6: Finding sigmoidal transitions")
    print("=" * 60)

    cond = params["transition_condition"]
    inf = params["transition_infusion"]
    df_dep_45 = x_array.query("condition == @cond & infusiontype == @inf").copy()

    all_fits = []
    for rat in df_dep_45.id.unique():
        sig = df_dep_45.loc[df_dep_45.id == rat, "cluster_photo"].to_numpy()
        y = np.logical_not(sig).astype(int)
        x = np.arange(len(y), dtype=float)

        fit = fit_logistic_per_series(y, x=x, prefer_4p=True, direction="decreasing",
                                      maxfev=params["logistic_maxfev"])
        all_fits.append({
            "id": rat,
            **fit["params"],
            "model": fit["model"],
            "x0_orig": fit["x0_orig"],
            "k": fit["k"],
            "r_squared": fit["r_squared"],
            "success": fit["success"],
            "note": fit["note"],
        })

    fits_df = pd.DataFrame(all_fits)
    fits_df = fits_df.query("success == True and x0_orig > 0").copy()

    print(f"  Successful fits: {len(fits_df)} / {len(df_dep_45.id.unique())} rats")
    print(f"  Transition points: {fits_df.x0_orig.round(1).tolist()}")
    return fits_df


# ──────────────────────────────────────────────────────────────────────
# STEP 7: Combine and realign
# ──────────────────────────────────────────────────────────────────────

def get_time_moving(snips, threshold=0.02, start_bin=50, end_bin=150):
    """Calculate proportion of time spent moving per trial."""
    moving = []
    for i in range(snips.shape[0]):
        snip = snips[i, start_bin:end_bin]
        tmp = len([x for x in snip if x > threshold]) / len(snip)
        moving.append(tmp)
    return np.array(moving)


def combine_and_realign(x_photo, snips_photo, snips_behav, fits_df, params, snips_behav_raw=None):
    """
    Step 7: Add AUCs and time_moving to x_array, create realigned deplete+45NaCl subset.
    """
    print("\n" + "=" * 60)
    print("STEP 7: Combining data and realigning trials")
    print("=" * 60)

    # Smooth velocity snips
    snips_vel_smooth = smooth_array(snips_behav, window_size=params["vel_smooth_window"])

    # Calculate AUCs using trapezoidal rule (true area under curve)
    s, e = params["auc_start_bin"], params["auc_end_bin"]
    auc_snips = np.array([np.trapz(snips_photo[i, s:e]) for i in range(len(snips_photo))])
    auc_vel = np.array([np.trapz(snips_vel_smooth[i, s:e]) for i in range(len(snips_vel_smooth))])

    # Calculate time moving (normalized)
    time_moving = get_time_moving(snips_behav, threshold=params["movement_threshold"],
                                   start_bin=s, end_bin=e)

    # Calculate time moving raw (if available)
    time_moving_raw = None
    if snips_behav_raw is not None:
        time_moving_raw = get_time_moving(snips_behav_raw, threshold=params["movement_threshold_raw"],
                                          start_bin=s, end_bin=e)

    x_combined = x_photo.assign(
        auc_snips=auc_snips,
        auc_vel=auc_vel,
        time_moving=time_moving,
    )
    
    # Add time_moving_raw if available
    if time_moving_raw is not None:
        x_combined = x_combined.assign(time_moving_raw=time_moving_raw)

    # Add trial_aligned column for ALL trials (NaN for rats without fits)
    trial_aligned = []
    for idx, row in x_combined.iterrows():
        rat_id = row['id']
        if rat_id not in fits_df.id.values:
            # No fit for this rat
            trial_aligned.append(np.nan)
        else:
            # Rat has a fit - calculate relative trial number
            transition = int(fits_df.query("id == @rat_id").x0_orig.values[0])
            trial_aligned.append(row['trial'] - transition)
    
    x_combined = x_combined.assign(trial_aligned=trial_aligned)

    # Now realign the deplete+45NaCl subset (with NaN values dropped)
    cond = params["transition_condition"]
    inf = params["transition_infusion"]
    z = (
        x_combined
        .query("condition == @cond & infusiontype == @inf")
        .dropna(subset=['trial_aligned'])
        .reset_index(drop=True)
    )

    print(f"  Combined x_array: {len(x_combined)} trials")
    print(f"  Added trial_aligned column (NaN for {(x_combined['trial_aligned'].isna()).sum()} trials without fits)")
    if time_moving_raw is not None:
        print(f"  Added time_moving_raw column (threshold={params['movement_threshold_raw']} pixels)")
    print(f"  Realigned deplete+45NaCl subset: {len(z)} trials with valid alignments")

    return x_combined, z


# ──────────────────────────────────────────────────────────────────────
# MAIN PIPELINE
# ──────────────────────────────────────────────────────────────────────

def run_pipeline(params=None):
    """Run the full data assembly pipeline."""
    if params is None:
        params = PARAMS

    print("=" * 60)
    print("BAZZINO DATA ASSEMBLY PIPELINE")
    print("=" * 60)
    print(f"  Behaviour metric: {params['behav_metric']}")
    print(f"  DLC folder: {params['dlc_folder']}")
    print(f"  Tank folder: {params['tank_folder']}")
    print(f"  Caching: behav={params['cache_behav']}, photo={params['cache_photo']}, "
          f"clustering={params['cache_clustering']}, transitions={params['cache_transitions']}")

    data_folder = params["data_folder"]

    # Step 1: Behaviour
    behav_cache_path = data_folder / params["cache_behav_file"]
    if params["cache_behav"]:
        cached = _load_cache(behav_cache_path, "behaviour")
    else:
        cached = None
    if cached is not None:
        snips_behav = cached["snips_behav"]
        snips_behav_raw = cached.get("snips_behav_raw", None)
        x_behav = cached["x_behav"]
        print(f"  Behaviour from cache: {snips_behav.shape[0]} trials")
    else:
        snips_behav, snips_behav_raw, x_behav = assemble_behaviour(params)
        cache_dict = {"snips_behav": snips_behav, "x_behav": x_behav}
        if snips_behav_raw is not None:
            cache_dict["snips_behav_raw"] = snips_behav_raw
        _save_cache(cache_dict, behav_cache_path, "behaviour", params)

    # Step 2: Photometry
    photo_cache_path = data_folder / params["cache_photo_file"]
    if params["cache_photo"]:
        cached = _load_cache(photo_cache_path, "photometry")
    else:
        cached = None
    if cached is not None:
        snips_photo, x_photo = cached["snips_photo"], cached["x_photo"]
        print(f"  Photometry from cache: {snips_photo.shape[0]} trials")
    else:
        snips_photo, x_photo = assemble_photometry(params)
        _save_cache({"snips_photo": snips_photo, "x_photo": x_photo},
                    photo_cache_path, "photometry", params)

    # Step 3: Equalize
    x_photo, snips_photo, x_behav, snips_behav, snips_behav_raw = equalize_datasets(
        x_photo, snips_photo, x_behav, snips_behav, snips_behav_raw
    )

    # Step 4: Clustering
    clustering_cache_path = data_folder / params["cache_clustering_file"]
    if params["cache_clustering"]:
        cached = _load_cache(clustering_cache_path, "clustering")
    else:
        cached = None
    if cached is not None:
        x_combined, pca_transformed = cached["x_combined"], cached["pca_transformed"]
        print(f"  Clustering from cache: {len(x_combined)} trials")
    else:
        x_combined, pca_transformed = cluster_photometry(snips_photo, x_photo, params)
        # Step 5 always runs with fresh clustering
        x_combined = compute_cluster_distances(x_combined, pca_transformed, params)
        _save_cache({"x_combined": x_combined, "pca_transformed": pca_transformed},
                    clustering_cache_path, "clustering", params)

    # Step 5: Cluster distances (only if loaded from cache, otherwise already done above)
    if cached is not None:
        x_combined = compute_cluster_distances(x_combined, pca_transformed, params)

    # Step 6: Sigmoidal transitions
    # NEW: Always load from sigmoidal_fits.pickle (the authoritative source)
    sigmoidal_fits_path = data_folder / "sigmoidal_fits.pickle"
    if sigmoidal_fits_path.exists():
        print(f"\n  Loading transitions from sigmoidal_fits.pickle (authoritative source)")
        with open(sigmoidal_fits_path, "rb") as f:
            fits_data = dill.load(f)
        fits_df = fits_data.get("fits_df_cluster_raw")
        if fits_df is not None:
            print(f"  Loaded {len(fits_df)} transition fits from sigmoidal_fits.pickle")
        else:
            print(f"  WARNING: fits_df_cluster_raw not found, calculating from scratch")
            fits_df = find_sigmoidal_transitions(x_combined, params)
    else:
        # Fallback: use cache or calculate
        transitions_cache_path = data_folder / params["cache_transitions_file"]
        if params["cache_transitions"]:
            cached = _load_cache(transitions_cache_path, "transitions")
        else:
            cached = None
        if cached is not None:
            fits_df = cached["fits_df"]
            print(f"  Transitions from cache: {len(fits_df)} fits")
        else:
            print(f"  WARNING: sigmoidal_fits.pickle not found, calculating transitions")
            fits_df = find_sigmoidal_transitions(x_combined, params)
            _save_cache({"fits_df": fits_df}, transitions_cache_path, "transitions", params)

    # Step 7: Combine and realign
    x_combined, z_dep45 = combine_and_realign(
        x_combined, snips_photo, snips_behav, fits_df, params, snips_behav_raw
    )

    # Create metadata about data processing
    metadata = {
        "behav_smoothed": params["dlc_smooth_method"] is not None,
        "behav_smooth_method": params["dlc_smooth_method"],
        "behav_smooth_window": params["dlc_smooth_window"] if params["dlc_smooth_method"] is not None else None,
        "behav_zscored": params["dlc_zscore_to_baseline"],
        "behav_bodyparts": params.get("dlc_bodyparts"),
        "dlc_likelihood_threshold": params.get("dlc_likelihood_threshold", 0.6),
        "photo_smoothed": False,  # Photometry is NOT smoothed during assembly
        "photo_zscored": True,  # Photometry is z-scored by trompy during processing
        "behav_metric": params["behav_metric"],
        # Movement analysis parameters
        "normalize_movement": params.get("normalize_movement", True),
        "movement_threshold": params.get("movement_threshold", 0.02),
        "calculate_raw_movement": params.get("calculate_raw_movement", False),
        "movement_threshold_raw": params.get("movement_threshold_raw", 0.5) if params.get("calculate_raw_movement", False) else None,
    }

    # Save output
    output = {
        "x_array": x_combined,
        "snips_photo": snips_photo,
        "snips_behav": snips_behav,
        "pca_transformed": pca_transformed,
        "fits_df": fits_df,
        "z_dep45": z_dep45,
        "metadata": metadata,
        "params": params,
    }
    
    # Add optional raw movement data
    if snips_behav_raw is not None:
        output["snips_behav_raw"] = snips_behav_raw

    output_path = params["data_folder"] / params["output_filename"]
    with open(output_path, "wb") as f:
        dill.dump(output, f)

    print("\n" + "=" * 60)
    print(f"DONE! Saved to {output_path}")
    print("=" * 60)
    print(f"  x_array shape:     {x_combined.shape}")
    print(f"  snips_photo shape: {snips_photo.shape}")
    print(f"  snips_behav shape: {snips_behav.shape}")
    print(f"  z_dep45 shape:     {z_dep45.shape}")
    print(f"  fits_df shape:     {fits_df.shape}")
    print(f"\nData processing metadata:")
    print(f"  Behaviour metric:    {metadata['behav_metric']}")
    print(f"  Behaviour smoothed:  {metadata['behav_smoothed']} (method: {metadata['behav_smooth_method']}, window: {metadata['behav_smooth_window']})")
    print(f"  Behaviour z-scored:  {metadata['behav_zscored']}")
    print(f"  Photometry smoothed: {metadata['photo_smoothed']}")
    print(f"  Photometry z-scored: {metadata['photo_zscored']}")

    return output


if __name__ == "__main__":
    # Change to project root so relative paths work
    project_root = SRC_DIR.parent
    os.chdir(project_root)

    run_pipeline(PARAMS)
