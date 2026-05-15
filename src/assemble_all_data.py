"""
Assemble all data for Bazzino & Roitman sodium appetite project.

This script runs the full data pipeline:
  1. Assemble behavioural (DLC) data — movement metric snips + x_array
  2. Assemble photometry data — photometry snips + x_array
    2b. Assemble Simba appetitive probability snips + x_array
    3. Equalize photometry, DLC, and Simba datasets so they match row-for-row
    4. Spectral clustering of photometry data (labels 0/1 added to x_array)
    5. Compute cluster distances (clusterness, euclidean diff)
    6. Sigmoidal transition fitting (deplete + 45NaCl only)
        7. Combine behaviour and photometry x_arrays, realign trials

Outputs a single pickle: data/assembled_data.pickle
containing: x_array, snips_photo, snips_simba_zscore,
            snips_simba_shifted_baseline, snips_simba_raw,
            snips_movement, fits_df, fits_df_da, fits_df_behav, z_dep45

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

from scipy.spatial.distance import cdist
from scipy.ndimage import uniform_filter1d
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
from extract_simba import (
    read_simba_probability,
    get_cam1_onsets_for_stub,
    make_simba_snips,
    get_shifted_snip_means,
    get_shifted_snip_means_multi_shift,
    get_time_above_simba_ci,
)
from model_fit_helpers import (
    fit_logistic_transitions_by_id,
    fit_continuous_sigmoid_transitions_by_id,
    build_realigned_trials,
    binarize_series_threshold,
    quantize_series_threshold,
)

# ──────────────────────────────────────────────────────────────────────
# CONFIGURABLE PARAMETERS — edit these to change defaults
# ──────────────────────────────────────────────────────────────────────
PARAMS = {
    # ── Paths ──
    "data_folder": Path("data"),
    "results_folder": Path("results"),
    # "tank_folder": Path("D:/TestData/bazzino/from_paula"),
    "tank_folder": Path("C:/Users/jmc010/Data/bazzino/tanks"), #laptop
    "dlc_folder": Path("D:/TestData/bazzino/output_csv_shuffle4"), #office
    # "dlc_folder": Path("C:/Users/jmc010/Data/bazzino/Output DLC shuffle 4 csv files"), #laptop
    # "simba_folder": Path("D:/TestData/bazzino/simba_preds/Output_all_animals_appetitive"), #office
    "simba_folder": Path("C:/Users/jmc010/Data/bazzino/simba"), #laptop

    # ── Behavioural metric ──
    # NOTE: Both movement and angular_velocity are now always calculated
    "behav_metric": "movement",  # Deprecated - kept for reference only

    # ── DLC parameters ──
    "dlc_likelihood_threshold": 0.6,
    "dlc_bodyparts": ["r_ear", "l_ear", "head_base"],  # Ears + head_base (matches old pipeline exactly)
    "dlc_smooth_method": "gaussian",  # "gaussian", "moving_avg", "savgol", or None
    "dlc_smooth_window": 5,
    "dlc_zscore_to_baseline": False,
    
    # ── SIMBA parameters ──
    "simba_probability_column": "Probability_Appetitive",
    "simba_fps": 10,
    "simba_use_cam1_timestamps": True,
    "simba_pre_bins": 50,
    "simba_post_bins": 150,
    "simba_use_shifted_baseline": True,
    "simba_shift_frames": 300,
    "simba_null_method": "multi_shift",  # "single_shift" or "multi_shift"
    # Multi-shift pool chosen from empirically safe offsets (~18-30s at 10 fps)
    # to avoid preserving trial-locked structure while keeping realistic dynamics.
    "simba_shift_pool_seconds": [18, 20, 23, 25, 28, 30],
    "simba_null_random_state": 42,
    "simba_n_shuffles": 100,
    "simba_ci_percentile": 95,
    # Exclude edge frames when estimating session mean/std for z-scoring.
    # This helps avoid start/end artifacts dominating normalization.
    "simba_zscore_edge_exclude_frames": 1800,

    # ── Photometry parameters ──
    "photo_baseline_seconds": 5,
    "photo_triallength_seconds": 20,
    "photo_bins": 200,

    # ── Conditions to exclude ──
    "conditions_to_exclude": ["thirsty", "replete_exp"],

    # ── Snip parameters for AUC calculation ──
    "auc_start_bin": 50,   # bin index for start of infusion window
    "auc_end_bin": 150,    # bin index for end of infusion window
    "auc_early_start_bin": 50,
    "auc_early_end_bin": 80,
    "auc_late_start_bin": 80,
    "auc_late_end_bin": 150,

    # ── Velocity smoothing ──
    "vel_smooth_window": 5,

    # ── Clustering ──
    "num_retained_pcs": 3,
    "n_clusters": None,           # None = auto-select via silhouette sweep; or set to int to force
    "max_n_clusters": 9,          # Maximum number of clusters to test in sweep
    "clustering_affinity": "sigmoid",
    "clustering_assign_labels": "discretize",

    # ── Movement analysis parameters ──
    "normalize_movement": False,              # Normalize to [0,1] by default
    "movement_threshold": 3,              # Threshold for normalized movement
    "calculate_raw_movement": True,         # Optional: also analyze raw pixels
    "movement_threshold_raw": 0.5,           # Threshold for raw pixel movement (pixels/frame)

    # ── Angular velocity analysis parameters ──
    "angvel_threshold": 1.0,                 # Threshold for angular velocity (deg/frame)

    # ── Sigmoidal transition fit ──
    "transition_condition": "deplete",
    "transition_infusion": "45NaCl",
    "logistic_maxfev": 60000,
    # Behavioral transition mode for realignment:
    #   - "median_balance": continuous sigmoid on simba_median_balance
    #   - "binarized": logistic fit on thresholded simba_median_balance (middle values ignored)
    #   - "quantized": continuous sigmoid on quantized (-1/0/+1) simba_median_balance
    "behavior_transition_mode": "median_balance",
    "behavior_transition_signal_col": "simba_median_balance",
    "behavior_transition_low_threshold": -0.7,
    "behavior_transition_high_threshold": 0.7,
    "behavior_transition_maxfev": 30000,
    "behavior_transition_require_negative_k": True,
    "behavior_transition_max_x0": 50,

    # ── Behavioral transition method ──
    # Selects how per-rat behavioral transition trials are determined:
    #   - "sigmoid": fit a continuous sigmoid to the signal (existing default)
    #   - "max_velocity": use the trial of maximum negative derivative (steepest decline)
    "behavior_transition_method": "sigmoid",
    # Smoothing window (must be odd) for the derivative when using "max_velocity".
    "behavior_velocity_smoothing_window": 3,

    # ── Output ──
    "output_filename": "assembled_data.pickle",

    # ── Caching ──
    # Set these to True to load from cached pickle instead of re-extracting.
    # Cache files are saved automatically after each step.
    # If the cache file doesn't exist, the step runs from scratch regardless.
    "cache_behav": True,          # Skip DLC extraction, load from cache
    "cache_photo": True,          # Skip TDT extraction, load from cache
    "cache_simba": True,          # Skip Simba extraction, load from cache
    "cache_clustering": True,     # Skip PCA + spectral clustering, load from cache
    "cache_transitions": True,    # Skip sigmoidal fitting, load from cache

    # Cache filenames (in data_folder)
    "cache_behav_file": "_cache_behav.pickle",
    "cache_photo_file": "_cache_photo.pickle",
    "cache_simba_file": "_cache_simba.pickle",
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
        normalize=False,
        calibration_factor=1.0,  # No scaling applied to movement metric
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
    df = calc_angular_velocity(df, rightear="r_ear", leftear="l_ear", absolute=True)
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
        baseline_length=params["photo_baseline_seconds"],
        trial_length=params["photo_triallength_seconds"],
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
    """Step 1: Assemble behavioural (DLC) snips and x_array. Always calculates both movement and angular velocity."""
    print("\n" + "=" * 60)
    print("STEP 1: Assembling behavioural data (movement + angular velocity)")
    print("=" * 60)

    data_folder = params["data_folder"]
    dlc_folder = params["dlc_folder"]

    meta_df = pd.concat([
        pd.read_csv(data_folder / "10NaCl_FileKey.csv"),
        pd.read_csv(data_folder / "45NaCl_FileKey.csv"),
    ])

    snips_movement_list, snips_angvel_list, snips_movement_raw_list, x_list = [], [], [], []
    behav_stats = {}  # Store mean/std for each file's behavior vectors
    
    for _, row in meta_df.iterrows():
        stub = row["Folder"]
        print(f"  Behaviour: {stub}", end=" ... ")
        try:
            # Get movement data
            movement_result = get_movement_vector(stub, dlc_folder, params)
            if isinstance(movement_result, tuple):
                movement_vec, movement_vec_raw = movement_result
            else:
                movement_vec = movement_result
                movement_vec_raw = None
            
            # Get angular velocity data
            angvel_vec = get_angular_velocity_vector(stub, dlc_folder, params)
            
            # Calculate mean and std for behavior vectors, excluding first and last 1800 frames
            exclude_frames = 1800
            valid_range = slice(exclude_frames, -exclude_frames if exclude_frames > 0 else None)
            
            movement_mean = np.nanmean(movement_vec[valid_range])
            movement_std = np.nanstd(movement_vec[valid_range])
            angvel_mean = np.nanmean(angvel_vec[valid_range])
            angvel_std = np.nanstd(angvel_vec[valid_range])
            
            # Store with file identifier (for reference/metadata)
            file_id = f"{row['Subject']}_{stub}"
            behav_stats[file_id] = {
                "subject": row["Subject"],
                "folder": stub,
                "treatment": row["TreatNum"],
                "condition": row["Physiological state"],
                "movement_mean": movement_mean,
                "movement_std": movement_std,
                "angvel_mean": angvel_mean,
                "angvel_std": angvel_std,
            }
            
            # Create snips from both
            solenoid_ts = get_ttls(stub, data_folder)
            snips_movement = get_behav_snips(
                solenoid_ts=solenoid_ts,
                behav_vector=movement_vec,
                zscore_to_baseline=params["dlc_zscore_to_baseline"],
            )
            snips_angvel = get_behav_snips(
                solenoid_ts=solenoid_ts,
                behav_vector=angvel_vec,
                zscore_to_baseline=False,  # Don't zscore angular velocity
            )
            
            n = len(snips_movement)
            print(f"{n} trials")
            snips_movement_list.append(snips_movement)
            snips_angvel_list.append(snips_angvel)
            
            # Also create snips for raw movement if available
            if movement_vec_raw is not None:
                snips_movement_raw = get_behav_snips(
                    solenoid_ts=solenoid_ts,
                    behav_vector=movement_vec_raw,
                    zscore_to_baseline=False,  # Don't zscore raw pixels
                )
                snips_movement_raw_list.append(snips_movement_raw)

            infusion_label = "45NaCl" if row["TreatNum"] == 45 else "10NaCl"
            x_list.append(pd.DataFrame({
                "trial": np.arange(n),
                "id": row["Subject"],
                "condition": row["Physiological state"],
                "infusiontype": infusion_label,
                "movement_mean": movement_mean,
                "movement_std": movement_std,
                "angvel_mean": angvel_mean,
                "angvel_std": angvel_std,
            }))
        except Exception as e:
            print(f"ERROR: {e}")

    snips_movement_all = np.concatenate(snips_movement_list)
    snips_angvel_all = np.concatenate(snips_angvel_list)
    snips_movement_raw_all = np.concatenate(snips_movement_raw_list) if snips_movement_raw_list else None
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
    snips_movement_all = snips_movement_all[mask.values]
    snips_angvel_all = snips_angvel_all[mask.values]
    if snips_movement_raw_all is not None:
        snips_movement_raw_all = snips_movement_raw_all[mask.values]
    x_all = x_all[mask].reset_index(drop=True)

    print(f"  Behaviour assembled: {snips_movement_all.shape[0]} trials, {snips_movement_all.shape[1]} bins")
    print(f"    Movement snips: {snips_movement_all.shape}")
    print(f"    Angular velocity snips: {snips_angvel_all.shape}")
    if snips_movement_raw_all is not None:
        print(f"  Raw movement snips also assembled")
    return snips_movement_all, snips_angvel_all, snips_movement_raw_all, x_all, behav_stats


def assemble_simba(params):
    """Assemble Simba prediction snips and x_array metadata."""
    print("\n" + "=" * 60)
    print("STEP 2b: Assembling Simba prediction snips")
    print("=" * 60)

    data_folder = params["data_folder"]
    simba_folder = params["simba_folder"]

    meta_df = pd.concat([
        pd.read_csv(data_folder / "10NaCl_FileKey.csv"),
        pd.read_csv(data_folder / "45NaCl_FileKey.csv"),
    ])

    snips_raw_list, snips_baseline_list, snips_zscore_list, x_list = [], [], [], []

    def _resolve_shift_pool_frames(local_params):
        fps_local = float(local_params["simba_fps"])
        seconds_pool = np.asarray(local_params.get("simba_shift_pool_seconds", []), dtype=float)
        seconds_pool = seconds_pool[np.isfinite(seconds_pool)]
        frames_pool = np.unique(np.round(seconds_pool * fps_local).astype(int))
        frames_pool = frames_pool[frames_pool != 0]
        if frames_pool.size == 0:
            raise ValueError("simba_shift_pool_seconds must provide at least one non-zero shift.")
        return frames_pool

    for _, row in meta_df.iterrows():
        stub = row["Folder"]
        print(f"  Simba: {stub}", end=" ... ")
        try:
            prob_vec, src_file = read_simba_probability(
                stub,
                simba_folder,
                probability_column=params["simba_probability_column"],
            )

            solenoid_ts = get_ttls(stub, data_folder)
            cam1_onsets = None
            alignment_method = "fps"
            if params.get("simba_use_cam1_timestamps", True):
                try:
                    cam1_onsets = get_cam1_onsets_for_stub(stub, params["tank_folder"])
                    alignment_method = "cam1"
                except Exception as cam_exc:
                    print(f"Cam1 unavailable, falling back to fps ({cam_exc})", end="; ")

            snips = make_simba_snips(
                prob_vec,
                solenoid_ts,
                fps=params["simba_fps"],
                pre_bins=params["simba_pre_bins"],
                post_bins=params["simba_post_bins"],
                cam1_onsets=cam1_onsets,
            )

            snips_raw = np.asarray(snips, dtype=float)

            null_method = str(params.get("simba_null_method", "single_shift")).lower()
            if null_method == "multi_shift":
                shift_pool_frames = _resolve_shift_pool_frames(params)
                shifted_means = get_shifted_snip_means_multi_shift(
                    prob_vec,
                    solenoid_ts,
                    fps=params["simba_fps"],
                    pre_bins=params["simba_pre_bins"],
                    post_bins=params["simba_post_bins"],
                    shift_frames_pool=shift_pool_frames,
                    n_shuffles=params["simba_n_shuffles"],
                    cam1_onsets=cam1_onsets,
                    random_state=params.get("simba_null_random_state", None),
                )
            else:
                shifted_means = get_shifted_snip_means(
                    prob_vec,
                    solenoid_ts,
                    fps=params["simba_fps"],
                    pre_bins=params["simba_pre_bins"],
                    post_bins=params["simba_post_bins"],
                    shift_frames=params["simba_shift_frames"],
                    n_shuffles=params["simba_n_shuffles"],
                    cam1_onsets=cam1_onsets,
                )

            _, simba_ci_upper = get_time_above_simba_ci(
                snips,
                shifted_means,
                percentile=params["simba_ci_percentile"],
                start_bin=params["auc_start_bin"],
                end_bin=params["auc_end_bin"],
            )

            # Median-balance score: centered around zero by construction.
            # score = 2 * p(above null median) - 1, in [-1, 1].
            simba_pct_time_above_50ci, _ = get_time_above_simba_ci(
                snips,
                shifted_means,
                percentile=50,
                start_bin=params["auc_start_bin"],
                end_bin=params["auc_end_bin"],
            )
            simba_median_balance = (2.0 * simba_pct_time_above_50ci) - 1.0

            # Session-level z-scoring for Simba snips, excluding start/end frames
            # to reduce edge artifacts in the normalization reference.
            edge_exclude = int(params.get("simba_zscore_edge_exclude_frames", 0))
            if edge_exclude > 0 and prob_vec.size > (2 * edge_exclude):
                prob_reference = np.asarray(prob_vec, dtype=float)[edge_exclude:-edge_exclude]
            else:
                prob_reference = np.asarray(prob_vec, dtype=float)

            prob_mean = float(np.nanmean(prob_reference))
            prob_std = float(np.nanstd(prob_reference))

            # ── Variant 1: null-median baseline subtracted (units: Δprobability) ──
            # Subtract the grand mean of all shuffled null snips per bin.
            null_grand_mean = float(np.nanmean(shifted_means))
            snips_baseline = snips_raw - null_grand_mean

            # ── Variant 2: z-scored to entire session file ──
            # Standardise using the full-session probability vector mean and SD.
            if prob_std > 0 and not np.isnan(prob_std):
                snips_zs = (snips_raw - prob_mean) / prob_std
            else:
                snips_zs = np.full_like(snips_raw, np.nan, dtype=float)

            simba_raw_mean = np.nanmean(
                snips_raw[:, params["auc_start_bin"]:params["auc_end_bin"]],
                axis=1,
            )

            # Mean of z-scored snips in infusion window (5-15 s by default).
            simba_zscore_mean = np.nanmean(
                snips_zs[:, params["auc_start_bin"]:params["auc_end_bin"]],
                axis=1,
            )

            n = len(snips)
            print(
                f"{n} trials ({src_file.name}, upper CI={simba_ci_upper:.3f}, "
                f"align={alignment_method}, null={null_method})"
            )
            snips_raw_list.append(snips_raw)
            snips_baseline_list.append(snips_baseline)
            snips_zscore_list.append(snips_zs)

            infusion_label = "45NaCl" if row["TreatNum"] == 45 else "10NaCl"
            x_list.append(pd.DataFrame({
                "trial": np.arange(n),
                "id": row["Subject"],
                "condition": row["Physiological state"],
                "infusiontype": infusion_label,
                "simba_raw_mean": simba_raw_mean,
                "simba_zscore_mean": simba_zscore_mean,
                "simba_median_balance": simba_median_balance,
            }))
        except Exception as e:
            print(f"ERROR: {e}")

    if not snips_baseline_list:
        raise RuntimeError("No Simba snips were assembled. Check simba_folder and input files.")

    snips_raw_all = np.concatenate(snips_raw_list)
    snips_baseline_all = np.concatenate(snips_baseline_list)
    snips_zscore_all = np.concatenate(snips_zscore_list)
    x_all = (
        pd.concat(x_list, ignore_index=True)
        .replace({"condition": _condition_map()})
    )

    subjects_df = (
        pd.concat([
            pd.read_csv(data_folder / "10NaCl_SubjectKey.csv").iloc[:, :2],
            pd.read_csv(data_folder / "45NaCl_SubjectKey.csv").iloc[:, :2],
        ])
        .drop_duplicates()
        .rename(columns={"Subject": "id", "Sex": "sex"})
    )
    x_all = pd.merge(x_all, subjects_df, on="id", how="left")

    mask = ~x_all.condition.isin(params["conditions_to_exclude"])
    snips_raw_all = snips_raw_all[mask.values]
    snips_baseline_all = snips_baseline_all[mask.values]
    snips_zscore_all = snips_zscore_all[mask.values]
    x_all = x_all[mask].reset_index(drop=True)

    print(f"  Simba assembled: {snips_baseline_all.shape[0]} trials, {snips_baseline_all.shape[1]} bins")
    print("    snips_simba_raw: uncorrected raw Simba probabilities")
    print("    snips_simba_shifted_baseline: null grand-mean subtracted (delta probability)")
    print("    snips_simba_zscore: session z-score (edge-trimmed reference)")
    return snips_raw_all, snips_baseline_all, snips_zscore_all, x_all


# ──────────────────────────────────────────────────────────────────────
# STEP 3: Equalize and combine
# ──────────────────────────────────────────────────────────────────────

def equalize_datasets(
    x_photo,
    snips_photo,
    x_behav,
    snips_movement,
    snips_angvel,
    snips_movement_raw=None,
    x_simba=None,
    snips_simba_raw=None,
    snips_simba=None,
    snips_simba_zscore=None,
):
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

    merged = pd.merge(
        df_p.assign(_idx_p=df_p.index),
        df_b.assign(_idx_b=df_b.index),
        on=merge_cols,
        how="inner",
    )

    if x_simba is not None and snips_simba is not None:
        df_s = x_simba[merge_cols].reset_index(drop=True)
        merged = pd.merge(
            merged,
            df_s.assign(_idx_s=df_s.index),
            on=merge_cols,
            how="inner",
        )

    idx_p = merged["_idx_p"].values
    idx_b = merged["_idx_b"].values

    x_photo = x_photo.iloc[idx_p].reset_index(drop=True)
    snips_photo = snips_photo[idx_p]
    x_behav = x_behav.iloc[idx_b].reset_index(drop=True)
    snips_movement = snips_movement[idx_b]
    snips_angvel = snips_angvel[idx_b]

    if snips_movement_raw is not None:
        snips_movement_raw = snips_movement_raw[idx_b]

    if x_simba is not None and snips_simba is not None:
        idx_s = merged["_idx_s"].values
        x_simba = x_simba.iloc[idx_s].reset_index(drop=True)
        if snips_simba_raw is not None:
            snips_simba_raw = snips_simba_raw[idx_s]
        snips_simba = snips_simba[idx_s]
        if snips_simba_zscore is not None:
            snips_simba_zscore = snips_simba_zscore[idx_s]

    print(f"  After equalization: {len(x_photo)} trials (photo), {len(x_behav)} trials (behav)")
    if x_simba is not None:
        print(f"                     {len(x_simba)} trials (simba)")
    assert len(x_photo) == len(x_behav), "Datasets still not aligned!"
    if x_simba is not None and snips_simba is not None:
        assert len(x_photo) == len(x_simba), "Simba dataset is not aligned!"
    return x_photo, snips_photo, x_behav, snips_movement, snips_angvel, snips_movement_raw, x_simba, snips_simba_raw, snips_simba, snips_simba_zscore


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
                random_state=123,
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
        # n_clusters=n_clusters,
        n_clusters=2, # Force 2 clusters for better interpretability, silhouette is almost identical to n=4
        affinity=params["clustering_affinity"],
        assign_labels=params["clustering_assign_labels"],
        random_state=123,
    )
    model.fit(transformed[:, :num_pcs])
    sil = silhouette_score(transformed[:, :num_pcs], model.labels_, metric="cosine")
    print(f"  Final clustering: k={n_clusters}, silhouette={sil:.3f}")

    # Reorder clusters so cluster 0 = most positive response during infusion
    # This matches the reorder_clusters function from spectral_clustering_all_trials.ipynb
    pre_window = int(params["photo_baseline_seconds"] * 10)  # 50 bins
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
    fits_df, _ = fit_logistic_transitions_by_id(
        df_dep_45,
        signal_col="cluster_photo",
        value_transform=lambda sig: np.logical_not(sig).astype(int),
        direction="decreasing",
        maxfev=params["logistic_maxfev"],
        min_x0=0,
    )

    print(f"  Successful fits: {len(fits_df)} / {len(df_dep_45.id.unique())} rats")
    print(f"  Transition points: {fits_df.x0_orig.round(1).tolist()}")
    return fits_df


def find_behavioral_transitions(x_array, params, signal_col="simba_median_balance"):
    """Fit per-rat behavioral transitions for deplete + 45NaCl.

    Controlled by params["behavior_transition_mode"]:
      - "median_balance": continuous sigmoid on signal_col
      - "binarized": logistic on thresholded signal_col
      - "quantized": continuous sigmoid on {-1,0,+1} transformed signal_col
    """
    mode = params.get("behavior_transition_mode", "median_balance")
    signal_col = params.get("behavior_transition_signal_col", signal_col)
    low_thr = float(params.get("behavior_transition_low_threshold", -0.7))
    high_thr = float(params.get("behavior_transition_high_threshold", 0.7))
    maxfev = int(params.get("behavior_transition_maxfev", 30000))
    require_negative_k = bool(params.get("behavior_transition_require_negative_k", True))
    max_x0 = params.get("behavior_transition_max_x0", 50)

    print(f"  Calculating behavioral transitions using mode='{mode}' from {signal_col}")

    if signal_col not in x_array.columns:
        print(f"  Behavioral transition column '{signal_col}' not found; skipping behavior fits.")
        return pd.DataFrame(columns=["id", "x0_orig", "success", "note"])

    cond = params["transition_condition"]
    inf = params["transition_infusion"]
    df_dep_45 = x_array.query("condition == @cond & infusiontype == @inf").copy()
    bin_col = f"{signal_col}_binarized"
    quant_col = f"{signal_col}_quantized"

    if mode == "median_balance":
        fits_df_behav, _ = fit_continuous_sigmoid_transitions_by_id(
            df_dep_45,
            signal_col=signal_col,
            maxfev=maxfev,
            require_negative_k=require_negative_k,
            min_x0=0,
            max_x0=max_x0,
        )
    elif mode == "binarized":
        if bin_col in df_dep_45.columns:
            fits_df_behav, _ = fit_logistic_transitions_by_id(
                df_dep_45,
                signal_col=bin_col,
                value_transform=None,
                direction="decreasing",
                maxfev=max(int(params.get("logistic_maxfev", 60000)), maxfev),
                min_x0=0,
                max_x0=max_x0,
            )
        else:
            fits_df_behav, _ = fit_logistic_transitions_by_id(
                df_dep_45,
                signal_col=signal_col,
                value_transform=lambda values: binarize_series_threshold(values, low=low_thr, high=high_thr),
                direction="decreasing",
                maxfev=max(int(params.get("logistic_maxfev", 60000)), maxfev),
                min_x0=0,
                max_x0=max_x0,
            )
    elif mode == "quantized":
        if quant_col not in df_dep_45.columns:
            df_dep_45 = df_dep_45.assign(
                **{quant_col: quantize_series_threshold(df_dep_45[signal_col].to_numpy(dtype=float), low=low_thr, high=high_thr)}
            )
        fits_df_behav, _ = fit_continuous_sigmoid_transitions_by_id(
            df_dep_45,
            signal_col=quant_col,
            maxfev=maxfev,
            require_negative_k=require_negative_k,
            min_x0=0,
            max_x0=max_x0,
        )
    else:
        raise ValueError(
            f"Unknown behavior_transition_mode '{mode}'. "
            "Expected one of: median_balance, binarized, quantized"
        )

    if not fits_df_behav.empty:
        fits_df_behav = fits_df_behav.copy()
        fits_df_behav["behavior_transition_mode"] = mode

    print(f"  Successful behavioral fits: {len(fits_df_behav)} / {len(df_dep_45.id.unique())} rats")
    if not fits_df_behav.empty:
        print(f"  Behavioral transition points: {fits_df_behav.x0_orig.round(1).tolist()}")
    return fits_df_behav


def find_velocity_transitions(x_array, params, signal_col=None):
    """Find per-rat behavioral transition using maximum negative velocity (derivative).

    Computes the smoothed first derivative of *signal_col* for each rat in the
    deplete + 45NaCl subset and returns the trial at which the most negative
    velocity (steepest behavioural decline) occurs.

    Parameters
    ----------
    x_array : pd.DataFrame
        Full assembled trial-level DataFrame.
    params : dict
        Pipeline parameters.  Reads:
          - ``behavior_transition_signal_col`` (default ``"simba_median_balance"``)
          - ``behavior_velocity_smoothing_window`` (default ``3``, must be odd)
          - ``transition_condition`` / ``transition_infusion`` for subsetting.
    signal_col : str or None
        Override the signal column; falls back to *params* if None.

    Returns
    -------
    pd.DataFrame
        Columns: ``id``, ``x0_orig``, ``velocity_value``, ``n_trials``, ``success``.
        Compatible with :func:`model_fit_helpers.build_realigned_trials`.
    """
    if signal_col is None:
        signal_col = params.get("behavior_transition_signal_col", "simba_median_balance")
    window_size = int(params.get("behavior_velocity_smoothing_window", 3))
    if window_size % 2 == 0:
        window_size += 1  # enforce odd for symmetric smoothing

    cond = params["transition_condition"]
    inf = params["transition_infusion"]
    df_dep_45 = x_array.query("condition == @cond & infusiontype == @inf").copy()

    print(f"  Finding velocity-based transitions for {len(df_dep_45.id.unique())} rats")
    print(f"  Signal: {signal_col}, smoothing window: {window_size}")

    if signal_col not in df_dep_45.columns:
        print(f"  Column '{signal_col}' not found; returning empty transitions.")
        return pd.DataFrame(columns=["id", "x0_orig", "velocity_value", "n_trials", "success"])

    results = []
    for rat in df_dep_45.id.unique():
        rat_data = df_dep_45.loc[df_dep_45.id == rat].copy().sort_values("trial")
        signal = rat_data[signal_col].to_numpy(dtype=float)

        if len(signal) < 2 or np.all(np.isnan(signal)):
            results.append({
                "id": rat, "x0_orig": np.nan, "velocity_value": np.nan,
                "n_trials": len(signal), "success": False,
            })
            continue

        derivative = np.diff(signal, prepend=signal[0])
        smoothed = uniform_filter1d(derivative, size=window_size, mode="nearest")
        max_vel_idx = int(np.argmin(smoothed))  # most negative velocity
        max_vel = float(smoothed[max_vel_idx])
        transition_trial = int(rat_data.iloc[max_vel_idx]["trial"])

        results.append({
            "id": rat,
            "x0_orig": float(transition_trial),
            "velocity_value": max_vel,
            "n_trials": len(signal),
            "success": True,
        })

    fits_df = pd.DataFrame(results)
    success_count = int(fits_df["success"].sum()) if not fits_df.empty else 0
    print(f"  Velocity transitions: {success_count}/{len(fits_df)} rats")
    if not fits_df.empty:
        print(f"  Transition trials: {fits_df.loc[fits_df['success'], 'x0_orig'].round(1).tolist()}")
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


def get_time_above_angvel_threshold(snips, threshold=1.0, start_bin=50, end_bin=150):
    """Calculate proportion of bins above angular velocity threshold per trial."""
    angvel_above = []
    for i in range(snips.shape[0]):
        snip = snips[i, start_bin:end_bin]
        tmp = len([x for x in snip if x > threshold]) / len(snip)
        angvel_above.append(tmp)
    return np.array(angvel_above)


def get_auc_by_window(snips, start_bin=50, end_bin=150):
    """Calculate per-trial AUC using trapezoidal integration over a bin window."""
    snips = np.asarray(snips, dtype=float)
    return np.array([np.trapezoid(snips[i, start_bin:end_bin])/10 for i in range(len(snips))])


def get_mean_by_window(snips, start_bin=50, end_bin=150):
    """Calculate per-trial mean over a bin window."""
    snips = np.asarray(snips, dtype=float)
    return np.nanmean(snips[:, start_bin:end_bin], axis=1)


def sync_aligned_columns(target_df, source_df, merge_cols=None):
    """Copy aligned metadata columns from source_df into target_df.

    This keeps freshly recomputed columns available even when downstream cached
    dataframes were created before those columns existed.
    """
    if merge_cols is None:
        merge_cols = ["trial", "id", "condition", "infusiontype"]

    source_cols = [col for col in source_df.columns if col not in merge_cols]
    if not source_cols:
        return target_df

    if len(target_df) == len(source_df) and target_df[merge_cols].equals(source_df[merge_cols]):
        for col in source_cols:
            target_df[col] = source_df[col].values
        return target_df

    aligned_source = (
        source_df[merge_cols + source_cols]
        .drop_duplicates(subset=merge_cols)
    )
    target_df = (
        target_df
        .assign(_row_order=np.arange(len(target_df)))
        .merge(aligned_source, on=merge_cols, how="left", suffixes=("", "_fresh"))
        .sort_values("_row_order")
        .drop(columns=["_row_order"])
        .reset_index(drop=True)
    )

    for col in source_cols:
        fresh_col = f"{col}_fresh"
        if fresh_col in target_df.columns:
            target_df[col] = target_df[fresh_col]
            target_df = target_df.drop(columns=[fresh_col])

    return target_df


def _drop_legacy_simba_columns(df):
    """Remove deprecated Simba columns so stale cache columns do not persist."""
    legacy_cols = [
        "mean_simba",
        "simba_pct_time_above_95ci",
        "simba_median_balance_score",
        "simba_median_balance_score_pct",
        "simba_mean_zscore_auc",
        "simba_primary_metric_value",
        "simba_primary_metric_name",
        "simba_trial_metric_value",
        "simba_trial_metric_name",
        "simba_null_method",
        "simba_alignment_method",
    ]
    cols_to_drop = [c for c in legacy_cols if c in df.columns]
    if cols_to_drop:
        df = df.drop(columns=cols_to_drop)
        print(f"  Dropped legacy Simba columns: {cols_to_drop}")
    return df


def _add_behavior_transition_columns(df, params):
    """Add binarized/quantized behavior columns so transition modes can reuse them."""
    signal_col = params.get("behavior_transition_signal_col", "simba_median_balance")
    low_thr = float(params.get("behavior_transition_low_threshold", -0.7))
    high_thr = float(params.get("behavior_transition_high_threshold", 0.7))

    if signal_col not in df.columns:
        print(f"  Could not create transition helper columns: '{signal_col}' not in x_array")
        return df

    bin_col = f"{signal_col}_binarized"
    quant_col = f"{signal_col}_quantized"

    signal = df[signal_col].to_numpy(dtype=float)
    df[bin_col] = binarize_series_threshold(signal, low=low_thr, high=high_thr)
    df[quant_col] = quantize_series_threshold(signal, low=low_thr, high=high_thr)
    print(f"  Added transition helper columns: {bin_col}, {quant_col}")
    return df


def combine_and_realign(x_photo, snips_photo, snips_movement, snips_angvel, fits_df_da, fits_df_behav, params, snips_movement_raw=None):
    """
    Step 7: Add AUCs and time_moving to x_array, create realigned deplete+45NaCl subset.
    """
    print("\n" + "=" * 60)
    print("STEP 7: Combining data and realigning trials")
    print("=" * 60)

    # Smooth velocity snips
    snips_movement_smooth = smooth_array(snips_movement, window_size=params["vel_smooth_window"])
    snips_angvel_smooth = smooth_array(snips_angvel, window_size=params["vel_smooth_window"])

    # Calculate AUCs using trapezoidal rule (true area under curve)
    s, e = params["auc_start_bin"], params["auc_end_bin"]
    s_early = params.get("auc_early_start_bin", 50)
    e_early = params.get("auc_early_end_bin", 80)
    s_late = params.get("auc_late_start_bin", 80)
    e_late = params.get("auc_late_end_bin", 150)
    auc_snips = get_auc_by_window(snips_photo, start_bin=s, end_bin=e)
    auc_snips_early = get_auc_by_window(snips_photo, start_bin=s_early, end_bin=e_early)
    auc_snips_late = get_auc_by_window(snips_photo, start_bin=s_late, end_bin=e_late)
    auc_movement = get_auc_by_window(snips_movement_smooth, start_bin=s, end_bin=e)
    auc_angvel = get_auc_by_window(snips_angvel_smooth, start_bin=s, end_bin=e)

    # Calculate time moving (normalized)
    time_moving = get_time_moving(snips_movement, threshold=params["movement_threshold"],
                                   start_bin=s, end_bin=e)

    # Calculate time moving raw (if available)
    time_moving_raw = None
    if snips_movement_raw is not None:
        time_moving_raw = get_time_moving(snips_movement_raw, threshold=params["movement_threshold_raw"],
                                          start_bin=s, end_bin=e)

    # Calculate proportion of time above angular velocity threshold
    time_above_angvel_threshold = get_time_above_angvel_threshold(
        snips_angvel, threshold=params["angvel_threshold"],
        start_bin=s, end_bin=e)

    x_combined = x_photo.assign(
        auc_snips=auc_snips,
        auc_snips_early=auc_snips_early,
        auc_snips_late=auc_snips_late,
        auc_movement=auc_movement,
        auc_angvel=auc_angvel,
        time_moving=time_moving,
        time_above_angvel_threshold=time_above_angvel_threshold,
    )
    
    # Add time_moving_raw if available
    if time_moving_raw is not None:
        x_combined = x_combined.assign(time_moving_raw=time_moving_raw)

    realigned_trials_da = build_realigned_trials(x_combined, fits_df_da, "realigned_trials_da")
    realigned_trials_behav = build_realigned_trials(x_combined, fits_df_behav, "realigned_trials_behav")

    x_combined = x_combined.assign(
        realigned_trials_da=realigned_trials_da,
        realigned_trials_behav=realigned_trials_behav,
        # Backward-compat alias while notebooks migrate to the clearer DA-specific name.
        trial_aligned=realigned_trials_da,
    )

    # Now realign the deplete+45NaCl subset (with NaN values dropped)
    cond = params["transition_condition"]
    inf = params["transition_infusion"]
    z = (
        x_combined
        .query("condition == @cond & infusiontype == @inf")
        .dropna(subset=['realigned_trials_da'])
        .reset_index(drop=True)
    )

    print(f"  Combined x_array: {len(x_combined)} trials")
    print(
        f"  Added realigned_trials_da column "
        f"(NaN for {(x_combined['realigned_trials_da'].isna()).sum()} trials without DA fits)"
    )
    print(
        f"  Added realigned_trials_behav column "
        f"(NaN for {(x_combined['realigned_trials_behav'].isna()).sum()} trials without behavioral fits)"
    )
    if time_moving_raw is not None:
        print(f"  Added time_moving_raw column (threshold={params['movement_threshold_raw']} pixels)")
    print(f"  Realigned deplete+45NaCl subset (DA-aligned): {len(z)} trials with valid alignments")

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
    print(f"  DLC folder: {params['dlc_folder']}")
    print(f"  Tank folder: {params['tank_folder']}")
    print(f"  Simba folder: {params['simba_folder']}")
    print(f"  Calculating: movement + angular velocity (both always calculated)")
    print(f"  Caching: behav={params['cache_behav']}, photo={params['cache_photo']}, simba={params['cache_simba']}, "
          f"clustering={params['cache_clustering']}, transitions={params['cache_transitions']}")

    data_folder = params["data_folder"]

    # Step 1: Behaviour
    behav_cache_path = data_folder / params["cache_behav_file"]
    if params["cache_behav"]:
        cached = _load_cache(behav_cache_path, "behaviour")
    else:
        cached = None
    if cached is not None:
        snips_movement = cached["snips_movement"]
        snips_angvel = cached["snips_angvel"]
        snips_movement_raw = cached.get("snips_movement_raw", None)
        x_behav = cached["x_behav"]
        behav_stats = cached.get("behav_stats", {})
        print(f"  Behaviour from cache: {snips_movement.shape[0]} trials (movement + angvel)")
    else:
        snips_movement, snips_angvel, snips_movement_raw, x_behav, behav_stats = assemble_behaviour(params)
        cache_dict = {"snips_movement": snips_movement, "snips_angvel": snips_angvel, "x_behav": x_behav, "behav_stats": behav_stats}
        if snips_movement_raw is not None:
            cache_dict["snips_movement_raw"] = snips_movement_raw
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

    # Step 2b: Simba
    simba_cache_path = data_folder / params["cache_simba_file"]
    if params["cache_simba"]:
        cached = _load_cache(simba_cache_path, "simba")
    else:
        cached = None
    if cached is not None:
        snips_simba_raw = cached.get("snips_simba_raw")
        snips_simba_shifted_baseline = cached.get("snips_simba_shifted_baseline", cached.get("snips_simba_baseline"))
        snips_simba_zscore = cached["snips_simba_zscore"]
        x_simba = cached["x_simba"]
        if snips_simba_raw is None or "simba_raw_mean" not in x_simba.columns:
            print("  Simba cache missing raw outputs; recomputing Simba step.")
            cached = None
        else:
            print(f"  Simba from cache: {snips_simba_shifted_baseline.shape[0]} trials")

    if cached is None:
        snips_simba_raw, snips_simba_shifted_baseline, snips_simba_zscore, x_simba = assemble_simba(params)
        _save_cache(
            {"snips_simba_raw": snips_simba_raw,
             "snips_simba_shifted_baseline": snips_simba_shifted_baseline,
             "snips_simba_baseline": snips_simba_shifted_baseline,
             "snips_simba_zscore": snips_simba_zscore,
             "x_simba": x_simba},
            simba_cache_path, "simba", params)

    # Step 3: Equalize
    x_photo, snips_photo, x_behav, snips_movement, snips_angvel, snips_movement_raw, x_simba, snips_simba_raw, snips_simba_shifted_baseline, snips_simba_zscore = equalize_datasets(
        x_photo,
        snips_photo,
        x_behav,
        snips_movement,
        snips_angvel,
        snips_movement_raw,
        x_simba=x_simba,
        snips_simba_raw=snips_simba_raw,
        snips_simba=snips_simba_shifted_baseline,
        snips_simba_zscore=snips_simba_zscore,
    )
    
    # Add behavioral stats columns from x_behav to x_photo (they're now aligned after equalization)
    behav_cols = ["movement_mean", "movement_std", "angvel_mean", "angvel_std"]
    for col in behav_cols:
        if col in x_behav.columns:
            x_photo[col] = x_behav[col]

    if x_simba is not None:
        x_photo = _drop_legacy_simba_columns(x_photo)
        simba_cols = [
            "simba_raw_mean",
            "simba_zscore_mean",
            "simba_median_balance",
        ]
        for col in simba_cols:
            if col in x_simba.columns:
                x_photo[col] = x_simba[col].values
            else:
                print(f"  Warning: Simba metadata missing '{col}'. Rerun with cache_simba=False to populate it.")

        # Always make thresholded behavior transition columns available on x_array.
        x_photo = _add_behavior_transition_columns(x_photo, params)

    # Step 4: Clustering
    clustering_cache_path = data_folder / params["cache_clustering_file"]
    if params["cache_clustering"]:
        cached = _load_cache(clustering_cache_path, "clustering")
    else:
        cached = None
    if cached is not None:
        x_combined, pca_transformed = cached["x_combined"], cached["pca_transformed"]
        x_combined = sync_aligned_columns(x_combined, x_photo)
        x_combined = _drop_legacy_simba_columns(x_combined)
        print(f"  Clustering from cache: {len(x_combined)} trials")
    else:
        x_combined, pca_transformed = cluster_photometry(snips_photo, x_photo, params)
        # Step 5 always runs with fresh clustering
        x_combined = compute_cluster_distances(x_combined, pca_transformed, params)
        x_combined = _drop_legacy_simba_columns(x_combined)
        _save_cache({"x_combined": x_combined, "pca_transformed": pca_transformed},
                    clustering_cache_path, "clustering", params)

    # Step 5: Cluster distances (only if loaded from cache, otherwise already done above)
    if cached is not None:
        x_combined = compute_cluster_distances(x_combined, pca_transformed, params)

    # Step 6: Sigmoidal transitions
    # Uses deterministic clustering (random_state=0) to generate consistent results
    print("\nSTEP 6: Calculate sigmoidal transitions")
    print("=" * 60)
    transitions_cache_path = data_folder / params["cache_transitions_file"]
    if params["cache_transitions"]:
        cached = _load_cache(transitions_cache_path, "transitions")
    else:
        cached = None
    cache_mode = None
    if cached is not None:
        cache_mode = cached.get("behavior_transition_mode")
        if cache_mode is None and "_cached_params" in cached:
            cache_mode = cached["_cached_params"].get("behavior_transition_mode")

    # Cache is valid only when both the sigmoid mode AND transition method match.
    cached_method = None
    if cached is not None:
        cached_method = cached.get("behavior_transition_method")
        if cached_method is None and "_cached_params" in cached:
            cached_method = cached["_cached_params"].get("behavior_transition_method")

    current_method = params.get("behavior_transition_method", "sigmoid")
    current_mode = params.get("behavior_transition_mode", "median_balance")
    cache_valid = (
        cached is not None
        and "fits_df" in cached
        and "fits_df_behav" in cached
        and cache_mode == current_mode
        and cached_method == current_method
    )

    if cache_valid:
        fits_df = cached["fits_df"]
        fits_df_da = cached.get("fits_df_da", fits_df)
        fits_df_behav = cached["fits_df_behav"]
        print(
            f"  Transitions from cache: {len(fits_df_da)} DA fits, {len(fits_df_behav)} behavioral fits "
            f"(method={current_method}, mode={cache_mode})"
        )
    else:
        if cached is not None:
            print(
                f"  Transition cache mismatch (cached method={cached_method}/{cache_mode}, "
                f"current={current_method}/{current_mode}); recomputing transitions."
            )
        print(f"  Calculating transitions from deterministic clustering (random_state=0)")
        fits_df_da = find_sigmoidal_transitions(x_combined, params)
        if current_method == "max_velocity":
            fits_df_behav = find_velocity_transitions(x_combined, params)
        else:
            fits_df_behav = find_behavioral_transitions(x_combined, params)
        fits_df = fits_df_da
        print(f"  Calculated {len(fits_df_da)} DA transition fits and {len(fits_df_behav)} behavioral transition fits")
        _save_cache(
            {
                "fits_df": fits_df_da,
                "fits_df_da": fits_df_da,
                "fits_df_behav": fits_df_behav,
                "behavior_transition_mode": current_mode,
                "behavior_transition_method": current_method,
                "behavior_transition_signal_col": params.get("behavior_transition_signal_col", "simba_median_balance"),
                "behavior_transition_thresholds": (
                    params.get("behavior_transition_low_threshold", -0.7),
                    params.get("behavior_transition_high_threshold", 0.7),
                ),
            },
            transitions_cache_path,
            "transitions",
            params,
        )

    if cached is not None and ("fits_df" not in locals()):
        fits_df = fits_df_da

    # Step 7: Combine and realign
    x_combined, z_dep45 = combine_and_realign(
        x_combined, snips_photo, snips_movement, snips_angvel, fits_df_da, fits_df_behav, params, snips_movement_raw
    )

    # Create metadata about data processing
    metadata = {
        "behav_smoothed": params["dlc_smooth_method"] is not None,
        "behav_smooth_method": params["dlc_smooth_method"],
        "behav_smooth_window": params["dlc_smooth_window"] if params["dlc_smooth_method"] is not None else None,
        "behav_zscored": params["dlc_zscore_to_baseline"],
        "behav_bodyparts": params.get("dlc_bodyparts"),
        "dlc_likelihood_threshold": params.get("dlc_likelihood_threshold", 0.6),
        "simba_probability_column": params.get("simba_probability_column", "Probability_Appetitive"),
        "simba_shifted_baseline": params.get("simba_use_shifted_baseline", True),
        "simba_shift_frames": params.get("simba_shift_frames", 300),
        "simba_n_shuffles": params.get("simba_n_shuffles", 1000),
        "simba_ci_percentile": params.get("simba_ci_percentile", 97.5),
        "behavior_transition_mode": params.get("behavior_transition_mode", "median_balance"),
        "behavior_transition_method": params.get("behavior_transition_method", "sigmoid"),
        "behavior_velocity_smoothing_window": params.get("behavior_velocity_smoothing_window", 3),
        "behavior_transition_signal_col": params.get("behavior_transition_signal_col", "simba_median_balance"),
        "behavior_transition_low_threshold": params.get("behavior_transition_low_threshold", -0.7),
        "behavior_transition_high_threshold": params.get("behavior_transition_high_threshold", 0.7),
        "photo_smoothed": False,  # Photometry is NOT smoothed during assembly
        "photo_zscored": True,  # Photometry is z-scored by trompy during processing
        "behav_metrics": "movement + angular_velocity (both always calculated)",
        # Movement analysis parameters
        "normalize_movement": params.get("normalize_movement", True),
        "movement_threshold": params.get("movement_threshold", 0.02),
        "auc_start_bin": params.get("auc_start_bin", 50),
        "auc_end_bin": params.get("auc_end_bin", 150),
        "auc_early_start_bin": params.get("auc_early_start_bin", 50),
        "auc_early_end_bin": params.get("auc_early_end_bin", 80),
        "auc_late_start_bin": params.get("auc_late_start_bin", 80),
        "auc_late_end_bin": params.get("auc_late_end_bin", 150),
        "calculate_raw_movement": params.get("calculate_raw_movement", False),
        "movement_threshold_raw": params.get("movement_threshold_raw", 0.5) if params.get("calculate_raw_movement", False) else None,
        # Angular velocity parameters
        "angvel_threshold": params.get("angvel_threshold", 1.0),
    }

    # Save output
    output = {
        "x_array": x_combined,
        "snips_photo": snips_photo,
        "snips_simba_raw": snips_simba_raw,
        "snips_simba_shifted_baseline": snips_simba_shifted_baseline,
        "snips_simba_zscore": snips_simba_zscore,
        "snips_simba_baseline": snips_simba_shifted_baseline,  # backward-compat alias
        "snips_simba": snips_simba_shifted_baseline,  # backward-compat alias
        "snips_movement": snips_movement,
        "snips_angvel": snips_angvel,
        "pca_transformed": pca_transformed,
        "fits_df": fits_df,
        "fits_df_da": fits_df_da,
        "fits_df_behav": fits_df_behav,
        "z_dep45": z_dep45,
        "metadata": metadata,
        "params": params,
        "behav_stats": behav_stats,  # Mean/std for behavior vectors (movement & angvel) per file, excluding first/last 1800 frames
    }
    
    # Add optional raw movement data
    if snips_movement_raw is not None:
        output["snips_movement_raw"] = snips_movement_raw

    output_path = params["data_folder"] / params["output_filename"]
    with open(output_path, "wb") as f:
        dill.dump(output, f)

    print("\n" + "=" * 60)
    print(f"DONE! Saved to {output_path}")
    print("=" * 60)
    print(f"  x_array shape:          {x_combined.shape}")
    print(f"  snips_photo shape:      {snips_photo.shape}")
    print(f"  snips_simba_raw shape:      {snips_simba_raw.shape}")
    print(f"  snips_simba_shifted_baseline shape: {snips_simba_shifted_baseline.shape}")
    print(f"  snips_simba_zscore shape:   {snips_simba_zscore.shape}")
    print(f"  snips_movement shape:   {snips_movement.shape}")
    print(f"  snips_angvel shape:     {snips_angvel.shape}")
    print(f"  z_dep45 shape:     {z_dep45.shape}")
    print(f"  fits_df shape:     {fits_df.shape}")
    print(f"  fits_df_da shape:  {fits_df_da.shape}")
    print(f"  fits_df_behav shape: {fits_df_behav.shape}")
    print(f"  behavioral transition method: {metadata['behavior_transition_method']}")
    print(f"  behavioral transition mode: {metadata['behavior_transition_mode']}")
    print(f"\nData processing metadata:")
    print(f"  Behaviour metrics:    {metadata['behav_metrics']}")
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
