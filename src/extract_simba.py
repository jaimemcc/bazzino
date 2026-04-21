from pathlib import Path

import numpy as np
import pandas as pd
import tdt


def find_simba_file(stub, simba_folder):
    """Find a Simba predictions CSV that contains the session stub."""
    simba_folder = Path(simba_folder)
    matches = sorted(simba_folder.glob(f"*_{stub}_*.csv"))

    if not matches:
        raise FileNotFoundError(
            f"No Simba file found for stub '{stub}' in {simba_folder}"
        )

    # Prefer files tagged as Cam1 when multiple camera exports exist.
    cam1 = [p for p in matches if "Cam1" in p.name]
    return cam1[0] if cam1 else matches[0]


def read_simba_probability(stub, simba_folder, probability_column="Probability_Appetitive"):
    """Load the appetitive behavior probability vector for a session stub."""
    file_path = find_simba_file(stub, simba_folder)
    df = pd.read_csv(file_path)

    if probability_column not in df.columns:
        raise KeyError(
            f"Column '{probability_column}' not found in {file_path.name}. "
            f"Available columns include: {', '.join(df.columns[:10])}"
        )

    return pd.to_numeric(df[probability_column], errors="coerce").fillna(0.0).to_numpy(), file_path


def get_cam1_onsets_for_stub(stub, tank_folder):
    """Load Cam1 onset timestamps for a session stub from its TDT tank."""
    tank_path = Path(tank_folder) / stub
    data = tdt.read_block(str(tank_path), evtype=["epocs"])
    try:
        onsets = np.asarray(data.epocs["Cam1"].onset, dtype=float)
    except Exception as exc:
        available = []
        try:
            available = list(data.epocs.keys())
        except Exception:
            pass
        raise KeyError(
            f"Could not read Cam1 onsets for stub '{stub}'. "
            f"Available epocs: {available}"
        ) from exc

    if onsets.size < 2:
        raise ValueError(f"Cam1 onsets for stub '{stub}' have fewer than 2 timestamps.")

    return onsets


def _get_trial_ttls(solenoid_ts):
    """Return trial TTL timestamps, matching existing behavior of skipping final pulse."""
    ts = np.asarray(solenoid_ts, dtype=float)
    if ts.size < 2:
        return np.asarray([], dtype=float)
    ts = ts[:-1]
    return ts[np.isfinite(ts)]


def _map_tdt_times_to_nearest_frame_indices(event_times, cam1_onsets):
    """Map TDT event times (seconds) to nearest Cam1 frame index."""
    event_times = np.asarray(event_times, dtype=float)
    cam1_onsets = np.asarray(cam1_onsets, dtype=float)

    if cam1_onsets.size < 2:
        raise ValueError("cam1_onsets must have at least 2 timestamps.")

    idx_right = np.searchsorted(cam1_onsets, event_times, side="left")
    idx_right = np.clip(idx_right, 1, cam1_onsets.size - 1)
    idx_left = idx_right - 1

    dist_left = np.abs(event_times - cam1_onsets[idx_left])
    dist_right = np.abs(cam1_onsets[idx_right] - event_times)
    choose_right = dist_right < dist_left

    return np.where(choose_right, idx_right, idx_left).astype(int)


def _extract_single_snip(vector, center_frame, pre_bins, post_bins):
    """Extract one fixed-length snip and pad with NaN when near edges."""
    start = center_frame - pre_bins
    end = center_frame + post_bins

    out = np.full(pre_bins + post_bins, np.nan, dtype=float)
    src_start = max(start, 0)
    src_end = min(end, len(vector))

    if src_end <= src_start:
        return out

    dst_start = src_start - start
    dst_end = dst_start + (src_end - src_start)
    out[dst_start:dst_end] = vector[src_start:src_end]
    return out


def make_simba_snips(
    probability_vector,
    solenoid_ts,
    fps=10,
    pre_bins=50,
    post_bins=150,
    cam1_onsets=None,
):
    """Create TTL-aligned Simba snips, matching behavior snip convention."""
    trial_ttls = _get_trial_ttls(solenoid_ts)

    if cam1_onsets is None:
        center_frames = (trial_ttls * fps).astype(int)
    else:
        center_frames = _map_tdt_times_to_nearest_frame_indices(trial_ttls, cam1_onsets)

    # Clip frame centers to available Simba vector range.
    max_idx = max(len(probability_vector) - 1, 0)
    center_frames = np.clip(center_frames, 0, max_idx)

    snips = []
    for frame_on in center_frames:
        snips.append(
            _extract_single_snip(
                probability_vector,
                center_frame=int(frame_on),
                pre_bins=pre_bins,
                post_bins=post_bins,
            )
        )

    if not snips:
        return np.empty((0, pre_bins + post_bins), dtype=float)

    return np.asarray(snips, dtype=float)


def get_shifted_snip_means(
    probability_vector,
    solenoid_ts,
    fps=10,
    pre_bins=50,
    post_bins=150,
    shift_frames=300,
    n_shuffles=1000,
    cam1_onsets=None,
):
    """Generate null-distribution snips by circularly shifting the vector.

    Returns all individual shuffled snips stacked, shape
    (n_shuffles * n_trials, n_bins), so that per-bin percentile thresholds
    reflect single-trial variability rather than the variance of trial means.
    """
    shifted = np.asarray(probability_vector, dtype=float).copy()
    all_snips = []

    for _ in range(n_shuffles):
        shifted = np.roll(shifted, shift_frames)
        snips = make_simba_snips(
            shifted,
            solenoid_ts,
            fps=fps,
            pre_bins=pre_bins,
            post_bins=post_bins,
            cam1_onsets=cam1_onsets,
        )
        all_snips.append(snips)

    return np.concatenate(all_snips, axis=0)


def get_shifted_snip_means_multi_shift(
    probability_vector,
    solenoid_ts,
    fps=10,
    pre_bins=50,
    post_bins=150,
    shift_frames_pool=(180, 200, 230, 250, 280, 300),
    n_shuffles=1000,
    cam1_onsets=None,
    random_state=None,
):
    """Generate null-distribution snips using a pool of circular shifts.

    This is more robust than a single fixed shift when event timing has
    structure. For each shuffle, one shift is sampled from shift_frames_pool,
    the full probability vector is circularly shifted by that amount, and
    TTL-aligned snips are extracted.

    Returns stacked individual shuffled snips with shape
    (n_shuffles * n_trials, n_bins).
    """
    shifts = np.asarray(shift_frames_pool, dtype=int).ravel()
    shifts = shifts[np.isfinite(shifts)]
    shifts = shifts[shifts != 0]
    if shifts.size == 0:
        raise ValueError("shift_frames_pool must contain at least one non-zero shift.")

    rng = np.random.default_rng(random_state)
    base = np.asarray(probability_vector, dtype=float)
    all_snips = []

    for shift in rng.choice(shifts, size=int(n_shuffles), replace=True):
        shifted = np.roll(base, int(shift))
        snips = make_simba_snips(
            shifted,
            solenoid_ts,
            fps=fps,
            pre_bins=pre_bins,
            post_bins=post_bins,
            cam1_onsets=cam1_onsets,
        )
        all_snips.append(snips)

    return np.concatenate(all_snips, axis=0)


def baseline_simba_snips(real_snips, shifted_snips):
    """Subtract global shuffle mean from each real snip.

    Works with shifted_snips of any shape — the grand mean is the same
    whether the input is per-shuffle means (n_shuffles, n_bins) or all
    individual shuffled snips stacked (n_shuffles*n_trials, n_bins).
    """
    baseline = float(np.nanmean(shifted_snips))
    return real_snips - baseline


def get_time_above_simba_ci(
    real_snips,
    shifted_snips,
    percentile=97.5,
    start_bin=50,
    end_bin=150,
):
    """Calculate per-trial proportion of bins above the shuffled upper CI bound.

    The threshold is computed per-bin from the null distribution of individual
    shuffled snips (shape: n_shuffles*n_trials, n_bins), so the CI reflects
    single-trial variability at each timepoint in the infusion window.
    The returned scalar threshold is the mean across infusion-window bins
    (for logging only).
    """
    # Per-bin threshold computed only over the infusion window
    threshold_per_bin = np.nanpercentile(
        np.asarray(shifted_snips, dtype=float)[:, start_bin:end_bin],
        percentile,
        axis=0,
    )  # shape: (end_bin - start_bin,)

    proportions = []
    for snip in np.asarray(real_snips, dtype=float):
        window = np.asarray(snip[start_bin:end_bin], dtype=float)
        valid_mask = ~np.isnan(window)
        if not np.any(valid_mask):
            proportions.append(np.nan)
        else:
            proportions.append(
                float(np.mean(window[valid_mask] > threshold_per_bin[valid_mask]))
            )

    mean_threshold = float(np.nanmean(threshold_per_bin))
    return np.asarray(proportions, dtype=float), mean_threshold


def count_bins_above_threshold(snips, threshold, start_bin=50, end_bin=150):
    """Count bins above threshold in infusion window for each trial."""
    return np.array([
        np.sum(np.asarray(snip[start_bin:end_bin]) > threshold)
        for snip in snips
    ])
