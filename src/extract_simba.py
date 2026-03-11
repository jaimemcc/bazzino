from pathlib import Path

import numpy as np
import pandas as pd


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


def make_simba_snips(probability_vector, solenoid_ts, fps=10, pre_bins=50, post_bins=150):
    """Create TTL-aligned Simba snips, matching behavior snip convention."""
    snips = []
    for i in range(len(solenoid_ts) - 1):
        frame_on = int(solenoid_ts[i] * fps)
        snips.append(
            _extract_single_snip(
                probability_vector,
                center_frame=frame_on,
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
):
    """Generate null-distribution snip means by circularly shifting the vector."""
    shifted = np.asarray(probability_vector, dtype=float).copy()
    shifted_means = []

    for _ in range(n_shuffles):
        shifted = np.roll(shifted, shift_frames)
        snips = make_simba_snips(
            shifted,
            solenoid_ts,
            fps=fps,
            pre_bins=pre_bins,
            post_bins=post_bins,
        )
        shifted_means.append(np.nanmean(snips, axis=0))

    return np.asarray(shifted_means, dtype=float)


def baseline_simba_snips(real_snips, shifted_means):
    """Subtract global shuffle mean from each real snip (notebook behavior)."""
    baseline = float(np.nanmean(shifted_means))
    return real_snips - baseline


def count_bins_above_threshold(snips, threshold, start_bin=50, end_bin=150):
    """Count bins above threshold in infusion window for each trial."""
    return np.array([
        np.sum(np.asarray(snip[start_bin:end_bin]) > threshold)
        for snip in snips
    ])
