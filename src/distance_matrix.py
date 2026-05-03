import numpy as np
import pandas as pd

from scipy.ndimage import gaussian_filter1d
from scipy.spatial.distance import cdist
from sklearn.manifold import MDS

from figure_config import COLORS

DISTANCE_METRIC = "euclidean"
APPLY_TRIAL_SMOOTHING = True
SMOOTH_SIGMA = 1.2

GROUP_ORDER = [
    ("deplete", "10NaCl"),
    ("deplete", "45NaCl"),
    ("replete", "10NaCl"),
    ("replete", "45NaCl"),
]

GROUP_LABELS = {
    ("replete", "10NaCl"): "Replete 10NaCl",
    ("replete", "45NaCl"): "Replete 45NaCl",
    ("deplete", "10NaCl"): "Deplete 10NaCl",
    ("deplete", "45NaCl"): "Deplete 45NaCl",
}

GROUP_COLORS = {group: COLORS[idx] for idx, group in zip([2, 3, 0, 1], GROUP_ORDER)}


def make_metric_matrix(x_array, metric_col, GROUP_ORDER, GROUP_LABELS):
        
    scalar_df = x_array.loc[:, ["id", "condition", "infusiontype", "trial", metric_col]].copy()
    scalar_df = scalar_df.dropna(subset=[metric_col])
    scalar_df["group_key"] = list(zip(scalar_df["condition"], scalar_df["infusiontype"]))
    scalar_df = scalar_df[scalar_df["group_key"].isin(GROUP_ORDER)].copy()
    scalar_df["group_key"] = pd.Categorical(scalar_df["group_key"], categories=GROUP_ORDER, ordered=True)
    scalar_df = scalar_df.sort_values(["group_key", "id", "trial"])

    trial_counts = scalar_df.groupby(["group_key", "id"]).trial.nunique().sort_index()
    print("Trials per animal/group series:")
    print(trial_counts.groupby(level=0).agg(["count", "min", "max"]))

    metric_matrix_df = scalar_df.pivot_table(
        index=["group_key", "id"],
        columns="trial",
        values=metric_col,
        aggfunc="mean",
    ).sort_index()

    if metric_matrix_df.isna().any().any():
        missing_by_series = metric_matrix_df.isna().sum(axis=1)
        missing_by_series = missing_by_series[missing_by_series > 0]
        raise ValueError(
            "Some animal/group series are missing trial values. Missing counts:\n"
            f"{missing_by_series.to_string()}"
        )

    metric_matrix = metric_matrix_df.to_numpy(dtype=float)
    if APPLY_TRIAL_SMOOTHING:
        metric_matrix_for_distance = gaussian_filter1d(metric_matrix, sigma=SMOOTH_SIGMA, axis=1)
    else:
        metric_matrix_for_distance = metric_matrix.copy()

    series_index = list(metric_matrix_df.index)
    series_labels = [f"{GROUP_LABELS[group_key]} | {animal_id}" for group_key, animal_id in series_index]
    group_names = [GROUP_LABELS[group_key] for group_key, _ in series_index]
    group_to_row_indices = {
        group_key: [idx for idx, (series_group, _) in enumerate(series_index) if series_group == group_key]
        for group_key in GROUP_ORDER
    }
    group_sizes = {GROUP_LABELS[group_key]: len(indices) for group_key, indices in group_to_row_indices.items()}
    group_boundaries = np.cumsum([len(group_to_row_indices[group_key]) for group_key in GROUP_ORDER])

    print(f"Metric matrix shape: {metric_matrix.shape}")
    print(f"Distance input shape: {metric_matrix_for_distance.shape}")
    print("Series per group:")
    print(group_sizes)

    metric_matrix_df.head()
    return metric_matrix_df, series_labels, group_names, group_to_row_indices


def get_group_matrix_diffs(metric_matrix_for_distance, series_labels, group_to_row_indices):
    distances = cdist(metric_matrix_for_distance, metric_matrix_for_distance, metric=DISTANCE_METRIC)
    distances_df = pd.DataFrame(distances, index=series_labels, columns=series_labels)

    group_distance_matrix = pd.DataFrame(
        index=[GROUP_LABELS[group_key] for group_key in GROUP_ORDER],
        columns=[GROUP_LABELS[group_key] for group_key in GROUP_ORDER],
        dtype=float,
    )

    for row_group in GROUP_ORDER:
        row_idx = group_to_row_indices[row_group]
        for col_group in GROUP_ORDER:
            col_idx = group_to_row_indices[col_group]
            group_distance_matrix.loc[GROUP_LABELS[row_group], GROUP_LABELS[col_group]] = distances[np.ix_(row_idx, col_idx)].mean()

    print(f"Distance metric: {DISTANCE_METRIC}")
    print(f"Mean off-diagonal distance: {distances[~np.eye(distances.shape[0], dtype=bool)].mean():.4f}")
    
    return group_distance_matrix