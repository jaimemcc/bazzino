import warnings

import numpy as np
import pandas as pd

from scipy.ndimage import gaussian_filter1d
from scipy.spatial.distance import cdist
from sklearn.manifold import MDS
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.metrics import confusion_matrix, classification_report, accuracy_score

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


def prepare_metric_classification_data(
    x_array,
    metric_col,
    group_order=GROUP_ORDER,
    group_labels=GROUP_LABELS,
    smooth_features=False,
):
    """Build per-series trial features for classification.

    This uses the same representation as the distance/MDS workflow:
    rows are (group_key, id) series and columns are trial values.

    Returns
    -------
    X : np.ndarray
        Feature matrix with shape (n_series, n_trials).
    y : np.ndarray
        Group label per row (human-readable labels from group_labels).
    rat_groups : np.ndarray
        Rat IDs per row (useful for group-aware CV splitting).
    metric_matrix_df : pd.DataFrame
        DataFrame form of the trial feature matrix.
    """

    scalar_df = x_array.loc[:, ["id", "condition", "infusiontype", "trial", metric_col]].copy()
    scalar_df = scalar_df.dropna(subset=[metric_col])
    scalar_df["group_key"] = list(zip(scalar_df["condition"], scalar_df["infusiontype"]))
    scalar_df = scalar_df[scalar_df["group_key"].isin(group_order)].copy()
    scalar_df["group_key"] = pd.Categorical(scalar_df["group_key"], categories=group_order, ordered=True)
    scalar_df = scalar_df.sort_values(["group_key", "id", "trial"])

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

    X = metric_matrix_df.to_numpy(dtype=float)
    if smooth_features:
        X = gaussian_filter1d(X, sigma=SMOOTH_SIGMA, axis=1)

    series_index = list(metric_matrix_df.index)
    y = np.array([group_labels[group_key] for group_key, _ in series_index])
    rat_groups = np.array([animal_id for _, animal_id in series_index])

    return X, y, rat_groups, metric_matrix_df


def evaluate_group_classifier_leave_one_rat_out(
    X,
    y,
    rat_groups,
    labels=None,
    normalize="true",
):
    """Evaluate a logistic classifier with leave-one-rat-out CV.

    Parameters
    ----------
    X, y, rat_groups : array-like
        Feature matrix, class labels, and grouping labels for CV.
    labels : list-like or None
        Class order for confusion matrix/report. If None, inferred from y.
    normalize : {'true', 'pred', 'all', None}
        Normalization mode for confusion matrix.

    Returns
    -------
    dict
        Keys: accuracy, report_text, y_true, y_pred, cm, labels.
    """

    if labels is None:
        labels = list(np.unique(y))

    clf = Pipeline([
        ("imputer", SimpleImputer(strategy="median")),
        ("scaler", StandardScaler()),
        ("model", LogisticRegression(max_iter=5000, solver="lbfgs")),
    ])

    logo = LeaveOneGroupOut()
    y_true_all, y_pred_all = [], []

    for train_idx, test_idx in logo.split(X, y, groups=rat_groups):
        clf.fit(X[train_idx], y[train_idx])
        y_pred_fold = clf.predict(X[test_idx])
        y_true_all.extend(y[test_idx])
        y_pred_all.extend(y_pred_fold)

    y_true_all = np.array(y_true_all)
    y_pred_all = np.array(y_pred_all)

    acc = accuracy_score(y_true_all, y_pred_all)
    report_text = classification_report(y_true_all, y_pred_all, labels=labels, zero_division=0)
    cm = confusion_matrix(y_true_all, y_pred_all, labels=labels, normalize=normalize)

    return {
        "accuracy": acc,
        "report_text": report_text,
        "y_true": y_true_all,
        "y_pred": y_pred_all,
        "cm": cm,
        "labels": labels,
    }


def classify_groups_from_metric(
    x_array,
    metric_col,
    group_order=GROUP_ORDER,
    group_labels=GROUP_LABELS,
    smooth_features=False,
    expected_n_trials=None,
    normalize="true",
):
    """Prepare per-rat trial features and evaluate group classification.

    This is a convenience wrapper for reuse in notebooks.
    """

    X, y, rat_groups, metric_matrix_df = prepare_metric_classification_data(
        x_array=x_array,
        metric_col=metric_col,
        group_order=group_order,
        group_labels=group_labels,
        smooth_features=smooth_features,
    )

    labels = [group_labels[g] for g in group_order if group_labels[g] in set(y)]
    results = evaluate_group_classifier_leave_one_rat_out(
        X=X,
        y=y,
        rat_groups=rat_groups,
        labels=labels,
        normalize=normalize,
    )

    if expected_n_trials is not None and X.shape[1] != expected_n_trials:
        print(f"Warning: expected {expected_n_trials} trial features, found {X.shape[1]}")

    results["X"] = X
    results["y"] = y
    results["rat_groups"] = rat_groups
    results["metric_matrix_df"] = metric_matrix_df
    return results


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
    
    return distances, group_distance_matrix

def compute_mds_coordinates(distances):
    distance_scale = np.nanmax(distances)

    if distance_scale > 0:
        distances_for_mds = distances / distance_scale
    else:
        distances_for_mds = distances.copy()

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=FutureWarning, module="sklearn.manifold._mds")
        mds = MDS(
            n_components=2,
            metric=True,
            dissimilarity="precomputed",
            random_state=42,
            n_init=8,
            max_iter=500,
            init="random",
        )
        coords = mds.fit_transform(distances_for_mds)
        
        return coords