"""Load SHAP results and assemble feature-level summary dataframes."""

from pathlib import Path
import re

import numpy as np
import pandas as pd


_SHAP_DROP_COLUMNS = [
    "Unnamed: 0",
    "Prediction_probability",
    "Sum",
    "Expected_value",
]


def _strip_reverse_suffix(feature: str) -> str:
    """Strip trailing reverse-coding suffixes to recover the base feature name."""
    for suffix in ("_percentile_rank", "_deviation"):
        if feature.endswith(suffix):
            feature = feature[: -len(suffix)]
    return feature


def assign_group(feature: str) -> str:
    """Assign a feature to movement, geometry, tortuosity, or other.

    Reverse-coding suffixes (_deviation, _percentile_rank) are stripped first
    so a feature's semantic group reflects the underlying signal, not whether
    it has been reverse-coded (see `is_reverse_coded` for that).
    """
    base = _strip_reverse_suffix(feature)
    if re.search(r"^Tortuosity_", base):
        return "tortuosity"
    if re.search(
        r"^(Movement_mouse_(nose|tail_base|left_ear|right_ear|head_base)|"
        r"Total_movement_all_bodyparts(_M1)?|Total_movement_M1_|"
        r"Tail_base_movement_M1_|Head_base_movement_M1_|Nose_movement_M1_)",
        base,
    ):
        return "movement"
    if re.search(
        r"^(Mouse_nose_to_tail|Mouse_head_to_tail|Mouse_Ear_distance|"
        r"M1_|Mouse1_(smallest|largest|mean)_euclid_distances_)",
        base,
    ):
        return "geometry"
    return "other"


def assign_bodypart(feature: str) -> str:
    """Extract the body-part label encoded in a feature name."""
    feature_lower = feature.lower()
    if "all_bodyparts" in feature_lower:
        return "whole rat"
    for key, label in [
        ("nose", "nose"),
        ("tail_base", "tail base"),
        ("head_base", "head base"),
        ("ear", "ears"),
    ]:
        if key in feature_lower:
            return label
    if feature_lower.startswith(("mouse_", "mouse1_")) or "total" in feature_lower:
        return "whole rat"
    return "other"


def get_time_window(feature: str) -> str:
    """Extract the rolling-window value from a feature name.

    Looks at the very end of the name, allowing the number to be followed by
    a reverse-coding suffix, so e.g. `..._mean_10_percentile_rank` still
    resolves to window "10" instead of "none".
    """
    match = re.search(r"_(\d+(?:\.\d+)?)(?:_deviation|_percentile_rank)*$", feature)
    return match.group(1) if match else "none"


def is_reverse_coded(feature: str) -> bool:
    """Identify relative-to-average features using the project naming convention."""
    return feature.endswith(("_deviation", "_percentile_rank"))


def _find_duplicate_percentile_rank_columns(df_shap_raw: pd.DataFrame) -> list[str]:
    """Find _percentile_rank columns that exactly duplicate a _deviation column.

    Some rolling-window features compute _percentile_rank with the same
    "mean - current" formula as _deviation instead of a true percentile rank
    (see hybrid_feature_extractor.py), so they carry no new information.
    """
    duplicates = []
    for column in df_shap_raw.columns:
        if not column.endswith("_percentile_rank"):
            continue
        deviation_column = column[: -len("_percentile_rank")] + "_deviation"
        if deviation_column not in df_shap_raw.columns:
            continue
        if pd.to_numeric(df_shap_raw[column], errors="coerce").equals(
            pd.to_numeric(df_shap_raw[deviation_column], errors="coerce")
        ):
            duplicates.append(column)
    return duplicates


def _read_shap_tables(datafolder: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    shap_frames = []
    raw_frames = []
    folders = sorted(path for path in datafolder.iterdir() if path.is_dir())

    for folder in folders:
        shap_frames.append(pd.read_csv(folder / "SHAP_values_Appetitive.csv"))
        raw_frames.append(pd.read_csv(folder / "RAW_SHAP_feature_values_Appetitive.csv"))

    if not shap_frames:
        raise FileNotFoundError(f"No SHAP result folders found in {datafolder}")

    df_shap = pd.concat(shap_frames, axis=0, ignore_index=True)
    df_shap_raw = pd.concat(raw_frames, axis=0, ignore_index=True)
    if len(df_shap) != len(df_shap_raw):
        raise ValueError("SHAP and raw SHAP tables have different numbers of rows")

    df_shap = df_shap.drop(columns=_SHAP_DROP_COLUMNS, errors="ignore")
    df_shap_raw = df_shap_raw.drop(columns=["Unnamed: 0"], errors="ignore")

    # Drop _percentile_rank columns that are exact duplicates of _deviation.
    duplicate_columns = _find_duplicate_percentile_rank_columns(df_shap_raw)
    df_shap = df_shap.drop(columns=duplicate_columns, errors="ignore")
    df_shap_raw = df_shap_raw.drop(columns=duplicate_columns, errors="ignore")
    return df_shap, df_shap_raw


def _make_feature_summary(
    df_shap: pd.DataFrame,
    df_shap_raw: pd.DataFrame,
    target_column: str,
    target_value: int,
) -> pd.DataFrame:
    feature_columns = [column for column in df_shap.columns if column != target_column]
    shap_all = df_shap[feature_columns]
    raw_all = df_shap_raw[feature_columns].apply(
        pd.to_numeric,
        errors="coerce",
    )
    target_mask = df_shap[target_column] == target_value
    non_target_mask = ~target_mask
    shap_target = df_shap.loc[target_mask, feature_columns]
    shap_non_target = df_shap.loc[non_target_mask, feature_columns]

    # Correlation between each feature's raw value and its own SHAP value.
    # A negative correlation means a *lower* raw value pushes SHAP (and the
    # prediction) up, e.g. movement features where less movement means more
    # appetitive; a positive correlation means a *higher* raw value pushes
    # the prediction up.
    raw_shap_corr = shap_all.corrwith(raw_all).reindex(feature_columns)

    feature_summary = pd.DataFrame({
        "feature": feature_columns,
        "mean_shap": shap_all.mean().to_numpy(),
        "importance": shap_all.abs().mean().to_numpy(),
        "variability": shap_all.std().to_numpy(),
        "mean_shap_appetitive": shap_target.mean().to_numpy(),
        "mean_shap_non_appetitive": shap_non_target.mean().to_numpy(),
        "raw_mean": raw_all.mean().reindex(feature_columns).to_numpy(),
        "raw_sd": raw_all.std().reindex(feature_columns).to_numpy(),
        "raw_min": raw_all.min().reindex(feature_columns).to_numpy(),
        "raw_max": raw_all.max().reindex(feature_columns).to_numpy(),
        "raw_shap_corr": raw_shap_corr.to_numpy(),
    })
    feature_summary["class_difference"] = (
        feature_summary["mean_shap_appetitive"]
        - feature_summary["mean_shap_non_appetitive"]
    )
    # Signed importance: magnitude from mean |SHAP|, sign from the direction
    # of the raw-value relationship, so a feature that must go *down* to
    # drive the prediction up gets a negative value here.
    feature_summary["signed_importance"] = feature_summary["importance"] * np.sign(
        feature_summary["raw_shap_corr"]
    )
    feature_summary["group"] = feature_summary["feature"].apply(assign_group)
    feature_summary["bodypart"] = feature_summary["feature"].apply(assign_bodypart)
    feature_summary["timewindow"] = feature_summary["feature"].apply(get_time_window)
    feature_summary["reverse_coded"] = feature_summary["feature"].apply(is_reverse_coded)
    # Same as signed_importance, but with the sign flipped back for
    # reverse-coded features so it reflects the direction of the original
    # (non-reverse-coded) underlying signal, comparable across both.
    feature_summary["aligned_signed_importance"] = feature_summary["signed_importance"] * feature_summary[
        "reverse_coded"
    ].map({True: -1, False: 1})
    # Group split by reverse-coding, e.g. "movement" vs "movement, reverse coded".
    feature_summary["display_group"] = feature_summary["group"] + feature_summary["reverse_coded"].map(
        {True: ", reverse coded", False: ""}
    )
    feature_summary = feature_summary.sort_values("importance", ascending=False).reset_index(drop=True)
    feature_summary["cumulative_importance"] = (
        feature_summary["importance"].cumsum() / feature_summary["importance"].sum()
    )
    feature_summary["reverse_cumulative_importance"] = (
        1 - feature_summary["cumulative_importance"]
    )
    return feature_summary


def assemble_shap_dfs(
    datafolder: str | Path,
    target_column: str = "Appetitive",
    target_value: int = 1,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Load SHAP tables and return the raw tables plus a feature summary."""
    df_shap, df_shap_raw = _read_shap_tables(Path(datafolder))
    feature_summary = _make_feature_summary(
        df_shap,
        df_shap_raw,
        target_column,
        target_value,
    )
    return df_shap, df_shap_raw, feature_summary
