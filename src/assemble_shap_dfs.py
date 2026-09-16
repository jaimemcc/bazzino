"""Load SHAP results and assemble feature-level summary dataframes."""

from pathlib import Path
import re

import pandas as pd


_SHAP_DROP_COLUMNS = [
    "Unnamed: 0",
    "Prediction_probability",
    "Sum",
    "Expected_value",
]


def assign_group(feature: str) -> str:
    """Assign a feature to movement, geometry, or other."""
    if re.search(
        r"^(Movement_mouse_(nose|tail_base|left_ear|right_ear|head_base)|"
        r"Total_movement_all_bodyparts_M1|Total_movement_M1_|"
        r"Tail_base_movement_M1_|Head_base_movement_M1_|Nose_movement_M1_)"
        r"|(_deviation|_percentile_rank)$",
        feature,
    ):
        return "movement"
    if re.search(
        r"^(Mouse_nose_to_tail|Mouse_head_to_tail|Mouse_Ear_distance|"
        r"M1_|Mouse1_(smallest|largest|mean)_euclid_distances_)",
        feature,
    ):
        return "geometry"
    return "other"


def assign_bodypart(feature: str) -> str:
    """Extract the body-part label encoded in a feature name."""
    feature_lower = feature.lower()
    if "all_bodyparts" in feature_lower:
        return "all body parts"
    for key, label in [
        ("nose", "nose"),
        ("tail_base", "tail base"),
        ("head_base", "head base"),
        ("left_ear", "left ear"),
        ("ear_left", "left ear"),
        ("right_ear", "right ear"),
        ("ear_right", "right ear"),
    ]:
        if key in feature_lower:
            return label
    if feature_lower.startswith(("mouse_", "mouse1_")) or "total" in feature_lower:
        return "whole mouse"
    return "other"


def get_time_window(feature: str) -> str:
    """Extract the trailing rolling-window value from a feature name."""
    match = re.search(r"_(\d+(?:\.\d+)?)$", feature)
    return match.group(1) if match else "none"


def is_reverse_coded(feature: str) -> bool:
    """Identify relative-to-average features using the project naming convention."""
    return feature.endswith(("_deviation", "_percentile_rank"))


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
    })
    feature_summary["class_difference"] = (
        feature_summary["mean_shap_appetitive"]
        - feature_summary["mean_shap_non_appetitive"]
    )
    feature_summary["group"] = feature_summary["feature"].apply(assign_group)
    feature_summary["bodypart"] = feature_summary["feature"].apply(assign_bodypart)
    feature_summary["timewindow"] = feature_summary["feature"].apply(get_time_window)
    feature_summary["reverse_coded"] = feature_summary["feature"].apply(is_reverse_coded)
    feature_summary = feature_summary.sort_values("importance", ascending=False).reset_index(drop=True)
    feature_summary["cumulative_importance"] = (
        feature_summary["importance"].cumsum() / feature_summary["importance"].sum()
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
