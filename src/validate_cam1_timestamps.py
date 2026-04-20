"""Validate camera frame timing from TDT Cam1 onset timestamps.

This script is intentionally scoped to one tank first.
It loads Cam1 onset timestamps from a single TDT tank, computes np.diff,
and reports interframe interval statistics and potential dropped-frame gaps.

Examples:
    python src/validate_cam1_timestamps.py
    python src/validate_cam1_timestamps.py --tank-name <FolderName>
    python src/validate_cam1_timestamps.py --tank-path C:/path/to/tank
"""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import tdt


DEFAULT_TANK_FOLDER = Path("C:/Users/jmc010/Data/bazzino/tanks")
DEFAULT_DATA_FOLDER = Path("data")
DEFAULT_FILEKEY = "10NaCl_FileKey.csv"
DEFAULT_RESULTS_FOLDER = Path("results")

DEFAULT_FILEKEYS_ALL = [
    "10NaCl_FileKey.csv",
    "45NaCl_FileKey.csv",
]


def resolve_tank_path(
    tank_path: str | None,
    tank_name: str | None,
    tank_folder: Path,
    data_folder: Path,
    filekey_name: str,
) -> Path:
    """Resolve which single tank to analyze.

    Priority:
      1) explicit --tank-path
      2) --tank-name within --tank-folder
      3) first row from selected file key CSV
    """
    if tank_path:
        return Path(tank_path)

    if tank_name:
        return tank_folder / tank_name

    filekey_path = data_folder / filekey_name
    meta = pd.read_csv(filekey_path)
    if meta.empty:
        raise ValueError(f"No rows found in file key: {filekey_path}")

    first_tank_name = str(meta.iloc[0]["Folder"])
    return tank_folder / first_tank_name


def load_all_tank_rows(data_folder: Path, filekeys: list[str]) -> pd.DataFrame:
    """Load and merge all requested file keys, keeping one row per unique folder."""
    frames = []
    for filekey in filekeys:
        path = data_folder / filekey
        df = pd.read_csv(path).copy()
        df["filekey"] = filekey
        frames.append(df)

    meta = pd.concat(frames, ignore_index=True)
    if "Folder" not in meta.columns:
        raise KeyError("Expected column 'Folder' in file key(s).")

    # Keep first occurrence if the same folder appears across file keys.
    meta = meta.drop_duplicates(subset=["Folder"], keep="first").reset_index(drop=True)
    return meta


def get_cam1_onsets(block) -> np.ndarray:
    """Get Cam1 onset timestamps from a TDT block.

    Uses dict-like epoc access per project convention:
    data.epocs["Cam1"].onset
    """
    try:
        onsets = np.asarray(block.epocs["Cam1"].onset, dtype=float)
    except Exception as exc:
        available = []
        try:
            available = list(block.epocs.keys())
        except Exception:
            pass
        raise KeyError(
            "Could not read block.epocs['Cam1'].onset. "
            f"Available epocs: {available}"
        ) from exc

    if onsets.size < 2:
        raise ValueError("Cam1 onset array has fewer than 2 timestamps; cannot compute np.diff.")

    return onsets


def summarize_intervals(intervals: np.ndarray, long_factor: float) -> dict:
    """Compute summary stats and identify long interframe intervals."""
    expected = float(np.median(intervals))
    long_threshold = expected * long_factor
    long_mask = intervals > long_threshold
    long_idx = np.where(long_mask)[0]

    return {
        "n_intervals": int(intervals.size),
        "expected_interval_s": expected,
        "mean_interval_s": float(np.mean(intervals)),
        "std_interval_s": float(np.std(intervals)),
        "min_interval_s": float(np.min(intervals)),
        "max_interval_s": float(np.max(intervals)),
        "long_factor": long_factor,
        "long_threshold_s": float(long_threshold),
        "n_long_intervals": int(long_idx.size),
        "long_indices": long_idx,
    }


def summarize_single_tank(tank: Path, long_factor: float) -> tuple[np.ndarray, np.ndarray, dict]:
    """Load one tank and compute Cam1 onset interval summary."""
    if not tank.exists():
        raise FileNotFoundError(f"Tank path does not exist: {tank}")

    block = tdt.read_block(str(tank))
    onsets = get_cam1_onsets(block)
    intervals = np.diff(onsets)
    report = summarize_intervals(intervals, long_factor=long_factor)
    return onsets, intervals, report


def print_report(tank: Path, onsets: np.ndarray, intervals: np.ndarray, report: dict) -> None:
    """Print a concise human-readable timing report."""
    print("=" * 70)
    print("Cam1 Interframe Interval Validation")
    print("=" * 70)
    print(f"Tank: {tank}")
    print(f"Frames (onsets): {onsets.size}")
    print(f"Intervals (np.diff): {report['n_intervals']}")
    print()
    print("Interval summary (seconds)")
    print(f"  Expected (median): {report['expected_interval_s']:.6f}")
    print(f"  Mean:              {report['mean_interval_s']:.6f}")
    print(f"  Std:               {report['std_interval_s']:.6f}")
    print(f"  Min:               {report['min_interval_s']:.6f}")
    print(f"  Max:               {report['max_interval_s']:.6f}")
    print()
    print(
        "Long-interval rule: interval > "
        f"{report['long_factor']:.2f} x median = {report['long_threshold_s']:.6f} s"
    )
    print(f"Long intervals detected: {report['n_long_intervals']}")

    if report["n_long_intervals"] > 0:
        print()
        print("Top long intervals (up to first 20):")
        long_indices = report["long_indices"][:20]
        for idx in long_indices:
            # interval index i corresponds to onset[i] -> onset[i+1]
            dt = intervals[idx]
            t0 = onsets[idx]
            t1 = onsets[idx + 1]
            print(
                f"  i={idx:>7d}  dt={dt:.6f}s  "
                f"t0={t0:.6f}s  t1={t1:.6f}s"
            )


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Validate Cam1 interframe timing in one TDT tank by analyzing "
            "np.diff(epocs['Cam1'].onset)."
        )
    )
    parser.add_argument(
        "--tank-path",
        type=str,
        default=None,
        help="Full path to a specific tank to analyze.",
    )
    parser.add_argument(
        "--tank-name",
        type=str,
        default=None,
        help="Tank folder name under --tank-folder.",
    )
    parser.add_argument(
        "--tank-folder",
        type=Path,
        default=DEFAULT_TANK_FOLDER,
        help=f"Base directory containing tank folders (default: {DEFAULT_TANK_FOLDER}).",
    )
    parser.add_argument(
        "--data-folder",
        type=Path,
        default=DEFAULT_DATA_FOLDER,
        help=f"Project data folder containing file keys (default: {DEFAULT_DATA_FOLDER}).",
    )
    parser.add_argument(
        "--filekey",
        type=str,
        default=DEFAULT_FILEKEY,
        help=(
            "File key CSV used when neither --tank-path nor --tank-name is provided "
            f"(default: {DEFAULT_FILEKEY})."
        ),
    )
    parser.add_argument(
        "--long-factor",
        type=float,
        default=1.5,
        help=(
            "Flag interval as long if dt > median_dt * long_factor "
            "(default: 1.5)."
        ),
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help=(
            "Analyze all tanks found in 10NaCl_FileKey.csv and 45NaCl_FileKey.csv, "
            "and save CSV outputs in --results-folder."
        ),
    )
    parser.add_argument(
        "--results-folder",
        type=Path,
        default=DEFAULT_RESULTS_FOLDER,
        help=f"Where output CSV files are written (default: {DEFAULT_RESULTS_FOLDER}).",
    )
    parser.add_argument(
        "--summary-csv",
        type=str,
        default="cam1_interval_summary.csv",
        help="Filename for per-tank summary CSV.",
    )
    parser.add_argument(
        "--long-events-csv",
        type=str,
        default="cam1_long_intervals.csv",
        help="Filename for long-interval event-level CSV.",
    )
    return parser


def run_all_tanks(args: argparse.Namespace) -> None:
    """Analyze all tanks from file keys and store summary/event-level outputs."""
    print("Loading all tanks from file keys...")
    meta = load_all_tank_rows(args.data_folder, DEFAULT_FILEKEYS_ALL)
    print(f"Found {len(meta)} unique tank rows.")

    summary_rows: list[dict] = []
    event_rows: list[dict] = []

    for i, row in meta.iterrows():
        folder = str(row["Folder"])
        tank = args.tank_folder / folder
        subject = row.get("Subject", "")
        condition = row.get("Physiological state", "")
        treat = row.get("TreatNum", "")
        filekey = row.get("filekey", "")

        print(f"[{i + 1}/{len(meta)}] {folder}", end=" ... ")
        try:
            onsets, intervals, report = summarize_single_tank(tank, long_factor=args.long_factor)
            print(f"ok ({report['n_long_intervals']} long intervals)")

            summary_rows.append({
                "folder": folder,
                "tank_path": str(tank),
                "subject": subject,
                "condition": condition,
                "treat_num": treat,
                "filekey": filekey,
                "status": "ok",
                "error": "",
                "n_onsets": int(onsets.size),
                "n_intervals": report["n_intervals"],
                "expected_interval_s": report["expected_interval_s"],
                "mean_interval_s": report["mean_interval_s"],
                "std_interval_s": report["std_interval_s"],
                "min_interval_s": report["min_interval_s"],
                "max_interval_s": report["max_interval_s"],
                "long_factor": report["long_factor"],
                "long_threshold_s": report["long_threshold_s"],
                "n_long_intervals": report["n_long_intervals"],
            })

            for idx in report["long_indices"]:
                event_rows.append({
                    "folder": folder,
                    "subject": subject,
                    "condition": condition,
                    "treat_num": treat,
                    "filekey": filekey,
                    "interval_index": int(idx),
                    "t0_s": float(onsets[idx]),
                    "t1_s": float(onsets[idx + 1]),
                    "dt_s": float(intervals[idx]),
                })
        except Exception as exc:
            print(f"error ({exc})")
            summary_rows.append({
                "folder": folder,
                "tank_path": str(tank),
                "subject": subject,
                "condition": condition,
                "treat_num": treat,
                "filekey": filekey,
                "status": "error",
                "error": str(exc),
                "n_onsets": np.nan,
                "n_intervals": np.nan,
                "expected_interval_s": np.nan,
                "mean_interval_s": np.nan,
                "std_interval_s": np.nan,
                "min_interval_s": np.nan,
                "max_interval_s": np.nan,
                "long_factor": args.long_factor,
                "long_threshold_s": np.nan,
                "n_long_intervals": np.nan,
            })

    args.results_folder.mkdir(parents=True, exist_ok=True)
    summary_path = args.results_folder / args.summary_csv
    events_path = args.results_folder / args.long_events_csv

    summary_df = pd.DataFrame(summary_rows)
    events_df = pd.DataFrame(event_rows)

    run_ts = datetime.now().isoformat(timespec="seconds")
    summary_df.insert(0, "run_timestamp", run_ts)
    events_df.insert(0, "run_timestamp", run_ts)

    summary_df.to_csv(summary_path, index=False)
    events_df.to_csv(events_path, index=False)

    n_ok = int((summary_df["status"] == "ok").sum()) if not summary_df.empty else 0
    n_err = int((summary_df["status"] == "error").sum()) if not summary_df.empty else 0
    total_long = int(summary_df["n_long_intervals"].fillna(0).sum()) if not summary_df.empty else 0

    print()
    print("=" * 70)
    print("Batch Validation Complete")
    print("=" * 70)
    print(f"Tanks processed: {len(summary_df)}")
    print(f"OK: {n_ok}  Errors: {n_err}")
    print(f"Total long intervals: {total_long}")
    print(f"Summary CSV: {summary_path}")
    print(f"Long-events CSV: {events_path}")


def main() -> None:
    args = build_arg_parser().parse_args()

    if args.all:
        run_all_tanks(args)
        return

    tank = resolve_tank_path(
        tank_path=args.tank_path,
        tank_name=args.tank_name,
        tank_folder=args.tank_folder,
        data_folder=args.data_folder,
        filekey_name=args.filekey,
    )
    print(f"Loading TDT block: {tank}")
    onsets, intervals, report = summarize_single_tank(tank, long_factor=args.long_factor)
    print_report(tank=tank, onsets=onsets, intervals=intervals, report=report)


if __name__ == "__main__":
    main()
