"""Add a target column to behavior CSV files.

Usage examples:
    # TTL-derived target column, writing to <csv_dir>/target_files/<same_name>.csv
    python src/add_target_columns.py data/my_file.csv -t data/ttls.csv

    # Override timing and output column, writing to a custom directory
    python src/add_target_columns.py data/my_file.csv -d 5 -i 120 -c target -o data/custom_targets -w

    # Write all zeros (skip reading TTLs)
    python src/add_target_columns.py data/my_file.csv -z -c target -w
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def get_stub(ttls: pd.DataFrame, csv_path: Path) -> str:
    """Find the TTL column whose name appears in the input csv filename."""
    for col in ttls.columns:
        if col in csv_path.name:
            return col
    raise ValueError(
        f"No TTL column name from {list(ttls.columns)} was found in '{csv_path.name}'."
    )


def get_frames(
    ttls: pd.DataFrame,
    stub: str,
    delay_frames: int = 0,
    infusion_frames: int = 100,
) -> tuple[list[int], list[int]]:
    ttl_onsets = pd.to_numeric(ttls.loc[:, stub], errors="coerce").dropna()
    frames_on = [int(onset * 10) + delay_frames for onset in ttl_onsets]
    frames_off = [frame + infusion_frames for frame in frames_on]
    return frames_on, frames_off


def make_vector(frames_on: list[int], frames_off: list[int], total_frames: int) -> np.ndarray:
    vector = np.zeros(total_frames, dtype=int)
    for on, off in zip(frames_on, frames_off):
        clipped_on = max(0, on)
        clipped_off = min(total_frames, off)
        if clipped_on < clipped_off:
            vector[clipped_on:clipped_off] = 1
    return vector


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Add a binary target column to a behavior CSV using TTL onset times, "
            "or write all zeros."
        )
    )
    parser.add_argument(
        "csv_file",
        type=Path,
        help="Path to input behavior CSV file to augment.",
    )
    parser.add_argument(
        "-t",
        "--ttl-file",
        type=Path,
        default=Path("data/ttls.csv"),
        help="Path to TTL CSV file (ignored when --all-zeros is used).",
    )
    parser.add_argument(
        "-d",
        "--delay-frames",
        type=int,
        default=0,
        help="Frame delay to add to TTL onset frame index.",
    )
    parser.add_argument(
        "-i",
        "--infusion-frames",
        type=int,
        default=100,
        help="Number of frames to mark as 1 after each onset.",
    )
    parser.add_argument(
        "-c",
        "--output-column",
        type=str,
        default="target",
        help="Name of output target column.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=None,
        help=(
            "Directory where output CSV will be written. Default: a 'target_files' "
            "subdirectory next to the input csv_file."
        ),
    )
    parser.add_argument(
        "-w",
        "--overwrite",
        action="store_true",
        help="Overwrite output file if it already exists.",
    )
    parser.add_argument(
        "-z",
        "--all-zeros",
        action="store_true",
        help="Write all zeros to output column and skip reading TTLs.",
    )
    return parser.parse_args()


def resolve_default_output_dir(csv_path: Path) -> Path:
    return csv_path.parent / "target_files"


def main() -> int:
    args = parse_args()

    csv_path = args.csv_file
    if not csv_path.exists():
        raise FileNotFoundError(f"Input CSV not found: {csv_path}")

    output_dir = args.output_dir or resolve_default_output_dir(csv_path)
    if output_dir.exists() and not output_dir.is_dir():
        raise NotADirectoryError(f"Output path exists but is not a directory: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / csv_path.name

    if output_path.exists() and not args.overwrite:
        raise FileExistsError(
            f"Output already exists: {output_path}. Use --overwrite to replace it."
        )

    csv_df = pd.read_csv(csv_path)
    total_frames = len(csv_df)

    if args.all_zeros:
        vector = np.zeros(total_frames, dtype=int)
        frames_on: list[int] = []
        frames_off: list[int] = []
        stub = "N/A (all-zeros mode)"
    else:
        if args.infusion_frames < 0:
            raise ValueError("--infusion-frames must be >= 0")
        ttls = pd.read_csv(args.ttl_file)
        stub = get_stub(ttls, csv_path)
        frames_on, frames_off = get_frames(
            ttls,
            stub,
            delay_frames=args.delay_frames,
            infusion_frames=args.infusion_frames,
        )
        vector = make_vector(frames_on, frames_off, total_frames)

    csv_df[args.output_column] = vector
    csv_df.to_csv(output_path, index=False)

    positive_frames = int(vector.sum())
    n_infusions = len(frames_on)
    print("Target column generation complete")
    print(f"Input CSV:      {csv_path}")
    print(f"Output CSV:     {output_path}")
    print(f"Output dir:     {output_dir}")
    print(f"Output column:  {args.output_column}")
    print(f"Mode:           {'all-zeros' if args.all_zeros else 'ttl-derived'}")
    print(f"TTL column:     {stub}")
    print(f"Total frames:   {total_frames}")
    print(f"Positive frames:{positive_frames}")
    print(f"Infusions:      {n_infusions}")
    if not args.all_zeros:
        print(f"Delay frames:   {args.delay_frames}")
        print(f"Infusion frames:{args.infusion_frames}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
