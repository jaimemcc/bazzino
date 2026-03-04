"""Filter DLC CSV to keep only selected bodyparts.

Usage examples:
    # Keep nose and left ear, single file with _filtered suffix
    python src/filter_dlc_bodyparts.py data/dlc_file.csv -b nose,l_ear -s _filtered

    # Keep multiple bodyparts, all CSVs in a directory
    python src/filter_dlc_bodyparts.py data/dlc_folder -b nose,l_ear,r_ear -s _filtered

    # Process directory recursively with overwrite
    python src/filter_dlc_bodyparts.py data/dlc_folder -b nose -r -w

    # Single file, custom output path
    python src/filter_dlc_bodyparts.py data/dlc_file.csv -b nose,l_ear,r_ear -o data/custom.csv -w
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def keep_dlc_bodyparts(
    input_csv: Path,
    bodyparts_to_keep: list[str],
    output_csv: Path | None = None,
    bodypart_level: int = 1,
) -> None:
    """Keep only selected bodyparts in a DLC CSV and write back with the same header layout."""
    input_csv = Path(input_csv)
    output_csv = Path(output_csv) if output_csv is not None else input_csv

    # DLC exports typically use a 3-row column header: scorer/bodyparts/coords.
    df = pd.read_csv(input_csv, header=[0, 1, 2])

    keep = set(bodyparts_to_keep)
    cols_to_keep = []
    for col in df.columns:
        # Keep metadata-like header columns plus requested bodyparts.
        if len(col) > bodypart_level and (
            col[bodypart_level] in keep or col[bodypart_level] == "bodyparts"
        ):
            cols_to_keep.append(col)

    if not cols_to_keep:
        raise ValueError(f"No columns matched requested bodyparts: {sorted(keep)}")

    filtered = df.loc[:, cols_to_keep]
    filtered.to_csv(output_csv, index=False)

    kept = sorted(
        {
            c[bodypart_level]
            for c in cols_to_keep
            if len(c) > bodypart_level and c[bodypart_level] != "bodyparts"
        }
    )
    print("DLC filter complete")
    print(f"Input file:      {input_csv}")
    print(f"Output file:     {output_csv}")
    print(f"Bodyparts kept ({len(kept)}): {sorted(kept)}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Filter DLC CSV(s) to keep only selected bodyparts. Works on single files or entire directories."
    )
    parser.add_argument(
        "input_path",
        type=Path,
        help="Path to input DLC CSV file or directory.",
    )
    parser.add_argument(
        "-b",
        "--bodyparts",
        type=str,
        required=True,
        help="Comma-separated list of bodyparts to keep (e.g., 'nose,l_ear,r_ear').",
    )
    parser.add_argument(
        "-s",
        "--suffix",
        type=str,
        default="",
        help="Suffix to append to output filenames (e.g., '_filtered'). Default: none (overwrites).",
    )
    parser.add_argument(
        "-o",
        "--output-file",
        type=Path,
        default=None,
        help="Custom output file path (single file mode only). If provided, --suffix is ignored.",
    )
    parser.add_argument(
        "-r",
        "--recursive",
        action="store_true",
        help="Recursively search subdirectories for CSV files (directory mode only).",
    )
    parser.add_argument(
        "-w",
        "--overwrite",
        action="store_true",
        help="Allow overwriting existing output files.",
    )
    return parser.parse_args()


def resolve_output_path(input_csv: Path, suffix: str) -> Path:
    """Construct output path by inserting suffix before .csv extension."""
    stem = input_csv.stem
    parent = input_csv.parent
    return parent / f"{stem}{suffix}.csv"


def main() -> int:
    args = parse_args()

    input_path = args.input_path
    if not input_path.exists():
        raise FileNotFoundError(f"Input path not found: {input_path}")

    # Parse comma-separated bodyparts list
    bodyparts = [bp.strip() for bp in args.bodyparts.split(",")]
    if not bodyparts or not all(bodyparts):
        raise ValueError("--bodyparts must contain at least one non-empty bodypart name.")

    # Single file mode
    if input_path.is_file():
        if not input_path.suffix.lower() == ".csv":
            raise ValueError(f"Input file must be a CSV: {input_path}")
        
        output_csv = args.output_file or resolve_output_path(input_path, args.suffix)
        if output_csv.exists() and not args.overwrite:
            raise FileExistsError(
                f"Output file already exists: {output_csv}. Use --overwrite to replace it."
            )
        keep_dlc_bodyparts(input_path, bodyparts, output_csv=output_csv)
        return 0

    # Directory mode
    if input_path.is_dir():
        if args.output_file:
            raise ValueError("--output-file cannot be used in directory mode.")
        
        pattern = "*.csv" if not args.recursive else "**/*.csv"
        csv_files = list(input_path.glob(pattern))
        
        if not csv_files:
            raise FileNotFoundError(f"No CSV files found in {input_path}")
        
        print(f"Processing {len(csv_files)} CSV file(s) from {input_path}")
        processed = 0
        for csv_file in csv_files:
            try:
                output_csv = resolve_output_path(csv_file, args.suffix)
                if output_csv.exists() and not args.overwrite:
                    print(f"Skipping {csv_file.name} (output already exists)")
                    continue
                keep_dlc_bodyparts(csv_file, bodyparts, output_csv=output_csv)
                processed += 1
            except Exception as e:
                print(f"Error processing {csv_file}: {e}")
        
        print(f"\nProcessed {processed}/{len(csv_files)} files successfully.")
        return 0

    raise ValueError(f"Input path is neither a file nor a directory: {input_path}")


if __name__ == "__main__":
    raise SystemExit(main())
