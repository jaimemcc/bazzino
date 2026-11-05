"""Build the Zenodo archives for the photometry (TDT tank) data: one zip per rat.

Reads data/{10NaCl,45NaCl}_FileKey.csv, keeps only the sessions whose condition
is not in `conditions_to_exclude` (src/assemble_all_data.py), and zips the files
that tdt.read_block uses plus any videos (.avi/.mp4). *_extracted.mat and other
files are left out.
Originals are never modified.

Each zip contains tanks/<tank folder>/<files>, so extracting it into raw_data/
gives raw_data/tanks/<tank folder>, the default `tank_folder` in PARAMS.

Usage:
    pixi run python scripts/prepare_zenodo_tanks.py --source C:/path/to/tanks --out C:/path/to/staging
    pixi run python scripts/prepare_zenodo_tanks.py --source ... --out ... --dry-run
"""

import argparse
import hashlib
import json
import sys
import zipfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

TDT_SUFFIXES = {".tev", ".tsq", ".tbk", ".tdx", ".tin", ".tnt"}
# Not used by the pipeline, shipped as part of the raw data. Already compressed, so stored as is.
VIDEO_SUFFIXES = {".avi", ".mp4"}
SHIPPED_SUFFIXES = TDT_SUFFIXES | VIDEO_SUFFIXES
# Physiological-state labels in the file keys that the pipeline keeps
KEEP_STATES = {"Sodium Depleted", "Sodium Replete"}


def sessions_to_ship():
    keys = pd.concat(
        [pd.read_csv(ROOT / "data" / f"{n}NaCl_FileKey.csv") for n in (10, 45)],
        ignore_index=True,
    )
    keys = keys[keys["SessionInclude"] == "Y"]
    return keys[keys["Physiological state"].isin(KEEP_STATES)][["Subject", "Folder", "Physiological state"]]


def md5sum(path, chunk=1024 * 1024):
    digest = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--source", type=Path, required=True, help="Folder containing the tank folders")
    parser.add_argument("--out", type=Path, required=True, help="Folder to write the zips to")
    parser.add_argument("--dry-run", action="store_true", help="Report what would be zipped, write nothing")
    args = parser.parse_args()

    sessions = sessions_to_ship()
    missing = [f for f in sessions["Folder"] if not (args.source / f).is_dir()]
    if missing:
        sys.exit(f"Tank folders not found in {args.source}: {missing}")

    if not args.dry_run:
        args.out.mkdir(parents=True, exist_ok=True)

    manifest = []
    for subject, group in sessions.groupby("Subject"):
        files = []
        for folder in sorted(group["Folder"].unique()):
            tank_files = [p for p in (args.source / folder).iterdir() if p.suffix.lower() in SHIPPED_SUFFIXES]
            # Skip an .mp4 that duplicates an .avi of the same name
            avi_stems = {p.stem for p in tank_files if p.suffix.lower() == ".avi"}
            tank_files = [p for p in tank_files if not (p.suffix.lower() == ".mp4" and p.stem in avi_stems)]
            if not any(p.suffix.lower() == ".tev" for p in tank_files):
                sys.exit(f"No .tev file in {folder}")
            files += [(folder, p) for p in sorted(tank_files)]
        total = sum(p.stat().st_size for _, p in files)
        name = f"tanks_{subject}.zip"
        print(f"{name}: {group['Folder'].nunique()} tanks, {len(files)} files, {total / 1e6:.0f} MB uncompressed")
        if args.dry_run:
            continue
        zip_path = args.out / name
        with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=1) as zf:
            for folder, p in files:
                method = zipfile.ZIP_STORED if p.suffix.lower() in VIDEO_SUFFIXES else zipfile.ZIP_DEFLATED
                zf.write(p, f"tanks/{folder}/{p.name}", compress_type=method, compresslevel=1)
        manifest.append({
            "archive": name,
            "subject": subject,
            "tanks": sorted(group["Folder"].unique()),
            "size_bytes": zip_path.stat().st_size,
            "md5": md5sum(zip_path),
        })

    if not args.dry_run:
        (args.out / "tanks_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf8")
        print(f"\n{len(manifest)} archives, {sum(m['size_bytes'] for m in manifest) / 1e9:.2f} GB in {args.out}")


if __name__ == "__main__":
    main()
