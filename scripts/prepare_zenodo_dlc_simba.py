"""Build the Zenodo archives for the DeepLabCut and SimBA outputs.

Selects the files for the sessions that ship (same selection as
prepare_zenodo_tanks.py) using the pipeline's own file-matching rules:
  DLC:   PB_NAapp-<date>_<tank folder>*.csv   (extract_behav_parameters.read_DLC_csv)
  SimBA: *_<tank folder>_*.csv                (extract_simba.find_simba_file)

Writes dlc.zip and simba.zip, with internal prefixes dlc/ and simba/, so extracting
them into raw_data/ gives raw_data/dlc and raw_data/simba (the defaults in PARAMS).

Usage:
    pixi run python scripts/prepare_zenodo_dlc_simba.py --dlc C:/path/output_csv_shuffle4 \
        --simba C:/path/Output_all_animals_appetitive --out C:/path/staging
"""

import argparse
import hashlib
import json
import sys
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from prepare_zenodo_tanks import md5sum, sessions_to_ship


def dlc_files(folder, stub):
    date = stub.split("-")[1]
    return sorted(folder.glob(f"PB_NAapp-{date}_{stub}*.csv"))


def simba_files(folder, stub):
    matches = sorted(folder.glob(f"*_{stub}_*.csv"))
    cam1 = [p for p in matches if "Cam1" in p.name]
    return cam1[:1] if cam1 else matches[:1]


def build(name, prefix, folder, finder, stubs, out, dry_run):
    files = []
    for stub in stubs:
        found = finder(folder, stub)
        if len(found) != 1:
            sys.exit(f"{name}: expected 1 file for {stub} in {folder}, found {len(found)}")
        files.append(found[0])
    total = sum(p.stat().st_size for p in files)
    print(f"{prefix}.zip: {len(files)} files, {total / 1e6:.0f} MB uncompressed")
    if dry_run:
        return None
    zip_path = out / f"{prefix}.zip"
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
        for p in files:
            zf.write(p, f"{prefix}/{p.name}")
    return {
        "archive": zip_path.name,
        "files": [p.name for p in files],
        "size_bytes": zip_path.stat().st_size,
        "md5": md5sum(zip_path),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dlc", type=Path, required=True, help="Folder with the DLC csv files")
    parser.add_argument("--simba", type=Path, required=True, help="Folder with the SimBA csv files")
    parser.add_argument("--out", type=Path, required=True, help="Folder to write the zips to")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    stubs = sorted(sessions_to_ship()["Folder"].unique())
    if not args.dry_run:
        args.out.mkdir(parents=True, exist_ok=True)

    manifest = [
        build("DLC", "dlc", args.dlc, dlc_files, stubs, args.out, args.dry_run),
        build("SimBA", "simba", args.simba, simba_files, stubs, args.out, args.dry_run),
    ]
    if not args.dry_run:
        (args.out / "dlc_simba_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf8")
        print(f"\n{sum(m['size_bytes'] for m in manifest) / 1e9:.2f} GB in {args.out}")


if __name__ == "__main__":
    main()
