"""Build the Zenodo archive with the processed data: processed_data.zip.

Contains, with internal prefix data/ (so extracting it in the repository root puts
everything where the notebooks expect it):
  data/assembled/              dataset used by the figure notebooks
  data/shap_values/            SHAP outputs
  data/{10,45}NaCl_FileKey.csv, data/{10,45}NaCl_SubjectKey.csv, data/ttls.csv

Usage:
    pixi run python scripts/prepare_zenodo_processed.py --out C:/path/staging
"""

import argparse
import json
import sys
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from prepare_zenodo_tanks import ROOT, md5sum

FOLDERS = ["data/assembled", "data/shap_values"]
FILES = [
    "data/10NaCl_FileKey.csv", "data/10NaCl_SubjectKey.csv",
    "data/45NaCl_FileKey.csv", "data/45NaCl_SubjectKey.csv", "data/ttls.csv",
]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", type=Path, required=True, help="Folder to write the zip to")
    args = parser.parse_args()

    files = [ROOT / f for f in FILES]
    for folder in FOLDERS:
        files += sorted(p for p in (ROOT / folder).rglob("*") if p.is_file())
    missing = [p for p in files if not p.exists()]
    if missing:
        sys.exit(f"Missing: {missing}")

    args.out.mkdir(parents=True, exist_ok=True)
    zip_path = args.out / "processed_data.zip"
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
        for p in files:
            zf.write(p, p.relative_to(ROOT).as_posix())

    entry = {
        "archive": zip_path.name,
        "files": [p.relative_to(ROOT).as_posix() for p in files],
        "size_bytes": zip_path.stat().st_size,
        "md5": md5sum(zip_path),
    }
    (args.out / "processed_data_manifest.json").write_text(json.dumps(entry, indent=2), encoding="utf8")
    print(f"{zip_path.name}: {len(files)} files, {entry['size_bytes'] / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
