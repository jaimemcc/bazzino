"""Download the raw data for this project from Zenodo.

Files are listed from the Zenodo API, downloaded into <dest>, verified against
the MD5 checksums Zenodo provides, and archives (.zip, .tar, .tar.gz) are
extracted in place. Re-running skips files that are already present and valid.

The record has two sets of archives:
  processed_data.zip   small (~21 MB): data/assembled, data/shap_values and the key files.
                       Extracted into --repo-root (default: current folder), which is all
                       the figure notebooks need.
  tanks_*.zip, dlc.zip, simba.zip
                       raw data (~34 GB). Extracted into <dest>, giving raw_data/tanks,
                       raw_data/dlc and raw_data/simba, the defaults in src/assemble_all_data.py.

Usage:
    pixi run python scripts/download_data.py --list          # just list the files
    pixi run python scripts/download_data.py --only processed  # small download, enough for the notebooks
    pixi run python scripts/download_data.py --only raw        # raw data only
    pixi run python scripts/download_data.py                   # everything
    pixi run python scripts/download_data.py --dest D:/data --no-extract
    pixi run python scripts/download_data.py --record 21394487
"""

import argparse
import hashlib
import json
import sys
import tarfile
import urllib.error
import urllib.request
import zipfile
from pathlib import Path

DEFAULT_RECORD = "21394487"
API = "https://zenodo.org/api/records/{record}"
CHUNK = 1024 * 1024


def fetch_file_list(record):
    url = API.format(record=record)
    try:
        with urllib.request.urlopen(url) as response:
            meta = json.load(response)
    except urllib.error.HTTPError as err:
        if err.code == 404:
            sys.exit(
                f"Record {record} was not found. It may not be published yet "
                f"(https://doi.org/10.5281/zenodo.{record})."
            )
        raise

    files = []
    for f in meta.get("files", []):
        # Zenodo has used both "key"/"links.self" and "filename"/"links.download".
        name = f.get("key") or f.get("filename")
        links = f.get("links", {})
        link = links.get("self") or links.get("download")
        checksum = f.get("checksum", "")
        files.append({
            "name": name,
            "size": f.get("size", 0),
            "md5": checksum.split(":", 1)[-1] if checksum else None,
            "url": link,
        })
    return files


def md5sum(path):
    digest = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(CHUNK), b""):
            digest.update(block)
    return digest.hexdigest()


def is_valid(path, entry):
    if not path.exists():
        return False
    if entry["md5"]:
        return md5sum(path) == entry["md5"]
    return path.stat().st_size == entry["size"]


def download(entry, path):
    tmp = path.with_suffix(path.suffix + ".part")
    done = 0
    with urllib.request.urlopen(entry["url"]) as response, open(tmp, "wb") as out:
        while True:
            block = response.read(CHUNK)
            if not block:
                break
            out.write(block)
            done += len(block)
            if entry["size"]:
                print(f"\r  {done / 1e6:8.1f} / {entry['size'] / 1e6:.1f} MB", end="", flush=True)
    print()
    tmp.replace(path)


def is_processed(name):
    return name.startswith("processed_data")


def extract(path, dest):
    name = path.name.lower()
    if name.endswith(".zip"):
        with zipfile.ZipFile(path) as zf:
            zf.extractall(dest)
    elif name.endswith((".tar", ".tar.gz", ".tgz")):
        with tarfile.open(path) as tf:
            tf.extractall(dest, filter="data")
    else:
        return False
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--record", default=DEFAULT_RECORD, help="Zenodo record ID")
    parser.add_argument("--dest", type=Path, default=Path("raw_data"), help="Destination folder for downloads and raw data")
    parser.add_argument("--repo-root", type=Path, default=Path("."), help="Where processed_data.zip is extracted")
    parser.add_argument("--only", choices=["raw", "processed"], help="Download only one set of archives")
    parser.add_argument("--list", action="store_true", help="List files and exit")
    parser.add_argument("--no-extract", action="store_true", help="Do not extract archives")
    args = parser.parse_args()

    files = fetch_file_list(args.record)
    if not files:
        sys.exit("The record has no files.")

    total = sum(f["size"] for f in files) / 1e9
    print(f"Record {args.record}: {len(files)} file(s), {total:.2f} GB")
    for f in files:
        print(f"  {f['name']}  ({f['size'] / 1e6:.1f} MB)")
    if args.list:
        return

    if args.only:
        files = [f for f in files if is_processed(f["name"]) == (args.only == "processed")]
        if not files:
            sys.exit(f"No '{args.only}' files in the record.")

    args.dest.mkdir(parents=True, exist_ok=True)
    for entry in files:
        path = args.dest / entry["name"]
        path.parent.mkdir(parents=True, exist_ok=True)
        if is_valid(path, entry):
            print(f"[skip] {entry['name']} already downloaded and verified")
        else:
            print(f"[get ] {entry['name']}")
            download(entry, path)
            if not is_valid(path, entry):
                path.unlink()
                sys.exit(f"Checksum mismatch for {entry['name']}; file removed. Try again.")
            print("  checksum OK")
        target = args.repo_root if is_processed(entry["name"]) else args.dest
        if not args.no_extract and extract(path, target):
            print(f"  extracted into {target}")

    print("Done.")


if __name__ == "__main__":
    main()
