"""Read and write the assembled analysis data in portable formats.

Layout of an assembled-data folder:
    manifest.json    key order and storage kind of every entry
    arrays.npz       all numpy arrays
    <key>.parquet    one file per DataFrame
    <key>.json       dicts (metadata, params, behav_stats, ...)

load_assembled() returns the same dict structure that the notebooks used to
get from the legacy pickle.
"""

import json
from pathlib import Path, PurePath

import numpy as np
import pandas as pd

MANIFEST = "manifest.json"
ARRAYS = "arrays.npz"
ABSOLUTE_PATH_PLACEHOLDER = "<absolute path removed>"


def _json_default(obj):
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError(f"Cannot serialise {type(obj).__name__} to JSON")


def _sanitise(obj):
    """Make nested containers JSON-friendly; drop absolute (personal) paths."""
    if isinstance(obj, PurePath):
        return ABSOLUTE_PATH_PLACEHOLDER if obj.is_absolute() else obj.as_posix()
    if isinstance(obj, dict):
        return {str(k): _sanitise(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_sanitise(v) for v in obj]
    return obj


def _use_python_strings(df):
    """Return string columns with python-backed storage.

    The analysis was developed with python-backed strings. pandas 3 switches to
    pyarrow-backed strings when pyarrow is installed (needed for parquet), and
    those do not support some numpy-style indexing the notebooks rely on.
    """
    py_str = pd.StringDtype("python", na_value=pd.NA)
    cols = [c for c in df.columns if isinstance(df[c].dtype, pd.StringDtype)]
    return df.astype({c: py_str for c in cols}) if cols else df


def save_assembled(data, folder):
    """Write `data` (dict of arrays, DataFrames and plain dicts) to `folder`."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)

    manifest = {}
    arrays = {}
    for key, value in data.items():
        if isinstance(value, np.ndarray):
            arrays[key] = value
            manifest[key] = {"kind": "array", "file": ARRAYS}
        elif isinstance(value, pd.DataFrame):
            fname = f"{key}.parquet"
            value.to_parquet(folder / fname)
            manifest[key] = {"kind": "dataframe", "file": fname}
        elif isinstance(value, dict):
            fname = f"{key}.json"
            with open(folder / fname, "w", encoding="utf-8") as f:
                json.dump(_sanitise(value), f, indent=1, default=_json_default)
            manifest[key] = {"kind": "dict", "file": fname}
        else:
            raise TypeError(f"Unsupported type for key '{key}': {type(value).__name__}")

    if arrays:
        np.savez_compressed(folder / ARRAYS, **arrays)
    with open(folder / MANIFEST, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=1)


def load_assembled(folder):
    """Load a folder written by save_assembled() back into a dict."""
    folder = Path(folder)
    with open(folder / MANIFEST, encoding="utf-8") as f:
        manifest = json.load(f)

    data = {}
    npz = np.load(folder / ARRAYS) if any(m["kind"] == "array" for m in manifest.values()) else None
    try:
        for key, meta in manifest.items():
            if meta["kind"] == "array":
                data[key] = npz[key]
            elif meta["kind"] == "dataframe":
                data[key] = _use_python_strings(pd.read_parquet(folder / meta["file"]))
            else:
                with open(folder / meta["file"], encoding="utf-8") as f:
                    data[key] = json.load(f)
    finally:
        if npz is not None:
            npz.close()
    return data
