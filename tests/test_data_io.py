import sys
from pathlib import Path, PureWindowsPath

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from data_io import ABSOLUTE_PATH_PLACEHOLDER, load_assembled, save_assembled


def make_data():
    rng = np.random.default_rng(0)
    return {
        "x_array": pd.DataFrame({
            "trial": np.arange(6),
            "id": pd.array(["a", "a", "b", "b", "c", None], dtype="string"),
            "score": [0.5, np.nan, 1.5, 2.5, -1.0, 3.0],
            "ok": [True, False, True, True, False, True],
        }),
        "fits_df_behav": pd.DataFrame({"k": [1.0, 2.0]}, index=[0, 3]),
        "snips_photo": rng.normal(size=(6, 10)),
        "snips_movement": np.array([[1.0, np.nan], [2.0, 3.0]]),
        "metadata": {"smoothed": True, "window": 5, "cols": ["a", "b"]},
        "params": {
            "data_folder": Path("data"),
            "tank_folder": PureWindowsPath("D:/personal/tanks"),
            "n_shuffles": np.int64(100),
            "threshold": np.float64(0.6),
            "flag": np.bool_(True),
            "nothing": None,
            "pool": [1, 2, 3],
        },
    }


def test_round_trip(tmp_path):
    data = make_data()
    save_assembled(data, tmp_path)
    loaded = load_assembled(tmp_path)

    assert list(loaded) == list(data)
    for key in ("snips_photo", "snips_movement"):
        assert loaded[key].dtype == data[key].dtype
        np.testing.assert_array_equal(loaded[key], data[key])
    for key in ("x_array", "fits_df_behav"):
        pd.testing.assert_frame_equal(
            loaded[key], data[key], check_column_type=False, check_dtype=False
        )
        assert list(loaded[key].index) == list(data[key].index)
    assert loaded["metadata"] == data["metadata"]


def test_string_columns_use_python_storage(tmp_path):
    save_assembled(make_data(), tmp_path)
    dtype = load_assembled(tmp_path)["x_array"]["id"].dtype
    assert isinstance(dtype, pd.StringDtype) and dtype.storage == "python"


def test_params_are_json_safe_and_hide_absolute_paths(tmp_path):
    save_assembled(make_data(), tmp_path)
    params = load_assembled(tmp_path)["params"]
    assert params["data_folder"] == "data"
    assert params["tank_folder"] == ABSOLUTE_PATH_PLACEHOLDER
    assert params["n_shuffles"] == 100 and params["flag"] is True
    assert params["nothing"] is None and params["pool"] == [1, 2, 3]


def test_unsupported_type_fails_loudly(tmp_path):
    with pytest.raises(TypeError):
        save_assembled({"bad": object()}, tmp_path)
