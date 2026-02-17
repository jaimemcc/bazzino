import sys
import types
from pathlib import Path, PosixPath, PurePath, PurePosixPath, PureWindowsPath, WindowsPath

import dill


def enable_dill_pathlib_compat() -> None:
    if getattr(dill, "_pathlib_local_compat_enabled", False):
        return

    original_load = dill.load

    def _dill_load_with_pathlib_compat(file_obj, *args, **kwargs):
        try:
            return original_load(file_obj, *args, **kwargs)
        except ModuleNotFoundError as error:
            if "pathlib._local" not in str(error):
                raise

            local_mod = types.ModuleType("pathlib._local")
            local_mod.Path = Path
            local_mod.PosixPath = PosixPath
            local_mod.WindowsPath = WindowsPath
            local_mod.PurePath = PurePath
            local_mod.PurePosixPath = PurePosixPath
            local_mod.PureWindowsPath = PureWindowsPath
            sys.modules["pathlib._local"] = local_mod

            return original_load(file_obj, *args, **kwargs)

    dill.load = _dill_load_with_pathlib_compat
    dill._pathlib_local_compat_enabled = True
