import sys
import types
from pathlib import Path, PosixPath, PurePath, PurePosixPath, PureWindowsPath, WindowsPath

import dill


def enable_dill_pathlib_compat() -> None:
    if getattr(dill, "_pathlib_local_compat_enabled", False):
        return

    # Patch pandas StringDtype to handle old pickle format
    try:
        import pandas as pd
        from pandas.core.arrays.string_ import StringDtype as _OriginalStringDtype
        
        class StringDtypeCompat(_OriginalStringDtype):
            def __new__(cls, *args, **kwargs):
                # Old pandas StringDtype had a different signature
                # New version only accepts 'storage' parameter
                # Ignore extra positional arguments from old pickles
                if len(args) > 1:
                    args = args[:1]  # Keep only first argument (storage)
                return super().__new__(cls)
            
            def __init__(self, *args, **kwargs):
                # Old version accepted more arguments, new version accepts only storage
                if len(args) > 1:
                    # Ignore extra arguments from old pickle format
                    args = (args[0],) if args else ()
                super().__init__(*args, **kwargs)
        
        # Replace StringDtype in pandas modules
        pd.StringDtype = StringDtypeCompat
        pd.core.arrays.string_.StringDtype = StringDtypeCompat
        sys.modules['pandas'].StringDtype = StringDtypeCompat
        sys.modules['pandas.core.arrays.string_'].StringDtype = StringDtypeCompat
    except (ImportError, AttributeError):
        pass  # pandas not available or different structure

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
