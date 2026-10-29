# Data format migration: pickle -> parquet / npz / JSON

Decision: replace `assembled_data.pickle` (dill) with portable formats. This is a code change, and the only refactor in scope before release.

## What is in `assembled_data.pickle` (measured)
A dict with 20 keys:
- 2-D `float64` arrays, shape (1960, 200): `snips_photo`, `snips_simba`, `snips_simba_raw`, `snips_simba_shifted_baseline`, `snips_simba_zscore`, `snips_simba_baseline`, `snips_movement`, `snips_angvel`, `pca_transformed`
- DataFrames with no object columns: `x_array` (1960x30), `z_dep45` (490x30), `fits_df*` (5 small ones)
- Plain dicts: `metadata` (34 scalar and list entries), `params` (71 entries). `params` contains `WindowsPath` objects and personal paths like `D:\TestData\...`. These are why `pickle_compat.py` exists.
- Not all 20 keys have been listed yet; the audit must list the rest (the inspection output was cut at 90 lines).

## Proposed layout
```
data/assembled/
  arrays.npz         all ndarrays, keyed by name (compressed)
  x_array.parquet    one parquet file per DataFrame
  z_dep45.parquet
  fits_df*.parquet
  metadata.json      analysis settings used by the notebooks
  params.json        paths converted to strings; personal paths removed or redacted
```
- Parquet for tables (typed, portable, readable from R/Julia). npz for numeric arrays (they are not tables). JSON for settings (human-readable).
- A single module `src/data_io.py` with `save_assembled(data, folder)` and `load_assembled(folder)`. It returns the same dict as today, so notebook code after the load line stays unchanged.
- Each notebook's `dill.load` block, `import dill` and `enable_dill_pathlib_compat` lines are replaced by one `load_assembled(DATAFOLDER / "assembled")` call.

## Things to check
- **Round-trip fidelity:** old dict == new dict (dtypes, column order, index, NaNs, `metadata` types). Index and dtype preservation in parquet is the usual trap, as is pandas 3.0's string dtype.
- **`params` paths:** decide which paths are kept. Anything with personal paths is stripped before the data is committed, so nothing personal ends up in the new repo.
- **Dependencies:** `pyarrow` is not in the pixi environment (checked); it must be added. `dill` can be removed only if the caches also stop using it (see below).
- **Caches:** `_cache_*.pickle` are written and read by `src/assemble_all_data.py` (lines ~224, ~248). They are not shipped. Leave them as pickle (they are private, regenerable) and keep `dill` for that, or convert them too. Default: leave them, and decide after the audit.
- **Other readers:** `src/plots_for_spectral.py` uses `pickle`; `figure_Sx_analyse_shap_values` imports dill. Check what they load.
- **Writer:** `src/assemble_all_data.py` line ~1805 (`dill.dump(output, f)`) must write the new format (and, for a while, optionally the old one for comparison).

## Plan
1. Write `src/data_io.py`, plus a test that converts the current pickle and checks the round trip.
2. Convert the pickle with it and check the result.
3. Update `assemble_all_data.py` to write the new format.
4. Update the 11 figure notebooks to load it, re-run them, and compare the figures with the current ones.
5. Remove `pickle_compat.py` if nothing else needs it.
6. Record the sizes and reload time against the old pickle.

This is done on the working branch before the new repo is cut, so the new repo gets the final result.
