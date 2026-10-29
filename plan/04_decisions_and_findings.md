# Decisions (from answered questions) and audit findings

## Decisions
- **History**: Option A. New repo `mccutcheonlab/bazzino-paper` with a fresh single-commit history. `jaimemcc/bazzino` stays as is, as the backup.
- **Data in repo**: everything in `data/` except the caches (`_cache_*.pickle`). Raw data is not in the repo; it is kept locally and uploaded to Zenodo. A download script is wanted.
- **Out**: `make_snipped_videos` (notebook and `src/` script), `tmp_shap_analysis.py`, `notebooks/text/*.py`, `tooling/notebook_workflow/`, `setup_notebook_workflow.py`, `update_notebook_workflow.py`, `notebook_workflow_config.json`, the notebook-workflow tasks in `pixi.toml`, `paper/paper.*`, `paper/todo_list.md`, `.afdesign`.
- **Data format**: switch from pickle to parquet/npz/JSON (see 05). This is the one code change before release.
- **Kept scripts**: the `src/` scripts not imported by notebooks (`add_target_columns`, `extract_simba`, `filter_dlc_bodyparts`, `find_folders_without_avis`, `hybrid_feature_extractor`, `plots_for_spectral`, `validate_cam1_timestamps`, `color_schemes`) stay because other parts of the project use them. Each must be documented in `src/README.md` (purpose, inputs, usage). Personal paths in them still get removed. `make_snipped_videos` stays out, as decided earlier.
- **Figures**: remove `FigS5_optimal_lag.png` and the three `FigSx_*.png` files (FR1_offset, Prog_offset, epochs). Keep Fig1-5, FigS1-S4, `FigS5_transition_correlations.png` and `FigS6_trimmed_NaCl.png`. Rename `Fig5_NaCl integration.png` to remove the space.
- **Notebooks removed now (local files and their `notebooks/text/*.py`)**: `figure_Sx_analyse_shap_values`, `figure_Sx_plotted_by_NaCl`, `figure_Sx_transitions`, `figure_Sx_withintrial_correlations`, `make_snipped_videos`. Consequences to check in the audit: `data/shap_values/` is still needed if `assemble_shap_dfs` is imported by `figure_S1` (it is); results CSVs written only by removed notebooks (`transition_sigmoid_fits_deplete_45NaCl.csv`, `figure_Sx_withintrial_crosscorr_best_lag_r.csv`) are probably no longer needed.
- **Prototype code**: some notebooks end with prototype cells that will be cut in the final version. Not yet marked; to be listed per notebook later (e.g. `figure_S4`: cells from `rats[0,2,8]` on).
- **Important**: `notebooks/**/*.ipynb` is gitignored in this repo (line 186 of `.gitignore`), so the notebooks with outputs are NOT in git history; only the jupytext `.py` copies are. The new repo must un-ignore `.ipynb`, and the old repo does not back up the notebooks. A copy of the migrated notebooks is in the session folder; make your own backup too.
- **Figure panels**: notebooks write panels to `paper/figs/panels/` (PDF); the PNGs in `paper/figs/` are the assembled figures (Affinity file). Decide whether the panels folder ships.
- **In**: the figure notebooks including the `Sx` ones (with outputs kept), final figures in `paper/figs/`, `LICENSE` (unchanged).
- **To add**: README (citation, DOI, environment setup, per-figure reproduction), `CITATION.cff`, a Zenodo download script.
- **Environment**: pixi only; remove `setup.cfg` and `bazzino.egg-info` if unneeded.
- **Refactoring**: not now. The notebook first-cell noise is fixed.
- **Verification**: a clean clone must reproduce the figures before release.

## Findings

### The first-cell "chatter"
- The first code cell of the figure notebooks (e.g. `figure_1_behaviour_simba.ipynb`, `figure_2_photometry.ipynb`) has a saved stderr output. `%autoreload 2` (IPython's deduperreload) tries to read stdlib files such as `functools.py` and `shlex.py` with the Windows default cp1252 encoding. It fails with `UnicodeDecodeError` and prints a traceback for each file.
- Likely fixes:
  1. Set `PYTHONUTF8=1` in the pixi environment (`[activation.env]` in `pixi.toml`) so that Python reads files as UTF-8.
  2. Drop `%autoreload` from the release notebooks. It is a development convenience, not needed to reproduce figures.
  3. Re-run so the saved outputs are clean.
- Recommendation: do 2, and 1 as well if it is cheap. Also, `figure_1` uses `%load_ext IPython.extensions.autoreload` while the others use `%load_ext autoreload`, so make that consistent.

### Hard-coded paths in saved outputs
- Saved outputs in figures 1, 2, S4 and `Sx_individual_transitions` print `C:\Users\jmc010\...`. These get removed by re-running from the new repo.
- Code with personal paths: `src/assemble_all_data.py` (lines ~78-82, commented laptop paths) and `src/validate_cam1_timestamps.py` (`DEFAULT_TANK_FOLDER`).

### Data flow
- Every figure notebook loads `data/assembled_data.pickle` through `dill`, with the `pickle_compat` shim for pandas/pathlib.
- `src/assemble_all_data.py` builds it from the raw tanks, DLC and SimBA folders, using the `_cache_*.pickle` files as intermediate caches.
- Notebook outputs total ~1.4 MB, so keeping outputs is cheap.

## Resolved: `assembled_data.pickle`
- It is the only input of the figure notebooks and can only be rebuilt from raw data (not in the repo).
- Decision: ship it in the new repo, but converted to portable formats (parquet for tables, npz for arrays, JSON for settings) instead of a dill pickle. Details in 05. The `_cache_*.pickle` files stay out.
