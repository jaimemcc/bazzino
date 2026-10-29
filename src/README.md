# `src/` guide

Code is grouped by role. Figure notebooks (in `../notebooks/`) import from here by adding `../src` to `sys.path`.

## Data pipeline

| File | Purpose |
|---|---|
| `assemble_all_data.py` | Builds the analysis dataset from raw TDT photometry tanks, DeepLabCut (DLC) CSVs and SimBA predictions. Writes `data/assembled/` via `data_io.save_assembled`. Needs the raw data (see the top-level README); configure paths in `PARAMS` at the top of the file. Intermediate caches (`data/_cache_*.pickle`) are private and regenerable. |
| `data_io.py` | `save_assembled` / `load_assembled`: portable storage (parquet + npz + json) for the assembled data. Notebooks load data with `load_assembled(DATAFOLDER / "assembled")`. |
| `extract_behav_parameters.py` | TTLs, DLC reading, interpolation of low-likelihood points, movement/angular velocity and behaviour snippets. |
| `extract_simba.py` | Reads SimBA classifier probabilities and builds SimBA-aligned snippets and summaries. |
| `assemble_shap_dfs.py` | Loads SHAP results (`data/shap_values/`) and builds feature-level summary tables. |

## Analysis and plotting helpers (used by notebooks)

| File | Purpose |
|---|---|
| `figure_config.py` | Shared paths, colours and matplotlib style for all figures. |
| `figure_plotting.py` | Reusable plotting functions. |
| `color_schemes.py` | Colour palettes (used by `figure_config`). |
| `stats_helpers.py` | Wilcoxon, mixed-model and multiple-comparison helpers. |
| `model_fit_helpers.py` | Logistic/sigmoid and other trial-wise fits and model selection. |
| `bootstrap_and_shuffle_helpers.py` | Bootstrap and shuffle null distributions. |
| `realignment_helpers.py`, `utils.py` | Realign behaviour to the fitted dopamine transition. |
| `distance_matrix.py` | Distance matrices and MDS embeddings for trial similarity. |

## Standalone scripts (not used by the figure notebooks)

These were used in other stages of the project (video/behaviour preparation and
quality control) and are kept for completeness. They are not needed to
reproduce the figures from `data/assembled/`.

| File | What it does |
|---|---|
| `add_target_columns.py` | CLI. Adds a `target` column to behaviour CSVs for SimBA, derived from TTLs (`data/ttls.csv`) or all zeros. See the docstring for examples. |
| `filter_dlc_bodyparts.py` | CLI. Keeps only selected body parts in DLC CSVs (single file or folder, optionally recursive). |
| `find_folders_without_avis.py` | Helper that lists folders under a root that contain no `.avi` files. |
| `hybrid_feature_extractor.py` | Custom SimBA feature extractor. Requires the `simba` package, which is **not** in the pixi environment. |
| `plots_for_spectral.py` | Plotting for spectral-clustering exploration of behavioural features. Loads plain pickles and uses many scikit-learn tools. |
| `validate_cam1_timestamps.py` | Checks TDT Cam1 inter-frame intervals from a tank and reports possible dropped frames. Needs raw tanks (`--tank-path`). |
