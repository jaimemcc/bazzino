# Audit findings (first pass)

## Notebook dependencies
- All figure notebooks except `Sx_analyse_shap_values` load `data/assembled_data.pickle`. The SHAP notebook builds its input with `assemble_shap_dfs` from `data/shap_values/`.
- Extra data files read by notebooks: `results/figure{1,2}_mixed_anova_jasp_{long,wide}.csv`, `results/transition_aligned_auc_summary.csv`, `results/transition_sigmoid_fits_deplete_45NaCl.csv`, `results/figure_Sx_withintrial_crosscorr_best_lag_r.csv`. Several are notebook outputs that other notebooks read back, so the notebook run order matters. Document it in the README.
- `src/` modules used by notebooks: `distance_matrix`, `extract_behav_parameters`, `figure_config`, `figure_plotting`, `model_fit_helpers`, `stats_helpers`, `realignment_helpers`, `bootstrap_and_shuffle_helpers`, `utils`, `assemble_shap_dfs`, `pickle_compat`.
- `figure_1` and `figure_S1` import one threshold function from `assemble_all_data` (cell near the end; other mentions are only text). That module is about 1800 lines and needs raw-data paths, so move that function into a small shared module (e.g. `utils`) or keep it importable without side effects. `extract_simba` is used by `assemble_all_data` and `color_schemes` by `figure_config`, so both stay. `hybrid_feature_extractor` has no importers at all (not even `assemble_all_data`); candidate for removal.

## `src/` modules not imported by any notebook
`add_target_columns`, `color_schemes`, `extract_simba`, `filter_dlc_bodyparts`, `find_folders_without_avis`, `hybrid_feature_extractor`, `make_snipped_videos`, `plots_for_spectral`, `validate_cam1_timestamps`.
- `extract_simba` and `hybrid_feature_extractor` may be used by `assemble_all_data` (the pipeline from raw data). Check before deleting.
- `make_snipped_videos` is out (decided). The rest need a keep/remove decision, in particular `color_schemes`, `find_folders_without_avis`, `plots_for_spectral`, `validate_cam1_timestamps` (hard-coded `C:/Users/jmc010/...` tank path).

## Packaging
- `setup.cfg` lists `py_modules` (a flat install of `src/`). `pyproject.toml` has the jupytext setting (not needed once `.py` copies are removed). Notebooks import with `sys.path.insert(0, "../src")`, so the package install is probably unnecessary. Confirm with `pixi.toml` before removing.

## Figures
- `paper/figs/` has 16 PNGs. Fig 5 is named `Fig5_NaCl integration.png` (a space in the name). `FigS5_optimal_lag.png` and `FigS5_transition_correlations.png` share a number. Check which are final and rename consistently, after you've confirmed which supplements are in the paper.
- Notebooks also write PDFs (`legend_pies_*.pdf`, `logistic_fits_45NaCl.pdf`). Check where.

## Secrets
- A search for key/password/secret/token in `.py`, `.toml`, `.json` and `.md` files found nothing. Personal paths (not secrets) remain, listed in 04.

## Still to do
`.gitignore` for the new repo, `pixi.lock` check, the two remaining pickle keys, SHAP data usage, the remaining notes files in `paper/` and `src/`.
