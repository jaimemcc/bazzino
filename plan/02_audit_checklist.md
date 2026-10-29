# Audit checklist (to do before any changes)

Status key: [x] done, [ ] to do

- [x] First-cell notebook "chatter" cause found (see 04)
- [x] Structure of `assembled_data.pickle` inspected (see 05; list the remaining keys)
- [ ] Map each figure notebook to the `src/` modules and `data/` files it uses
- [ ] Find unused `src/` modules and dead code (candidates: `plots_for_spectral.py`, `find_folders_without_avis.py`, `filter_dlc_bodyparts.py`, `validate_cam1_timestamps.py`, `add_target_columns.py`: decide per file whether it belongs with the paper)
- [ ] List which files each notebook writes (`results/`, `paper/figs/`)
- [ ] Find hard-coded and personal paths in `src/`, notebooks and saved outputs (known: `assemble_all_data.py`, `validate_cam1_timestamps.py`, `params` in the pickle)
- [ ] Scan the tree for secrets, credentials and personal info
- [ ] Review `.gitignore` for the new repo (caches, `.pixi/`, `*.egg-info`, `__pycache__`)
- [ ] Check what `setup.cfg`, `pyproject.toml` and `bazzino.egg-info` are used for (is `src` imported as a package or by `sys.path`? Notebooks use `sys.path.insert`)
- [ ] Check that `pixi.lock` reproduces the environment (add `pyarrow`)
- [ ] Decide which `paper/figs/` files are final figures; look at `paper/feature_and_shap_notes.md`, `individual_heatmaps.*` and `src/feature_classification_notes.md` (in or out?)
- [ ] Check what `data/shap_values/` is produced by, and whether `Sx_analyse_shap_values` needs all five runs
