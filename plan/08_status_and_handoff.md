# Status and handoff (last updated 2026-10-29)

Branch: `plan/reformat-for-publication` in the old repo (`jaimemcc/bazzino`, left untouched on `main` as the backup).
Target: new repo `mccutcheonlab/bazzino-paper`, one fresh commit.

## Phase status (see 03_phases.md)

| Phase | Status |
|---|---|
| 1 Audit | Done (06_audit_findings.md) |
| 2 Data format | Done except end-to-end run of the writer (needs raw data). `data_io.py` + 4 passing tests (`pixi run test`) |
| 3 Notebook cleanup | Done: autoreload removed, personal paths scrubbed, all 9 notebooks run clean, outputs saved |
| 4 Remove out-of-scope | Done (tooling, caches, .afdesign files, notes, draft paper, dead figures untracked/removed) |
| 5 Docs and metadata | README, src/README.md, CITATION.cff done. **Zenodo download script not written** |
| 6 Cut new repo | Not started |
| 7 Verify from fresh clone | Not started |
| 8 Release | Not started |

## Decisions made (all)

- Fresh-history new repo; old repo is the archive.
- Ship: `data/assembled/` (parquet/npz/json), file/subject key CSVs, `data/ttls.csv`, `data/shap_values/`, 4 `results/figure{1,2}_mixed_anova_jasp_*.csv`, final figures `paper/figs/Fig*.png`, notebooks with outputs, LICENSE.
- Do not ship: caches (`data/_cache_*.pickle`, stay `dill`, private), `assembled_data.pickle`, `.afdesign`, `paper/figs/panels/`, notes files, draft paper, notebook-workflow tooling.
- All unused `src/` scripts kept and documented in `src/README.md`.
- Notebook `figure_Sx_individual_transitions` was renamed `figure_S5_individual_transitions`. `transition_aligned_auc_summary.csv` not needed.
- Zenodo DOI: 10.5281/zenodo.21394487 (record NOT yet published as of today; raw data to be uploaded). Paper still a bioRxiv preprint; update README + CITATION.cff `preferred-citation` once accepted.

## Gotchas

- `data/assembled/` is **untracked** in the old repo on purpose (19.7 MB). It must be copied into the new repo. Check `git status` shows it before the cut. It is not gitignored.
- Run Python via `pixi run python ...`; the env python crashes outside activation. For notebooks set `PYTHONUTF8=1`.
- pandas 3 + pyarrow gives pyarrow-backed strings; `load_assembled` casts back to python storage. Do not reintroduce `pickle_compat`.
- Scratch scripts are in the session folder `files/` (not in the repo).
- Possible clash: `paper/figs/FigS5_transition_correlations.png` vs the notebook now called S5. Confirm numbering with the paper.
- `CITATION.cff` license field is "See LICENSE" (not valid SPDX); replace with the real ID. Repo URL is assumed.
- Some functions in `src/figure_plotting.py` are only used by cut notebook cells; left for a later refactor.
- `figure_S4` cell 3 is a possible duplicate (user hasn't ruled).

## Tomorrow: once raw data are on Zenodo

1. Publish the Zenodo record; confirm DOI resolves and note file names/layout and checksums.
2. Write `scripts/download_data.py` (or similar): download from the record, verify checksums, extract to `raw_data/{tanks,dlc,simba}` (the defaults in `PARAMS` of `src/assemble_all_data.py`). Mention it in the README (currently the README just says download and extract).
3. Download the raw data and run `pixi run python src/assemble_all_data.py`. This exercises the new `save_assembled` writer for the first time on real data. Expect cache files to be created (slow). Then compare the new `data/assembled/` with the current one (arrays should be equal, DataFrames equal; `params` will differ only in paths). Fix any failure in the writer (`_sanitise` may meet types not seen before).
4. Re-run the 9 notebooks (`PYTHONUTF8=1`, `jupyter nbconvert --execute --inplace`), scrub absolute paths from outputs again (the scrub script is in the session `files/scrub.py`), compare figures with `paper/figs`.
5. Phase 6: export the tree with `git archive` (tracked files) plus `data/assembled/`, into a new directory; `git init`; one commit; check size and that no caches/`.afdesign` are present; push to `mccutcheonlab/bazzino-paper`.
6. Phase 7: fresh clone, `pixi install`, `pixi run test`, run notebooks.
7. Phase 8: tag, update the Zenodo record to link the code, archive the old repo, update citation after acceptance.
