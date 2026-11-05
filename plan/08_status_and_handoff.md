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

## Addendum
- Unused imports pruned from all 9 notebooks (AST check across all cells; re-executed clean). Scripts in session files/: unused_imports.py (report), prune_imports.py (apply). If notebooks are edited again, re-run the report.
- Notebook outputs no longer carry absolute paths: jasp CSV prints are repo-relative, and expected statsmodels warnings are silenced in src/figure_config.py. Re-running in VS Code should stay clean; still grep tracked content for 'jmc010' before the final cut.
- Zenodo tanks (2026-11-05): `scripts/prepare_zenodo_tanks.py` builds one zip per rat (`tanks_PBxx.zip`, internal prefix `tanks/<folder>/`, so extracting into `raw_data/` gives `raw_data/tanks/<folder>`). Only the TDT files read by `tdt.read_block` (.tev .tsq .Tbk .Tdx .tin .tnt) plus videos (.avi/.mp4, unused by the pipeline but part of the raw data) are included; `_extracted.mat` and `StoresListing.txt` are left out. Only Depleted + Replete sessions are shipped (20 rats, 40 tanks); Replete Experienced and Thirsty are dropped by `conditions_to_exclude` anyway. Output: 20 zips, 22.8 GB (videos exist for only 23 of the 40 tanks in the source folder; videos are stored, not compressed), plus `tanks_manifest.json` (md5s), in `C:\Users\jmc010\Data\bazzino\zenodo_staging` (outside the repo). `PB64-221027-105535` (Replete Experienced) is in the key but missing on disk, harmless.
- Zenodo DLC + SimBA (2026-11-05): `scripts/prepare_zenodo_dlc_simba.py` builds `dlc.zip` (40 csv, prefix `dlc/`) and `simba.zip` (40 csv, prefix `simba/`), same 40 sessions, with `dlc_simba_manifest.json`. Sources: DLC = `Data\bazzino\output_csv_shuffle4` (its `behav_stats` reproduce the stored ones exactly; `from_paula\Sodium_Appetite_DLC`, the older shuffle1 model, is unused and the pipeline cannot read it), SimBA = `simba_preds\Output_all_animals_appetitive` (`quick_test` is a strict subset, not needed). Staging folder total: 24 files, 22.05 GB (the 20 tank zips are ~21.2 GB).
- Rebuild test (2026-11-05, scratch folder, since deleted): extracted the staged zips (without videos) and ran `run_pipeline` with data/raw paths pointed at it. Completed in ~400 s. All 20 keys match the existing `data/assembled`: arrays bit-identical, DataFrames equal, `metadata` equal; `params` differs only in `data_folder`/`output_*` paths; `behav_stats` has 40 vs 48 entries (the 8 missing are Replete Experienced sessions; common entries equal). The 18 logged errors are all for the dropped sessions. This also exercised `save_assembled` on real data. (An earlier note in this file about ~1e-4 photometry differences was wrong: that was comparing raw snips with the z-scored ones in the dataset.)
- Videos complete (2026-11-05): the 18 missing `.avi` files were copied from `D:\Missing_files_copy` into the tank folders (all 40 shipped tanks now have one), the duplicate `_Cam1.mp4` in `PB71-221123-113609` is skipped by `prepare_zenodo_tanks.py`, and the 20 tank zips were rebuilt: 32.81 GB (staging total 33.66 GB, 25 files incl. dlc/simba zips and manifests; under Zenodo's 50 GB limit). Each zip checked: 2 tanks, 2 avi, no mp4/mat, integrity OK. Four SimBA-like CSVs in `D:\Missing_files_copy` (PB39 220627, PB23 220601, PB62 221022/221024) were not used.
- Remaining before upload: write the Zenodo description / README for the record (file list, layout, how to use `download_data.py`).
- Keys trimmed (2026-11-05): the Replete Experienced and Thirsty rows were removed from `data/{10,45}NaCl_FileKey.csv` (30->20, 28->20 rows), their columns from `data/ttls.csv` (57->40), and their entries from `data/assembled/behav_stats.json` (48->40; `data/assembled` is untracked, local edit). The full keys are in the old repo's history (`main`). A rebuild from the staged zips with the trimmed keys matches `data/assembled` on all 20 keys with zero errors in the log. `prepare_zenodo_tanks.py` still filters by state, now a no-op.
- Added scripts/download_data.py (tested on a public record; the project record 21394487 was still 404 on 2026-10-29, so run it with --list first tomorrow and check archive names/layout match raw_data/{tanks,dlc,simba}) and scripts/precut_check.py (run before the cut: personal paths, forbidden files, notebook errors, sizes, required files; currently passes).
- Panel outputs renamed figSx_* -> figS5_* (individual transitions) and figS6_nacl_da_trimmed. CITATION.cff license set to GPL-3.0-only (LICENSE is plain GPL-3.0 text; change to -or-later if intended).
