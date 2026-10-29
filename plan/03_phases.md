# Phases

Work happens on this branch in the old repo until the result is ready; then the new repo is cut from it.

1. **Audit**: complete 02_audit_checklist.md and record the findings.
2. **Data format**: implement 05_data_format_migration.md (`data_io.py`, convert data, update writer and notebooks, verify round trip).
3. **Notebook cleanup**: remove `%autoreload`, remove personal paths, make cell 1 consistent, re-run all notebooks so the saved outputs are clean, confirm the figures are unchanged.
4. **Remove out-of-scope files**: video snippets, `tmp_shap_analysis.py`, `notebooks/text/`, notebook-workflow tooling, draft paper files, `.afdesign`, caches, and `setup.cfg` / `egg-info` if unused. Strip the matching tasks from `pixi.toml`.
5. **Docs and metadata**: README (citation, DOI, environment, how to reproduce each figure), `CITATION.cff`, Zenodo download script. `LICENSE` unchanged.
6. **Cut the new repo**: export the clean tree to a new directory, `git init`, one commit, check size (no cache pickles, no `.afdesign`), push to `mccutcheonlab/bazzino-paper`.
7. **Verify**: clone fresh, `pixi install`, download data from Zenodo, run every figure notebook, compare the figures with the paper.
8. **Release**: tag, Zenodo DOI, link from the article, mark the old repo as archived/superseded.
