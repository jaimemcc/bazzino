# Publication cleanup: overview

Branch: `plan/reformat-for-publication`
Status: planning. Questions are answered (see 04); nothing is implemented yet.

## Goal
Release a small, clean, reproducible repository alongside the published article, as a new repo (`mccutcheonlab/bazzino-paper`) with a fresh history. The current repo (`jaimemcc/bazzino`) stays untouched as the full backup.

## Current state (measured)
- 267 commits; packed `.git` is about 361 MiB.
- Large tracked files: `paper/figs/bazzino_figs_v03.afdesign` (~48 MB), `data/assembled_data.pickle` (~26 MB), `data/_cache_*.pickle` (3-10 MB each), `paper/*.pdf`, `paper/figs/*.png`.
- Clutter: `tmp_shap_analysis.py`, notebook-workflow tooling and scripts, `notebooks/text/`.

## History strategy: decided
Fresh single-commit history in a new repo (Option A). The old repo keeps the full history, so the history-rewriting options are dropped.

## Planning documents
- [01_open_questions.md](01_open_questions.md): the questions and your answers
- [02_audit_checklist.md](02_audit_checklist.md): what to inventory before changing anything
- [03_phases.md](03_phases.md): order of work
- [04_decisions_and_findings.md](04_decisions_and_findings.md): decisions, plus findings so far (notebook chatter, paths)
- [05_data_format_migration.md](05_data_format_migration.md): pickle to parquet/npz/JSON
