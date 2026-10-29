# Notebook run status

Run on 2026-10-29 from a fresh kernel with `PYTHONUTF8=1`, using `jupyter nbconvert --execute` without `--allow-errors`.

**All 9 remaining notebooks run top to bottom without errors.**

| Notebook | Status | Cut (prototype / old code removed) |
|---|---|---|
| `figure_1_behaviour_simba` | clean | cells 28-42 removed: "Organization" scratch, the "Bonus" movement vs angular velocity cells, per-rat time series, per-snip cross-correlation. Figures dropped: `fig1_deplete45_perrat_timeseries`, `fig1_deplete45_perrat_angvel_thresholds` |
| `figure_2_photometry` | clean | none |
| `figure_3_crosscorrelation` | clean | none |
| `figure_4_realignment` | clean | cells 14-68 removed (everything below "Old stuff below here"): old full-epoch and single-parameter plots, normalised realignment, old bootstrap blocks, Figure 4D shuffle control (`fig4d_realignment_shuffle_control_*`), `fig4b_auc_before_after_transition`, velocity-based transition detection, and the export of `results/transition_aligned_auc_summary.csv` |
| `figure_5_plotted_by_NaCl` | clean | none |
| `figure_S1_behaviour` | clean | cells 28-40 removed (same "Organization" / "Bonus" block as figure 1) |
| `figure_S2_and_S3_ratbyrat_heatmaps` | clean | none |
| `figure_S4_individual_rats` | clean | cells 5-9 removed: talk figure (`figSx_three_rep_rats`) and the scratch cells (`rats[0,2,8]`, `figure_S3a-e_individual_rats`). Cell 3 ("from realignment notebook, maybe useful for normalization") is still there; it duplicates `snips_to_plot` from cells 2 and 4 and could go too |
| `figure_Sx_individual_transitions` | clean | none (keep or remove is still open) |

The pickle-vs-parquet control run showed that the earlier errors had nothing to do with the data format change.

## Open items from the cuts
- `results/transition_aligned_auc_summary.csv` is no longer written by any notebook. Check whether it is a supplementary table of the paper; if so it ships as-is from `results/` or needs a notebook cell again.
- `src/figure_plotting.py` has plotting functions that only the cut cells used (e.g. `make_realignment_panel_behav_and_da`, `make_normalized_realignment_plot`, `make_realignment_panel_only_one_parameter`). Remove them when the code is cleaned up later; no change now.
- Saved outputs in the notebooks are from the old runs. Re-running in place (phase 3) will refresh them, drop the personal paths, and write new panel PDFs in `paper/figs/panels/`.
