# Notebook run status (top-to-bottom, fresh kernel)

Run on 2026-10-29 after the `data_io` migration, with `PYTHONUTF8=1`.
Control: the same notebooks with the old pickle loader gave identical error lists for `figure_1`, `figure_4` and `figure_S1`, so none of these errors come from the data format change.

| Notebook | Result | Notes |
|---|---|---|
| `figure_1_behaviour_simba` | 3 errors (cells 33, 35, 42, `NameError`) | `rats` is never defined before the per-rat time-series cells. The saved outputs came from an older run state. |
| `figure_2_photometry` | clean | |
| `figure_3_crosscorrelation` | clean | |
| `figure_4_realignment` | 29 errors, first at cell 15 (`KeyError: 'k_da_orig'`), the rest follow from it | The notebook has the `k_da_orig` line commented out, while `src/figure_plotting.py` (line 1237) still reads it. One of the two is out of date. Needs a decision on which is right before release. |
| `figure_5_plotted_by_NaCl` | clean | |
| `figure_S1_behaviour` | 6 errors (cells 33, 35-39, `NameError`) | Same kind of problem as `figure_1`; probably the end-of-notebook prototype cells. |
| `figure_S2_and_S3_ratbyrat_heatmaps` | clean | |
| `figure_S4_individual_rats` | 3 errors (cells 6-8) | Prototype cells: `rats[0,2,8]` (error also in the original saved output) and cells using a column `auc_simba` that is not in the data. |
| `figure_Sx_individual_transitions` | clean | |

## Consequences
- Release blocker: every notebook must run top-to-bottom without errors from a clean clone (verification phase). The 4 notebooks above need the prototype or stale cells fixed or cut. Your note about prototype code at the end of some notebooks probably covers `figure_1`, `figure_S1` and `figure_S4`; `figure_4` looks like a real mismatch between notebook and `src/`.
- `figure_Sx_individual_transitions` is clean but is a "Sx" notebook; keep or remove is still open.
- The saved outputs (embedded figures) in the old notebooks were produced by older code, so re-running to refresh outputs may change figures slightly. Compare against the paper figures when the notebooks are fixed.
