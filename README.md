# bazzino-paper

Code and data for the figures in:

> Bazzino P, Russchen DCR, Keinath AT*, Roitman MF*, McCutcheon JE*. (2026). *Rapidly decreasing phasic dopamine responses precede a valence switch during satiation of sodium appetite*. bioRxiv (preprint; update on acceptance). https://doi.org/10.64898/2026.07.16.738888

Fibre photometry of dopamine and SimBA-classified behaviour during sodium appetite experiments.

- Archived data and code: Zenodo DOI [10.5281/zenodo.21394487](https://doi.org/10.5281/zenodo.21394487)
- Cite this work: see [CITATION.cff](CITATION.cff)

## Quick start

1. Install [pixi](https://pixi.sh), then in this folder run:
   ```
   pixi install
   ```
2. The processed data needed for every figure is in `data/assembled/` (parquet, npz, json). It is also on Zenodo (`pixi run python scripts/download_data.py --only processed`) if you want to refresh or restore it.
3. Open a notebook and run it top to bottom:
   ```
   pixi run jupyter lab
   ```
   Notebooks write panels to `paper/figs/panels/` and statistics tables to `results/`.

Run the tests with `pixi run test`.

## Notebooks to figures

| Notebook | Figure |
|---|---|
| `notebooks/figure_1_behaviour_simba.ipynb` | Fig. 1 |
| `notebooks/figure_2_photometry.ipynb` | Fig. 2 |
| `notebooks/figure_3_crosscorrelation.ipynb` | Fig. 3 |
| `notebooks/figure_4_realignment.ipynb` | Fig. 4 |
| `notebooks/figure_5_plotted_by_NaCl.ipynb` | Fig. 5 |
| `notebooks/figure_S1_behaviour.ipynb` | Fig. S1 |
| `notebooks/figure_S2_and_S3_ratbyrat_heatmaps.ipynb` | Figs. S2, S3 |
| `notebooks/figure_S4_individual_rats.ipynb` | Fig. S4 |
| `notebooks/figure_S5_individual_transitions.ipynb` | Fig. S5 |

Notebook outputs are saved, so figures can be viewed on GitHub without running anything.
Final composed figures are in `paper/figs/`.

## Repository layout

| Path | Contents |
|---|---|
| `notebooks/` | One notebook per figure |
| `src/` | Pipeline and helper code, see [src/README.md](src/README.md) |
| `data/assembled/` | Processed dataset used by the notebooks |
| `data/*.csv`, `data/shap_values/` | File keys, subject keys, TTLs and SHAP outputs |
| `results/` | Statistics tables exported for JASP |
| `paper/figs/` | Final figures |
| `scripts/` | `download_data.py` (fetch data from Zenodo), `prepare_zenodo_*.py` (build the Zenodo archives), `precut_check.py` (pre-release checks) |
| `tests/` | Tests for the data I/O module |

## Rebuilding the processed data from raw data (optional)

Raw recordings (TDT tanks, DeepLabCut CSVs, SimBA predictions) are on Zenodo (DOI [10.5281/zenodo.21394487](https://doi.org/10.5281/zenodo.21394487)).
Download and extract them with `pixi run python scripts/download_data.py --only raw` (about 34 GB, writes to `raw_data/`), then, if you changed the folders, set the paths in `PARAMS` at the top of `src/assemble_all_data.py`
(defaults: `raw_data/tanks`, `raw_data/dlc`, `raw_data/simba`), then run:

```
pixi run python src/assemble_all_data.py
```

This rewrites `data/assembled/` and creates private caches (`data/_cache_*.pickle`, git-ignored).

## License

See [LICENSE](LICENSE).
