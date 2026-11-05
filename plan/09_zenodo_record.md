# Zenodo record: draft description and checklist

Record DOI: 10.5281/zenodo.21394487 (reserved, not published). Staging folder: `C:\Users\jmc010\Data\bazzino\zenodo_staging`.

## Files to upload (27 in staging, 23 archives go to Zenodo; total ~33.7 GB, limit 50 GB)

| File | Size | Contents |
|---|---|---|
| `processed_data.zip` | 21 MB | `data/assembled/`, `data/shap_values/`, file/subject keys, `ttls.csv` |
| `tanks_PB23.zip` ... `tanks_PB75.zip` (20 files) | 1.4-1.8 GB each, 32.8 GB total | TDT tanks and videos, one zip per rat |
| `dlc.zip` | 180 MB | 40 DeepLabCut CSVs |
| `simba.zip` | 675 MB | 40 SimBA prediction CSVs |

Upload as well (optional, helpful for checking downloads): `tanks_manifest.json`, `dlc_simba_manifest.json`, `processed_data_manifest.json` (MD5s; Zenodo shows its own checksums, so these are optional).
Do NOT upload `missing_videos.csv` (working file).

## Description (paste into the Zenodo form)

**Title:** Data for "Rapidly decreasing phasic dopamine responses precede a valence switch during satiation of sodium appetite"

Data accompanying the paper by Bazzino, Russchen, Keinath, Roitman and McCutcheon (bioRxiv, https://doi.org/10.64898/2026.07.16.738888). Code to reproduce all figures: https://github.com/mccutcheonlab/bazzino-paper.

The record contains two sets of files. You only need the first to reproduce the figures.

**1. Processed data (`processed_data.zip`, 21 MB).** The assembled dataset used by every figure notebook (`data/assembled/`: parquet, npz and json files), SHAP outputs (`data/shap_values/`), and the session/subject key files and TTL times (`data/*.csv`). Extract it in the root of the code repository so that the files end up in `data/`.

**2. Raw data (about 34 GB).** Everything needed to rebuild the processed dataset with `src/assemble_all_data.py`:
- `tanks_<rat>.zip` (20 files): fibre photometry recordings (TDT tank files) for the sodium-depleted and sodium-replete sessions of each rat, plus the session video (`*_Cam1.avi`). The videos are not used by the analysis code.
- `dlc.zip`: DeepLabCut pose-estimation output (one CSV per session, model "bazzino-dlcNov9shuffle4").
- `simba.zip`: SimBA behaviour classifier output (one CSV per session, column `Probability_Appetitive`).

Only the 40 sessions used in the paper are included (20 rats x sodium-depleted and sodium-replete). Other recorded sessions (sodium-replete experienced, thirsty) were not analysed and are not part of this record.

**Download and use.** With the code repository cloned, `python scripts/download_data.py --only processed` fetches and extracts the processed data; `--only raw` fetches the raw data into `raw_data/` (tanks, dlc, simba); no option fetches both. Archives are checked against the Zenodo MD5 checksums. See the repository README.

**Data layout after extraction.**
```
data/assembled/ ...                 processed dataset
data/shap_values/shap_run_1..5/     SHAP outputs
data/{10,45}NaCl_{File,Subject}Key.csv, data/ttls.csv
raw_data/tanks/<rat>-<date>-<time>/ TDT files (.tev .tsq .Tbk .Tdx .tin .tnt) and *_Cam1.avi
raw_data/dlc/*.csv
raw_data/simba/*.csv
```

## Decisions and details to check before publishing

- **Data licence:** Zenodo asks for one. CC BY 4.0 is the usual choice for data (the code repo is GPL-3.0). Needs your decision.
- **Zenodo version vs code:** the record DOI is already in README and CITATION.cff. Once the GitHub repo exists and is released, link it as a related identifier ("is supplemented by" the GitHub repo URL), or use the Zenodo-GitHub integration to archive the code in the same concept record.
- **Related identifier:** add the bioRxiv DOI 10.64898/2026.07.16.738888 as "is supplement to".
- **Rats and sessions:** the description says 20 rats and 40 sessions (verified from the keys: 20 + 20 rows, each rat has one Depleted and one Replete session).
- **Video content:** I have not looked at the videos. Check that sharing them is fine (ethics/data policy, nothing identifying in frame).
- **Authors and affiliations** on the Zenodo form should match `CITATION.cff` (an encoding error in the Tromsø affiliations in that file was fixed on 2026-11-05).
- **The 4 CSVs in `D:\Missing_files_copy`** (PB39 220627, PB23 220601, PB62 221022 and 221024) were not included; they are not needed.
