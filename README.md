# bazzino
Analysis of data in sodium appetite experiment by Bazzino and Roitman

## Data processing
The analysis uses photometry data and head position/movement data.

Photometry data were assembled by reading in TDT files with assemble_data.ipynb, extracting "snips" of each trial aligned with solenoid onset.
Data were saved as .pickle files that included an accompanying "x" array describing chaacteristics of each trial.

After assembling photometry and DLC data, important to run, equalize_data notebook to ensure trials match up and to save final data .pickle file for all further analysis
## Setup on a New Machine

### Initial Setup
When you first clone this repository on a new machine, run:

```bash
pixi install
```

This installs all dependencies, including `nbstripout` which is used to strip notebook outputs from git commits.

### Configure nbstripout with Git
After the initial pixi install, configure the git filter for nbstripout:

```bash
pixi run nbstripout --install --attributes .gitattributes
```

This one-time setup configures git to automatically strip Jupyter notebook outputs and metadata when committing, preventing merge conflicts caused by cell execution counts and IDs changing locally.

**Why is this needed?**
The git filter configuration is stored in `.git/config` which is not version-controlled. Each machine needs to configure it locally. Without this, you'll see notebooks showing as modified even after syncing with the remote.