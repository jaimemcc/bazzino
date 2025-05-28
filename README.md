# bazzino
Analysis of data in sodium appetite experiment by Bazzino and Roitman

## Data processing
The analysis uses photometry data and head position/movement data.

Photometry data were assembled by reading in TDT files with assemble_data.ipynb, extracting "snips" of each trial aligned with solenoid onset.
Data were saved as .pickle files that included an accompanying "x" array describing chaacteristics of each trial.

After assembling photometry and DLC data, important to run, equalize_data notebook to ensure trials match up and to save final data .pickle file for all further analysis
