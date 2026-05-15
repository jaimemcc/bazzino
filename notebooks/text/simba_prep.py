# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.2
#   kernelspec:
#     display_name: default
#     language: python
#     name: python3
# ---

# %%
from pathlib import Path

import numpy as np
import pandas as pd

# %%
DATAFOLDER = Path("../data")

ttl_file = DATAFOLDER / "ttls.csv"
csv_file = DATAFOLDER / "PB_NAapp-220629_PB39-220629-105034_Cam1.csv"
file_with_preds = DATAFOLDER / "PB_NAapp-221024_PB62-221024-120029_Cam1.csv"

output_file = DATAFOLDER / "target_files" / csv_file.name


# %%
## Prototype for writing binary target behaviours

# %%

def get_stub(ttls, csv_file):
    for col in ttls.columns:
        if col in csv_file.name:
            print(f"Found {col} in {csv_file.name}")
            return col
    return None

def get_frames(ttls, stub, delay_frames=0, infusion_frames=100):
    ttl_onsets = ttls.loc[:, stub]
    frames_on = [int(onset * 10) + delay_frames for onset in ttl_onsets]
    frames_off = [frame + infusion_frames for frame in frames_on]
    return frames_on, frames_off

def make_vector(frames_on, frames_off, total_frames):
    vector = np.zeros(total_frames)
    for on, off in zip(frames_on, frames_off):
        vector[on:off] = 1
    return vector

def write_output(csv_file, vector, output_file, output_column="target"):
    csv = pd.read_csv(csv_file)
    csv[output_column] = vector
    csv.to_csv(output_file, index=False)

ttls = pd.read_csv(ttl_file)
stub = get_stub(ttls, csv_file)
frames_on, frames_off = get_frames(ttls, stub)

csv = pd.read_csv(csv_file)
vector = make_vector(frames_on, frames_off, len(csv))

write_output(csv_file, vector, output_file, output_column="target")


# %%
dlc_file = DATAFOLDER / "PB_NAapp-220614_PB31-220614-113720_Cam1DLC_Resnet50_bazzino-dlcNov9shuffle4_snapshot_best-200.csv"

dlc_csv = pd.read_csv(dlc_file)

# %%
dlc_csv.head()


# %%
def keep_dlc_bodyparts(
    input_csv,
    bodyparts_to_keep,
    output_csv=None,
    bodypart_level=1,
):
    """Keep only selected bodyparts in a DLC CSV and write back with the same 3-row header layout."""
    input_csv = Path(input_csv)
    output_csv = Path(output_csv) if output_csv is not None else input_csv

    # DLC exports typically use a 3-row column header: scorer/bodyparts/coords.
    df = pd.read_csv(input_csv, header=[0, 1, 2])

    keep = set(bodyparts_to_keep)
    cols_to_keep = []
    for col in df.columns:
        # Keep metadata-like header columns plus requested bodyparts.
        if len(col) > bodypart_level and (col[bodypart_level] in keep or col[bodypart_level] == "bodyparts"):
            cols_to_keep.append(col)

    if not cols_to_keep:
        raise ValueError(f"No columns matched requested bodyparts: {sorted(keep)}")

    filtered = df.loc[:, cols_to_keep]
    filtered.to_csv(output_csv, index=False)

    kept = sorted({c[bodypart_level] for c in cols_to_keep if len(c) > bodypart_level and c[bodypart_level] != "bodyparts"})
    print("Saved filtered DLC CSV")
    print(f"Input:  {input_csv}")
    print(f"Output: {output_csv}")
    print(f"Bodyparts kept ({len(kept)}): {kept}")


# Example: keep only nose + left ear and overwrite the same file.
keep_dlc_bodyparts(dlc_file, ["nose", "l_ear"], output_csv=dlc_file)

# %%
import dill
with open("../data/_cache_simba.pickle", "rb") as f:
    simba = dill.load(f)

# %%
x_array = simba["x_simba"]

# %%
print(len(x_array))
for id in x_array.id.unique():
    tmp = x_array[x_array.id == id]
    if len(tmp) < 98:
        print(id, len(tmp), tmp.condition.unique())

# %%
x_array.id.unique()

# %%
import dill
with open("../data/_cache_behav.pickle", "rb") as f:
    dlc = dill.load(f)

x_dlc = dlc["x_behav"]

# %%
x_dlc.id.unique()

# %%
