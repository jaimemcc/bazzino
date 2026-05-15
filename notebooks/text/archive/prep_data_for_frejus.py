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

# %% [markdown]
# # Preparation of data for Frejus talk.
# This simple script basically adds smoothing to angular velocity data and AUCs so that Frejus notebook can focus on plotting of figures etc. In theory, this could be used for preparing data for all figures etc.
#
# Maybe should be renamed - final_post_processing

# %%
from pathlib import Path
import numpy as np
import pandas as pd

import dill

DATAFOLDER = Path("..//data")

# %%
with open(DATAFOLDER / "bazzino_data_with_clusters.pickle", "rb") as f:
    data = dill.load(f)

x_array = data["x_array"]
snips_photo = data["snips_photo"]
snips_vel = data["snips_vel"]


# %%
def smooth_array(arr, window_size=5):
    """
    Smooth a 2D array along one dimension using a moving average.
    
    :param arr: 2D NumPy array
    :param window_size: Size of the smoothing window
    :return: Smoothed 2D array
    """
    kernel = np.ones(window_size) / window_size
    smoothed = np.apply_along_axis(lambda m: np.convolve(m, kernel, mode='same'), axis=1, arr=arr)
    return smoothed

# Example usage
snips_vel = smooth_array(snips_vel, window_size=5)

# %%
auc_snips = snips_photo[:, 50:150].mean(axis=1)
auc_vel = snips_vel[:, 50:150].mean(axis=1)
auc_vel.shape, auc_snips.shape

df = pd.DataFrame({"auc_snips": auc_snips, "auc_vel": auc_vel}).dropna(axis='rows')
x_array = pd.concat([x_array, df], axis=1)

# %%
with open(DATAFOLDER / "bazzino_data_for_frejus.pickle", "wb") as f:
    dill.dump({
        "x_array": x_array,
        "snips_photo": snips_photo,
        "snips_vel": snips_vel
    }, f)

# %%
