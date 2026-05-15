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
import matplotlib.pyplot as plt
from matplotlib import rcParams
import seaborn as sns
import trompy as tp

from scipy.ndimage import gaussian_filter1d

import dill

rcParams['font.family'] = 'Arial'
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['savefig.transparent'] = True
colors = ["#67AFD2", "#016895", "#F4795B", "#C74632"]

savefigs = False

DATAFOLDER = Path("..//data")
RESULTSFOLDER = Path("..//results")
FIGSFOLDER = Path("C:/Users/jmc010/Dropbox/Publications in Progress/Bazzino Roitman_sodium/figs")

# %%
with open(DATAFOLDER / "bazzino_data_with_clusters_and_dists.pickle", "rb") as f:
    data = dill.load(f)
x_array = data["x_array"]
snips_photo = data["snips_photo"]
snips_vel = data["snips_vel"]
pca = data["pca"]

# %%
# get reduced x_array and list of rats
x_array_red = x_array.query("infusiontype == '45NaCl' & condition == 'deplete'")
rats = x_array_red.id.unique()

# %%
# make smoothed clusterness column in case needed
for rat in rats:
    mask = x_array.id == rat
    clusterness = x_array.loc[mask, 'clusterness_photo'].values
    
    # Apply Gaussian smoothing
    smoothed = gaussian_filter1d(clusterness, sigma=2)  # Adjust sigma
    x_array.loc[mask, 'clusterness_photo_smoothed'] = smoothed

# %%

x_rat = x_array_red.query("id == @rat")

# x = x_rat.clusterness_photo.values
x = x_rat.euclidean_diff.values
y = x_rat.auc_vel.values

# Compute cross-correlation
cross_corr = np.correlate(x, y, mode='full')
# Lags
lags = np.arange(-len(x) + 1, len(x))
# Find the lag with the maximum correlation
max_corr_index = np.argmax(cross_corr)
best_lag = lags[max_corr_index]
print("Cross-correlation values:", cross_corr)
print("Best lag:", best_lag)
# Plot the cross-correlation
plt.stem(lags, cross_corr)
plt.title("Cross-Correlation")
plt.xlabel("Lag")
plt.ylabel("Cross-Correlation")
plt.show()

# %%
x_rat.columns

# %%
all_lags = []

for rat in rats:
    x_rat = x_array_red.query("id == @rat")

    x = x_rat.euclidean_diff.values
    # x = x_rat.clusterness_photo.values
    y = x_rat.auc_vel.values

    cross_corr = np.correlate(x, y, mode='full')
    # Lags
    lags = np.arange(-len(x) + 1, len(x))
    # Find the lag with the maximum correlation
    max_corr_index = np.argmax(cross_corr)
    min_corr_index = np.argmin(cross_corr)
    best_lag = lags[min_corr_index]

    all_lags.append(best_lag)

print("All best lags:", all_lags)

# lags are probably due to edge effects, let's think about cross-correlating smoothed data or a better measure than angular veloscity

# %%
x_rat.columns

# %%
all
