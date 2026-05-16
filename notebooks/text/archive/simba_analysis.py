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
import dill

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.ndimage import gaussian_filter1d

import seaborn as sns
import umap

from figure_config import COLORS 

from trompy import save_figure_atomic

# %%
DATAFOLDER = Path("../data")
PREDFOLDER = Path("D:/TestData/bazzino/simba_preds/Output_all_animals_appetitive")
# PREDFOLDER = Path("D:/TestData/bazzino/simba_preds/quick_test")

ttl_file = DATAFOLDER / "ttls.csv"
# file_with_preds = DATAFOLDER / "PB_NAapp-221024_PB62-221024-120029_Cam1.csv"

# %%
def get_snips(data, ttls, stub):
    
    snips = []
    for inf in ttls[stub].values:
        frame_on = int(inf * 10)
        
        snips.append(prob_app[frame_on-50:frame_on+150].values)

    snips = np.array(snips[:-1])
    
    return snips

def get_shifted_snips(data, ttls, stub, shifted_frames=300):
    shifted_means = []

    shifted=prob_app.copy()
    for x in range(1000):
        shifted = np.roll(shifted, shifted_frames)

        snips = []
        for inf in ttls[stub].values:
            frame_on = int(inf * 10)
            
            snips.append(shifted[frame_on-50:frame_on+150])

        snips = np.array(snips[:-1])
        shifted_means.append(np.mean(snips, axis=0))

    shifted_means = np.array(shifted_means)
    
    return shifted_means


# %%
ttls = pd.read_csv(ttl_file)

for file_with_preds in PREDFOLDER.glob("*.csv"):
    stub = file_with_preds.name.split("_")[2]
    print(stub)
    data = pd.read_csv(file_with_preds)
    prob_app = data.Probability_Appetitive
    
    try:
        snips = get_snips(prob_app, ttls, stub)
        shifted_snips = get_shifted_snips(prob_app, ttls, stub)
        
        real_snips_baselined = np.subtract(snips, np.mean(shifted_snips))
        binarised_snips = (snips > np.percentile(shifted_snips, 97.5)).astype(float)

        f, ax = plt.subplots(ncols=5, figsize=(9, 2),
                            gridspec_kw={'wspace': 0.5})
        
        sns.heatmap(real_snips_baselined, ax=ax[0], cbar=False)
        sns.heatmap(binarised_snips, ax=ax[1], cbar=False)
        
        ax[2].plot(np.mean(snips, axis=0), label="Real Snips")
        ax[2].plot(np.mean(shifted_snips, axis=0), label="Shifted Snips")
        ax[2].fill_between(range(200), np.percentile(shifted_snips, 2.5), np.percentile(shifted_snips, 97.5), alpha=0.2, color='red', label="95% CI")
        
        app_by_trial = [np.sum(snip[50:150] > np.percentile(shifted_snips, 97.5))/100 for snip in snips]
        avr_by_trial = [np.sum(snip[50:150] < np.percentile(shifted_snips, 2.5))/100 for snip in snips]
        hedonic_index = np.array(app_by_trial) - np.array(avr_by_trial)
        hedonic_index = gaussian_filter1d(hedonic_index, sigma=1.2)
        
        ax[3].plot(app_by_trial, label="Appetitive")
        ax[3].plot(avr_by_trial, label="Aversive")
        ax[3].set_ylabel("Prop. of frames > 99th percentile")
        
        ax[4].plot(hedonic_index)
        ax[4].set_ylabel("Hedonic Index")
        ax[4].set_ylim(-1, 1)
        ax[4].axhline(0, color='grey', linestyle='--')
        
        sns.despine(ax=ax[2])
        sns.despine(ax=ax[3])
        sns.despine(ax=ax[4])
    except Exception as e:
        print(f"Error occurred while processing {file_with_preds.name}: {e}")
    
    

# %%

# %%
ttls = pd.read_csv(ttl_file)
data = pd.read_csv(file_with_preds)
prob_app = data.Probability_Appetitive

# %%
stub = 'PB62-221024-120029'

snips = []
for inf in ttls[stub].values:
    frame_on = int(inf * 10)
    
    snips.append(prob_app[frame_on-50:frame_on+150].values)

snips = np.array(snips[:-1])


# %%
f, ax = plt.subplots(nrows=2, figsize=(10, 5), sharex=True)
sns.heatmap(snips, cmap="viridis", ax=ax[0])

ax[1].plot(np.mean(snips, axis=0))
ax[1].set_xlabel("Time")
ax[1].set_ylabel("Mean Probability")

sns.despine(ax=ax[1])

# %%
# np.random.seed(0)
shifted_means = []

shifted=prob_app.copy()
for x in range(1000):
    shifted = np.roll(shifted, 300)

    snips = []
    for inf in ttls[stub].values:
        frame_on = int(inf * 10)
        
        snips.append(shifted[frame_on-50:frame_on+150])

    snips = np.array(snips[:-1])
    shifted_means.append(np.mean(snips, axis=0))

shifted_means = np.array(shifted_means)

real_snips = []
for inf in ttls[stub].values:
    frame_on = int(inf * 10)
    
    real_snips.append(prob_app[frame_on-50:frame_on+150].values)

real_snips = np.array(real_snips[:-1])

f, ax = plt.subplots(figsize=(5, 3))
ax.plot(np.mean(shifted_means, axis=0), label="Mean of shuffles")
ax.fill_between(range(200), np.percentile(shifted_means, 2.5, axis=0), np.percentile(shifted_means, 97.5, axis=0), alpha=0.3, label="95% CI")

ax.plot(np.mean(real_snips, axis=0), color="red", label="Real probability")

ax.axvline(50, color="black", linestyle="--")
ax.axvline(150, color="black", linestyle="--")
ax.set_xticks([0, 50, 100, 150, 200], labels=["-5s", "0s", "+5s", "+10s", "+15s"])

ax.set_ylabel("Prob of appetitive behaviour")

ax.legend()

sns.despine(ax=ax)


    

# %%
real_snips_baselined = np.subtract(real_snips, np.mean(shifted_means))

plt.plot(np.mean(real_snips_baselined, axis=0))

sns.heatmap(real_snips_baselined, cmap="viridis")

# %%
np.percentile(shifted_means, 97.5)

# %%
plt.plot([np.sum(snip[50:150] > np.percentile(shifted_means, 97.5)) for snip in real_snips])

# %%
import dill

# %%
#make list/array of flattened snips for all animals
# sort by condition then infusion
# use cdist to work out pairwise differences
# plot as heatmap/confusion matrix

# Load assembled data from figure_4 preprocessing
assembled_data_path = DATAFOLDER / "assembled_data.pickle"

with open(assembled_data_path, "rb") as f:
    data = dill.load(f)

x_array = data["x_array"]
snips_movement = data["snips_movement"]
snips_simba = data["snips_simba"]

# %%
snips_flattened = []
aucs_flattened = []
aucs_smoothed = []
group_index = []
for id in x_array.id.unique():
    x_id = x_array.query("id == @id")
    snips_id = snips_simba[x_id.index]
    x_id_newindex = x_id.reset_index(drop=True)
    for cond in x_id.condition.unique():
        inf = cond
        x_id_cond = x_id_newindex.query("condition == @cond")
        snips_id_cond = snips_id[x_id_cond.index, 50:150]
        snips_flattened.append(snips_id_cond.flatten())
        aucs_flattened.append(np.mean(snips_id_cond, axis=1))
        aucs_smoothed.append(gaussian_filter1d(np.mean(snips_id_cond, axis=1), sigma=1.2))
        group_index.append((id, cond, inf))

snips_flattened = np.array(snips_flattened)
aucs_flattened = np.array(aucs_flattened)
aucs_smoothed = np.array(aucs_smoothed)

        # f, ax = plt.subplots(figsize=(2, 2))
        # sns.heatmap(snips_id_cond, ax=ax)

# %%
aucs_smoothed.shape

# %%
sort_idx = sorted(
    range(len(group_index)),
    key=lambda i: (group_index[i][1], group_index[i][2])
)

snips_flattened_sorted = [snips_flattened[i] for i in sort_idx]
aucs_flattened_sorted = [aucs_flattened[i] for i in sort_idx]
aucs_smoothed_sorted = [aucs_smoothed[i] for i in sort_idx]
# hedonic_index_sorted = [hedonic_index_all[i] for i in sort_idx]
group_index_sorted = [group_index[i] for i in sort_idx]

# %%
from scipy.spatial.distance import cdist

metric = "euclidean"
# metric = "correlation"
# metric = "cosine"
# metric = "cityblock"

distances = cdist(snips_flattened_sorted[:], snips_flattened_sorted[:], metric=metric)
distances_auc = cdist(aucs_flattened_sorted[:], aucs_flattened_sorted[:], metric=metric)
# distances_hedonic = cdist(hedonic_index_sorted[:], hedonic_index_sorted[:], metric=metric)
distances_auc_smoothed = cdist(aucs_smoothed_sorted[:], aucs_smoothed_sorted[:], metric=metric)



# %%
f, ax = plt.subplots(ncols=3, figsize=(10,3))

sns.heatmap(distances, ax=ax[0],
            # xticklabels=group_index_sorted,
            # yticklabels=group_index_sorted,
            # vmin=0, vmax=15,
            cmap="viridis"
            )

sns.heatmap(distances_auc, ax=ax[1],
            # vmin=0, vmax=6,
            cmap="viridis"
)

sns.heatmap(distances_auc_smoothed, ax=ax[2],
            # vmin=15, vmax=45, # for cityblock metric
            vmin=2, vmax=8,
            cmap="viridis"
)

ax[2].set_yticks([5, 15, 25, 35], labels=["d.10", "d.45", "r.10", "r.45"])
ax[2].set_xticks([5, 15, 25, 35], labels=["d.10", "d.45", "r.10", "r.45"])

# %%
distances_auc_smoothed.shape

# %%



# %%
distance_matric_to_use = distances_auc_smoothed
cmap = plt.get_cmap("Oranges").reversed()

# heatmap with all animals
f, ax = plt.subplots(figsize=(2, 2))

sns.heatmap(distance_matric_to_use, ax=ax,
            # vmin=15, vmax=45, # for cityblock metric
            vmin=3, vmax=7,
            cmap=cmap,
            cbar=False,
)
ax.set_yticks([])
ax.set_xticks([])

save_figure_atomic(f, "simba_distance_heatmap_all_animals.png", folder=Path("../figures"))

# heatmap with groups averaged
arr = np.asarray(distance_matric_to_use)
n_groups = 4

if arr.shape != (40, 40):
    raise ValueError(f"Expected a 40x40 matrix, got {arr.shape}")

idx_groups = np.array_split(np.arange(arr.shape[0]), n_groups)
distances_auc_smoothed_4x4 = np.array(
    [[arr[np.ix_(r_idx, c_idx)].mean() for c_idx in idx_groups] for r_idx in idx_groups]
 )

distances_auc_smoothed_4x4

f, ax = plt.subplots(figsize=(2, 2))
f_, ax_ = plt.subplots(figsize=(0.2, 2))
sns.heatmap(distances_auc_smoothed_4x4, ax=ax,
            # vmin=15, vmax=45, # for cityblock metric
            vmin=3, vmax=7,
            cmap=cmap,
            linewidths=1,      # thickness of grid lines
            linecolor="white",
            cbar_ax=ax_
)

ax.set_yticks([])
ax.set_xticks([])

ax_.set_yticks([])

save_figure_atomic(f, "simba_distance_heatmap.png", folder=Path("../figures"))
save_figure_atomic(f_, "simba_distance_heatmap_colorbar.png", folder=Path("../figures"))

# %%
from sklearn.manifold import MDS

mds = MDS(
n_components=2,
metric=True, # use metric MDS
dissimilarity="euclidean",
random_state=42,
n_init=4,
max_iter=300
)

X_scaled = (distance_matric_to_use - np.min(distance_matric_to_use)) / (np.max(distance_matric_to_use) - np.min(distance_matric_to_use))
X_2d = mds.fit_transform(X_scaled) # shape: (n_samples, 2)

# %%
plt.figure(figsize=(5, 5))
for i, (id, cond, inf) in enumerate(group_index_sorted[:]):
    if i < 10:
        plt.scatter(X_2d[i, 0], X_2d[i, 1], s=40, label=f"{id}_{cond}", color=COLORS[2])
    elif 10 <= i < 20:
        plt.scatter(X_2d[i, 0], X_2d[i, 1], s=40, label=f"{id}_{cond}", color=COLORS[3])
    elif 20 <= i < 30:
        plt.scatter(X_2d[i, 0], X_2d[i, 1], s=40, label=f"{id}_{cond}", color=COLORS[0])
    elif 30 <= i < 40:
        plt.scatter(X_2d[i, 0], X_2d[i, 1], s=40, label=f"{id}_{cond}", color=COLORS[1])
        
# plt.scatter(X_2d[:, 0], X_2d[:, 1], s=40)
plt.xlabel("MDS1")
plt.ylabel("MDS2")
plt.title("MDS projection to 2D")

# %%
f, ax = plt.subplots(figsize=(2, 2),
                     gridspec_kw={'left': 0.15, "bottom": 0.15})

reordered_colors = [COLORS[2], COLORS[3], COLORS[0], COLORS[1]]

for i, (id, cond, inf) in enumerate(group_index_sorted[:]):
    if i < 10:
        plt.scatter(X_2d[i, 0], X_2d[i, 1], s=40, label=f"{id}_{cond}", edgecolor=COLORS[2], facecolor='none')
    elif 10 <= i < 20:
        plt.scatter(X_2d[i, 0], X_2d[i, 1], s=40, label=f"{id}_{cond}", edgecolor=COLORS[3], facecolor='none')
    elif 20 <= i < 30:
        plt.scatter(X_2d[i, 0], X_2d[i, 1], s=40, label=f"{id}_{cond}", edgecolor=COLORS[0], facecolor='none')
    elif 30 <= i < 40:
        plt.scatter(X_2d[i, 0], X_2d[i, 1], s=40, label=f"{id}_{cond}", edgecolor=COLORS[1], facecolor='none')

# plot average for each of the four groups
group_means = []
for group in range(4):
    group_indices = range(group*10, (group+1)*10)
    group_mean = X_2d[group_indices].mean(axis=0)
    group_means.append(group_mean)
    ax.scatter(group_mean[0], group_mean[1], s=100, label=f"Group {group+1}", color=reordered_colors[group])

# plot line connecting each dot to its group mean
for i in range(40):
    group = i // 10
    ax.plot([X_2d[i, 0], group_means[group][0]], [X_2d[i, 1], group_means[group][1]], color=reordered_colors[group], alpha=0.3, zorder=0)
    
ax.set_xlabel("MDS1")
ax.set_ylabel("MDS2")

ax.set_xticks([])
ax.set_yticks([])

sns.despine(ax=ax, offset=5)

save_figure_atomic(f, "simba_distance_mds.png", folder=Path("../figures"))

# %%

reducer = umap.UMAP(n_components=2,
                    metric="euclidean",
                    random_state=42
                    )
X_2d = reducer.fit_transform(aucs_smoothed_sorted)

# %%
plt.figure(figsize=(5, 5))
for i, (id, cond, inf) in enumerate(group_index_sorted[:]):
    if i < 10:
        plt.scatter(X_2d[i, 0], X_2d[i, 1], s=40, label=f"{id}_{cond}", color=COLORS[2])
    elif 10 <= i < 20:
        plt.scatter(X_2d[i, 0], X_2d[i, 1], s=40, label=f"{id}_{cond}", color=COLORS[3])
    elif 20 <= i < 30:
        plt.scatter(X_2d[i, 0], X_2d[i, 1], s=40, label=f"{id}_{cond}", color=COLORS[0])
    elif 30 <= i < 40:
        plt.scatter(X_2d[i, 0], X_2d[i, 1], s=40, label=f"{id}_{cond}", color=COLORS[1])
        
# plt.scatter(X_2d[:, 0], X_2d[:, 1], s=40)
plt.xlabel("MDS1")
plt.ylabel("MDS2")
plt.title("MDS projection to 2D")

# %%
from scipy.ndimage import gaussian_filter1d
        
hedonic_index_all = []

for id in x_array.id.unique():
    x_id = x_array.query("id == @id")
    for cond in x_id.condition.unique():
        snips = snips_simba[x_id.query("condition == @cond").index]
        
        print(snips.shape)
        
        f, ax = plt.subplots(ncols=5, figsize=(9, 2),
                            gridspec_kw={'wspace': 0.5})
        
        sns.heatmap(snips, ax=ax[0])
        
        ax[1].plot(np.mean(snips, axis=0), label="Real Snips")
        


#         binarised_snips = (snips > np.percentile(shifted_snips, 97.5)).astype(float)
#         sns.heatmap(binarised_snips, ax=ax[1], cbar=False)
        z_lims = (-0.2, 0.2)
        app_by_trial = [np.sum(snip[50:150] > z_lims[1])/100 for snip in snips]
        avr_by_trial = [np.sum(snip[50:150] < z_lims[0])/100 for snip in snips]
        hedonic_index = np.array(app_by_trial) - np.array(avr_by_trial)
        


        smoothed_hedonic = gaussian_filter1d(hedonic_index, sigma=1.2)
        hedonic_index_all.append(smoothed_hedonic)
        
        ax[2].plot(app_by_trial, label="Appetitive")
        ax[2].plot(avr_by_trial, label="Aversive")
        ax[2].set_ylabel("Prop. of frames > 99th percentile")
        
        ax[3].plot(smoothed_hedonic)
        ax[3].set_ylabel("Hedonic Index")
        ax[3].set_ylim(-1, 1)
        ax[3].axhline(0, color='grey', linestyle='--')
        
        ax[4].plot([np.mean(snip[50:150]) for snip in snips])
        ax[4].set_ylabel("Mean Z-score")
        ax[4].set_ylim(-2, 2)
        
        sns.despine(ax=ax[1])        
        sns.despine(ax=ax[2])
        sns.despine(ax=ax[3])
        sns.despine(ax=ax[4])

#     except Exception as e:
#         print(f"Error occurred while processing {file_with_preds.name}: {e}")
    
    

# %%
snips_simba.mean()

# %%
