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
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
import trompy as tp

from scipy import stats
from scipy.optimize import curve_fit

import dill

rcParams['font.family'] = 'Arial'
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['savefig.transparent'] = True
colors = ["#67AFD2", "#016895", "#F4795B", "#C74632"]
custom_cmap = LinearSegmentedColormap.from_list("custom_diverging", [colors[1], "white", colors[3]])

savefigs = True

DATAFOLDER = Path("..//data")
RESULTSFOLDER = Path("..//results")
FIGSFOLDER = Path("C:/Users/jmc010/Dropbox/Publications in Progress/Bazzino Roitman_sodium/figs")

# %%
with open(DATAFOLDER / "bazzino_data_for_frejus.pickle", "rb") as f:
    data = dill.load(f)

x_array = data["x_array"]
snips_photo = data["snips_photo"]
snips_vel = data["snips_vel"]


# %%
# commented out because loading in auc from frejus file but need to make sure
# that different onset doesn't change anything too much - or can adjust AUC
# calculation

# onset = 70

# auc_snips = snips_photo[:, onset:150].mean(axis=1)
# auc_vel = snips_vel[:, onset:150].mean(axis=1)
# auc_vel.shape, auc_snips.shape

# df = pd.DataFrame({"auc_snips": auc_snips, "auc_vel": auc_vel}).dropna(axis='rows')
# x_array = pd.concat([x_array, df], axis=1)

# %%
# Define a sigmoidal function (logistic function)
def sigmoid(x, L, x0, k, b):
    """
    Sigmoid function.
    L: Maximum value of the curve
    x0: Midpoint of the sigmoid
    k: Steepness of the curve
    b: Baseline offset
    """
    return L / (1 + np.exp(-k * (x - x0))) + b

df2_dep_45 = x_array.query("condition == 'deplete' & infusiontype == '45NaCl'").reset_index(drop=True)

snips_x0, snips_k, snips_r, vel_x0, vel_k, vel_r = [], [], [], [], [], []

for id in df2_dep_45.id.unique():

    tmp = df2_dep_45.query("id == @id")
    f, ax = plt.subplots(figsize=(2.5, 2.5))

    x = tmp.trial.values
    y = tmp.auc_snips.values # can change for tmp.auc_snips.values to fit dopamine signals

    ax.plot(x, y, color=colors[2], alpha=0.5)

    # Fit the sigmoid function to the data
    try:
        popt, pcov = curve_fit(sigmoid, x, y, p0=[max(y), np.median(x), -1, min(y)], maxfev=10000)
        y_fit = sigmoid(x, *popt)
        ax.plot(x, y_fit, color=colors[2], lw=2)
        
        print(f"id: {id}, params: {popt}")
        snips_x0.append(popt[1])  # x0 is the second parameter in popt
        snips_k.append(popt[2])  # k is the third parameter in popt

        y_fit = sigmoid(x, *popt)
        ss_res = np.sum((y - y_fit) ** 2)
        ss_tot = np.sum((y - np.mean(y)) ** 2)
        r_squared = 1 - (ss_res / ss_tot)
        snips_r.append(r_squared)

    except RuntimeError as e:
        print(f"Could not fit sigmoid for id {id}: {e}")
        snips_x0.append(np.nan)  # Append NaN if fitting fails
        snips_k.append(np.nan)
        snips_r.append(np.nan)

    y = tmp.auc_vel.values # can change for tmp.auc_snips.values to fit dopamine signals

    ax.plot(x, y, color=colors[0], alpha=0.5)

    # Fit the sigmoid function to the data
    try:
        popt, pcov = curve_fit(sigmoid, x, y, p0=[max(y), np.median(x), 1, min(y)], maxfev=10000)
        y_fit = sigmoid(x, *popt)
        ax.plot(x, y_fit, color=colors[0], lw=2)
        
        print(f"id: {id}, params: {popt}")
        vel_x0.append(popt[1])  # x0 is the second parameter in popt
        vel_k.append(popt[2])  # k is the third parameter in popt
        
        y_fit = sigmoid(x, *popt)
        ss_res = np.sum((y - y_fit) ** 2)
        ss_tot = np.sum((y - np.mean(y)) ** 2)
        r_squared = 1 - (ss_res / ss_tot)
        vel_r.append(r_squared)

    except RuntimeError as e:
        print(f"Could not fit sigmoid for id {id}: {e}")
        vel_x0.append(np.nan)  # Append NaN if fitting fails
        vel_k.append(np.nan)
        vel_r.append(np.nan)
        
snips_x0 = np.array(snips_x0)
snips_k = np.array(snips_k)
vel_x0 = np.array(vel_x0)
vel_k = np.array(vel_k)


# %%
df_fitted_params = pd.DataFrame({"id": df2_dep_45.id.unique(),
                                 "snips_x0": snips_x0, "snips_k": snips_k, "snips_r": snips_r,
                                 "vel_x0": vel_x0, "vel_k": vel_k, "vel_r": vel_r})

lower, upper = 0, 50
df_fitted_params = df_fitted_params[
    (df_fitted_params['snips_x0'] >= lower) & (df_fitted_params['snips_x0'] <= upper)
    & (df_fitted_params['vel_x0'] >= lower) & (df_fitted_params['vel_x0'] <= upper)
    & (df_fitted_params['snips_r'] >= 0.3) & (df_fitted_params['vel_r'] >= 0.01)
].reset_index(drop=True)

# %%
df_fitted_params

# %%
plt.scatter(df_fitted_params.snips_x0,
            df_fitted_params.vel_x0,
            # c=snips_k,
            cmap=custom_cmap, s=50, edgecolor='k', alpha=0.7)

# %%
snips_x0

# %%
x_id = []
aucs_id = []
aucs_vel_id = []

for i, id in enumerate(df2_dep_45.id.unique()):
    transition = snips_x0[i]

    if np.isnan(transition) == True or transition < 0 or transition > 50:
        continue

    transition = int(transition)
    tmp = df2_dep_45.query("id == @id")
    print(len(tmp))

    x_id.append(tmp.trial.values - transition)
    aucs_id.append(tmp.auc_snips.values)
    aucs_vel_id.append(tmp.auc_vel.values)

    

    

# %%
df_realigned = pd.DataFrame({
    "trial": np.concatenate(x_id),
    # "auc_snips": np.concatenate(aucs_id),
    "auc_vel": np.concatenate(aucs_vel_id)
})

# %%
vel_realigned = df_realigned.groupby("trial").mean()
vel_realigned = vel_realigned
f, ax = plt.subplots(figsize=(2.5, 2.5))
ax.bar(1, np.nanmean(vel_realigned.query("trial < 0").auc_vel.values))
ax.bar(2, np.nanmean(vel_realigned.query("trial > 0").auc_vel.values))

# %%
df_realigned = pd.DataFrame({
    "trial": np.concatenate(x_id),
    "auc_snips": np.concatenate(aucs_id),
    # "auc_vel": np.concatenate(aucs_vel_id)
})

vel_realigned = df_realigned.groupby("trial").mean()
vel_realigned = vel_realigned
f, ax = plt.subplots(figsize=(2.5, 2.5))
ax.bar(1, np.nanmean(vel_realigned.query("trial < 0").auc_snips.values))
ax.bar(2, np.nanmean(vel_realigned.query("trial > 0").auc_snips.values))

# %%
vel_realigned


# %%
def compare_before_and_after(df):
    
    deltas = []
    for i in range(1, 48):
        before_mean = df[:i].mean()
        after_mean = df[i:].mean()
        deltas.append(after_mean - before_mean)

    return deltas

grouped_df = df2_dep_45.groupby("trial").auc_vel.mean()
deltas_snips = compare_before_and_after(grouped_df)

deltas_snips = np.array(deltas_snips)
np.max(np.abs(deltas_snips))

# %%
deltas_snips

# %%
np.nanmean(vel_realigned.auc_vel.iloc[:0].values)

# %%
df_realigned.groupby("trial").mean().plot()

# %%
df2_dep_45.groupby("trial").auc_snips.mean().plot()

# %%
popt, pcov = curve_fit(sigmoid, x, y, p0=[max(y), np.median(x), 1, min(y)])

# %%
realigned_aucs = []

f, ax = plt.subplots(figsize=(2.5, 2.5))
for i, x in enumerate(x_id):
    y = aucs_id[i]
    ax.plot(x, y, color=colors[2], alpha=0.5)
    



# %%
