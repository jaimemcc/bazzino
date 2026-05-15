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

from scipy import stats

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
with open(DATAFOLDER / "bazzino_data.pickle", "rb") as f:
    data = dill.load(f)

x_array = data["x_array"]
snips_data = data["snips_photo"]
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
def get_heatmap_data(snips, x_array, condition, sex, infusiontype):
    if sex == "F" or sex == "M":
        query_string = "condition == @condition & sex == @sex & infusiontype == @infusiontype"
    else:
        query_string = "condition == @condition & infusiontype == @infusiontype"

    heatmap_data = []
    df = x_array.query(query_string)
    # print(df.head())
    for trial in df.trial.unique():
        tmp_snips = snips[x_array.query(query_string + " & trial == @trial").index]
        mean_snip = np.nanmean(tmp_snips, axis=0)
        heatmap_data.append(mean_snip)
        
    return np.array(heatmap_data)


def init_figure():
    f = plt.figure(figsize=(2, 3.5))
    gs = f.add_gridspec(2, 2, hspace=0.1, wspace=0.05, width_ratios=[10, 1])
    
    ax1 = f.add_subplot(gs[0, 0])
    ax2 = f.add_subplot(gs[1, 0])
    cbar_ax = f.add_subplot(gs[0, 1])
    
    return f, ax1, ax2, cbar_ax

def make_heatmap(data, ax, vlim, cbar_ax=None, inf_bar=False):
    
    if cbar_ax is not None:
        cbar=True
    else:
        cbar=False
        
    sns.heatmap(np.array(data),
                cmap="Blues",
                ax=ax,
                cbar=cbar,
                cbar_ax=cbar_ax,
                vmin=0, vmax=vlim
                )
    
    if cbar_ax is not None:
        cbar_ax.set_yticks([])
    
    if inf_bar:
        ax.plot([50,150], [-3, -3], color="black", lw=2, alpha=0.5, clip_on=False)
        
    ax.set_xticks([])
    ax.set_yticks([])

vlim = 5
# savefigs = False

# female, replete
f, ax1, ax2, cbar_ax = init_figure()

heatmap_data = get_heatmap_data(snips_vel, x_vel, "replete", "F", "10NaCl")
f_replete_10_inf = np.nanmean(heatmap_data[:,49:150], axis=1)
make_heatmap(heatmap_data, ax1, vlim, inf_bar=True)

heatmap_data = get_heatmap_data(snips_vel, x_vel, "replete", "F", "45NaCl")
f_replete_45_inf = np.nanmean(heatmap_data[:,49:150], axis=1)
make_heatmap(heatmap_data, ax2, vlim)

cbar_ax.remove()

if savefigs:
    f.savefig(FIGSFOLDER / "heatmap_vel_female_replete.pdf")

# # female, deplete
f, ax1, ax2, cbar_ax = init_figure()

heatmap_data = get_heatmap_data(snips_vel, x_vel, "deplete", "F", "10NaCl")
f_deplete_10_inf = np.nanmean(heatmap_data[:,49:150], axis=1)
make_heatmap(heatmap_data, ax1, vlim, inf_bar=True)

heatmap_data = get_heatmap_data(snips_vel, x_vel, "deplete", "F", "45NaCl")
f_deplete_45_inf = np.nanmean(heatmap_data[:,49:150], axis=1)
make_heatmap(heatmap_data, ax2, vlim)

cbar_ax.remove()

if savefigs:
    f.savefig(FIGSFOLDER / "heatmap_vel_female_deplete.pdf")

# # male, replete
f, ax1, ax2, cbar_ax = init_figure()

heatmap_data = get_heatmap_data(snips_vel, x_vel, "replete", "M", "10NaCl")
m_replete_10_inf = np.nanmean(heatmap_data[:,49:150], axis=1)
make_heatmap(heatmap_data, ax1, vlim, inf_bar=True)

heatmap_data = get_heatmap_data(snips_vel, x_vel, "replete", "M", "45NaCl")
m_replete_45_inf = np.nanmean(heatmap_data[:,49:150], axis=1)
make_heatmap(heatmap_data, ax2, vlim)

cbar_ax.remove()

if savefigs:
    f.savefig(FIGSFOLDER / "heatmap_vel_male_replete.pdf")

# # male, deplete
f, ax1, ax2, cbar_ax = init_figure()

heatmap_data = get_heatmap_data(snips_vel, x_vel, "deplete", "M", "10NaCl")
m_deplete_10_inf = np.nanmean(heatmap_data[:,49:150], axis=1)
make_heatmap(heatmap_data, ax1, vlim, inf_bar=True)

heatmap_data = get_heatmap_data(snips_vel, x_vel, "deplete", "M", "45NaCl")
m_deplete_45_inf = np.nanmean(heatmap_data[:,49:150], axis=1)
make_heatmap(heatmap_data, ax2, vlim, cbar_ax=cbar_ax)

if savefigs:
    f.savefig(FIGSFOLDER / "heatmap_vel_male_deplete.pdf")


# %%
def get_mean_snips(snips, x_array, sex, condition):
    if sex.upper() == "F" or sex.upper() == "M":
        query_string = "condition == @condition & sex == @sex"
    else:
        query_string = "condition == @condition"

    snips_10, snips_45 = [], []
    for id in x_array.query(query_string + " & infusiontype == '10NaCl'").id.unique():
        snips_10.append(np.nanmean(snips[x_array.query(query_string + " & id == @id").index], axis=0))
    for id in x_array.query(query_string + " & infusiontype == '45NaCl'").id.unique():
        snips_45.append(np.nanmean(snips[x_array.query(query_string + " & id == @id").index], axis=0))
        
    return np.array(snips_10), np.array(snips_45)

def init_figure():
    f = plt.figure(figsize=(2, 2))
    gs = f.add_gridspec(1, 2, hspace=0.1, wspace=0.05, width_ratios=[10, 1], bottom=0.2)
    ax = f.add_subplot(gs[0, 0])
    
    return f, ax
    
def plot_snips(snips_10, snips_45, ax, colors_10, colors_45, ylims, scalebar=False):
    for snips, col in zip([snips_10, snips_45], [colors_10, colors_45]):
        x = np.arange(snips.shape[1]) / 10
        mean = np.mean(snips, axis=0)
        sd = np.std(snips, axis=0)
        sem = sd / np.sqrt(snips.shape[0])
        
        ax.plot(x, mean, color=col, lw=1.5)
        ax.fill_between(x, mean-sem, mean+sem, alpha=0.3, color=col)
        
        sns.despine(ax=ax, top=True, right=True, left=True, bottom=True)
        ax.set_xticks([])
        ax.set_yticks([])
        
        ax.set_xlim([0, 20])
        ax.set_ylim(ylims)
        ax.plot([15,20], [ylims[0], ylims[0]], color="black", lw=2, alpha=0.5, clip_on=False)
        ax.text(17.5, ylims[0]-0.1, "5 s", ha="center", va="top", fontsize=10)
        
        if scalebar:
            ax.plot([0, 0], [1, 2], color="black", lw=2, alpha=0.5, clip_on=False)
            
        # print(ax.get_ylim(), ax.get_xlim())
        
    
def get_auc(snips):
    auc = []
    for snip in snips:
        auc.append(np.trapz(snip[50:150]))
    return np.array(auc)

# savefigs = False
ylims = [-5, 10]

# female, replete
f, ax = init_figure()

snips_10, snips_45 = get_mean_snips(snips_vel, x_vel, "F", "replete")
plot_snips(snips_10, snips_45, ax, colors[0], colors[1], ylims)
f_replete_aucs = [get_auc(snips_10), get_auc(snips_45)]

if savefigs:
    f.savefig(FIGSFOLDER / "snips_vel_female_replete.pdf")

# female, deplete
f, ax = init_figure()

snips_10, snips_45 = get_mean_snips(snips_vel, x_vel, "F", "deplete")
plot_snips(snips_10, snips_45, ax, colors[2], colors[3], ylims, scalebar=True)
f_deplete_aucs = [get_auc(snips_10), get_auc(snips_45)]

if savefigs:
    f.savefig(FIGSFOLDER / "snips_vel_female_deplete.pdf")
    
# male, replete
f, ax = init_figure()

snips_10, snips_45 = get_mean_snips(snips_vel, x_vel, "M", "replete")
plot_snips(snips_10, snips_45, ax, colors[0], colors[1], ylims)
m_replete_aucs = [get_auc(snips_10), get_auc(snips_45)]

if savefigs:
    f.savefig(FIGSFOLDER / "snips_vel_male_replete.pdf")
    
# male, deplete
f, ax = init_figure()

snips_10, snips_45 = get_mean_snips(snips_vel, x_vel, "M", "deplete")
plot_snips(snips_10, snips_45, ax, colors[2], colors[3], ylims, scalebar=True)
m_deplete_aucs = [get_auc(snips_10), get_auc(snips_45)]

if savefigs:
    f.savefig(FIGSFOLDER / "snips_vel_male_deplete.pdf")

# %%
## sexes pooled

f, ax1, ax2, cbar_ax = init_figure()

heatmap_data = get_heatmap_data(snips_vel, x_vel, "replete", "both", "10NaCl")
both_replete_10_inf = np.mean(heatmap_data[:,49:150], axis=1)
make_heatmap(heatmap_data, ax1, vlim, inf_bar=True)

heatmap_data = get_heatmap_data(snips_vel, x_vel, "replete", "both", "45NaCl")
both_replete_45_inf = np.mean(heatmap_data[:,49:150], axis=1)
make_heatmap(heatmap_data, ax2, vlim)

cbar_ax.remove()

# if savefigs:
#     f.savefig(FIGSFOLDER / "heatmap_female_replete.pdf", dpi=600, transparent=True)

# female, deplete
f, ax1, ax2, cbar_ax = init_figure()

heatmap_data = get_heatmap_data(snips_vel, x_vel, "deplete", "both", "10NaCl")
both_deplete_10_inf = np.mean(heatmap_data[:,49:150], axis=1)
make_heatmap(heatmap_data, ax1, vlim, inf_bar=True)

heatmap_data = get_heatmap_data(snips_vel, x_vel, "deplete", "both", "45NaCl")
both_deplete_45_inf = np.mean(heatmap_data[:,49:150], axis=1)
make_heatmap(heatmap_data, ax2, vlim)

cbar_ax.remove()

# %%
f, ax = plt.subplots(figsize=(2.2, 2.2),
                     gridspec_kw={"left": 0.25, "right": 0.95, "top": 0.95, "bottom": 0.2})


aucs = [f_replete_aucs, f_deplete_aucs]
barx = [1, 2]
barwidth=0.35

spacer = 0.2
ax.bar(barx[0] - spacer, np.mean(aucs[0][0]), color=colors[0], width=barwidth)
ax.bar(barx[0] + spacer, np.mean(aucs[0][1]), color=colors[1], width=barwidth)

ax.bar(barx[1] - spacer, np.mean(aucs[1][0]), color=colors[2], width=barwidth)
ax.bar(barx[1] + spacer, np.mean(aucs[1][1]), color=colors[3], width=barwidth)

ax.scatter([barx[0] - spacer]*len(aucs[0][0]), aucs[0][0], facecolors="white", edgecolors=colors[0], alpha=0.5)
ax.scatter([barx[0] + spacer]*len(aucs[0][1]), aucs[0][1], facecolors="white", edgecolors=colors[1], alpha=0.5)
ax.scatter([barx[1] - spacer]*len(aucs[1][0]), aucs[1][0], facecolors="white", edgecolors=colors[2], alpha=0.8)
ax.scatter([barx[1] + spacer]*len(aucs[1][1]), aucs[1][1], facecolors="white", edgecolors=colors[3], alpha=0.5)

sns.despine(ax=ax, top=True, right=True, left=False, bottom=True)
ax.set_xticks([])
ax.set_yticks([-0, 50, 100, 150])
ax.set_ylim([-290,200])

ax.set_ylabel("Ang. velocity (AUC)", fontsize=10)

if savefigs:
    f.savefig(FIGSFOLDER / "aucs_vel_female.pdf")

# %%
f, ax = plt.subplots(figsize=(2.2, 2.2),
                     gridspec_kw={"left": 0.25, "right": 0.95, "top": 0.95, "bottom": 0.2})


aucs = [m_replete_aucs, m_deplete_aucs]
barx = [1, 2]
barwidth=0.35

spacer = 0.2
ax.bar(barx[0] - spacer, np.mean(aucs[0][0]), color=colors[0], width=barwidth)
ax.bar(barx[0] + spacer, np.mean(aucs[0][1]), color=colors[1], width=barwidth)

ax.bar(barx[1] - spacer, np.mean(aucs[1][0]), color=colors[2], width=barwidth)
ax.bar(barx[1] + spacer, np.mean(aucs[1][1]), color=colors[3], width=barwidth)

ax.scatter([barx[0] - spacer]*len(aucs[0][0]), aucs[0][0], facecolors="white", edgecolors=colors[0], alpha=0.5)
ax.scatter([barx[0] + spacer]*len(aucs[0][1]), aucs[0][1], facecolors="white", edgecolors=colors[1], alpha=0.5)
ax.scatter([barx[1] - spacer]*len(aucs[1][0]), aucs[1][0], facecolors="white", edgecolors=colors[2], alpha=0.8)
ax.scatter([barx[1] + spacer]*len(aucs[1][1]), aucs[1][1], facecolors="white", edgecolors=colors[3], alpha=0.5)

sns.despine(ax=ax, top=True, right=True, left=False, bottom=True)
ax.set_xticks([])
ax.set_yticks([-0, 50, 100, 150])
ax.set_ylim([-290,200])

ax.set_ylabel("Ang. velocity (AUC)", fontsize=10)

if savefigs:
    f.savefig(FIGSFOLDER / "aucs_vel_male.pdf")

# %%
savefigs = True

def draw_regression_line(y, ax, color):
    x = np.arange(len(y))
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    y_fit = slope * x + intercept
    ax.plot(x, y_fit, color=color, lw=1.5)
    
    print(f"r = {r_value:.2f}, p = {p_value:.3f}")
    
    return r_value, p_value

def make_correlation_plot(inf10, inf45, col10, col45, yaxis=False):
    
    f, ax = plt.subplots(figsize=(1.8,1.8),
                         gridspec_kw={"left": 0.25, "right": 0.95, "top": 0.85, "bottom": 0.22})

    ax.scatter(np.arange(len(inf10)), inf10, color=col10, alpha=0.5)
    ax.scatter(np.arange(len(inf45)), inf45, color=col45, alpha=0.5)

    r, p = draw_regression_line(inf10, ax, col10)
    ax.text(-2, 4, f"0.10 M: r = {r:.2f}, p = {p:.3f}", color=col10, fontsize=8,
            va="bottom", ha="left")
    r, p = draw_regression_line(inf45, ax, col45)
    ax.text(-2, 3.5, f"0.45 M: r = {r:.2f}, p = {p:.3f}", color=col45, fontsize=8,
            va="bottom", ha="left")

    sns.despine(ax=ax)

    ax.set_ylim([0, 3.5])
  
    if yaxis:
        ax.set_yticks([0, 1, 2, 3])
        ax.set_ylabel("Ang. velocity (AUC)")
    else:
        ax.set_yticks([0, 1, 2, 3], labels=["", "", "", ""])

    ax.set_xticks([0, 10, 20, 30, 40, 49], labels=["0", "10", "20", "30", "40", "50"])
    ax.set_xlabel("Trial")
    
    return f
    
f = make_correlation_plot(f_replete_10_inf, f_replete_45_inf, colors[0], colors[1], yaxis=True)
if savefigs:
    f.savefig(FIGSFOLDER / "correlation_vel_f_replete.pdf")
    
f = make_correlation_plot(f_deplete_10_inf, f_deplete_45_inf, colors[2], colors[3])
if savefigs:
    f.savefig(FIGSFOLDER / "correlation_vel_f_deplete.pdf")
    
f = make_correlation_plot(m_replete_10_inf, m_replete_45_inf, colors[0], colors[1], yaxis=True)
if savefigs:
    f.savefig(FIGSFOLDER / "correlation_vel_m_replete.pdf")
    
f = make_correlation_plot(m_deplete_10_inf, m_deplete_45_inf, colors[2], colors[3])
if savefigs:
    f.savefig(FIGSFOLDER / "correlation_vel_m_deplete.pdf")


# %%
with open(DATAFOLDER / "snips_photo_forDLC.pickle", "rb") as f:
    snips = dill.load(f)


# %%
snips.shape

# %%
auc_snips = snips[:, 50:150].mean(axis=1)
auc_vel = snips_vel[:, 50:150].mean(axis=1)
auc_vel.shape, auc_snips.shape

# %%
df = pd.DataFrame({"auc_snips": auc_snips, "auc_vel": auc_vel}).dropna(axis='rows')

f, ax = plt.subplots()
ax.scatter(df.auc_snips, df.auc_vel, alpha=0.1)

x = df.auc_snips
y = df.auc_vel
color = "#D26769"

slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
y_fit = slope * x + intercept
ax.plot(x, y_fit, color=color, lw=2)

ax.set_ylim([-15, 15])

# %%
slope, intercept, r_value, p_value, std_err

# %%
df = pd.DataFrame({"auc_snips": auc_snips, "auc_vel": auc_vel})

# %%
df2 = pd.concat([x_vel, df], axis=1)

# %%
snips_vel.shape

# %%
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

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

df2_dep_45 = df2.query("condition == 'deplete' & infusiontype == '45NaCl'").reset_index(drop=True)

snips_x0 = []
snips_k = []
for id in df2_dep_45.id.unique():

    tmp = df2_dep_45.query("id == @id")
    f, ax = plt.subplots(figsize=(2.5, 2.5))
    ax.plot(tmp.trial, tmp.auc_snips, color=colors[2], alpha=0.5)

    x = tmp.trial.values
    y = tmp.auc_snips.values

    # Fit the sigmoid function to the data
    try:
        popt, pcov = curve_fit(sigmoid, x, y, p0=[max(y), np.median(x), 1, min(y)])
        y_fit = sigmoid(x, *popt)
        ax.plot(x, y_fit, color=colors[2], lw=2)
        
        print(f"id: {id}, params: {popt}")
        snips_x0.append(popt[1])  # x0 is the second parameter in popt
        snips_k.append(popt[2])  # k is the third parameter in popt
    except RuntimeError as e:
        print(f"Could not fit sigmoid for id {id}: {e}")
        snips_x0.append(np.nan)  # Append NaN if fitting fails
        snips_k.append(np.nan)
snips_x0 = np.array(snips_x0)
snips_k = np.array(snips_k)


# %%
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

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

df2_dep_45 = df2.query("condition == 'deplete' & infusiontype == '45NaCl'").reset_index(drop=True)

vel_x0 = []
vel_k = []
for id in df2_dep_45.id.unique():

    tmp = df2_dep_45.query("id == @id")
    f, ax = plt.subplots(figsize=(2.5, 2.5))
    ax.plot(tmp.trial, tmp.auc_vel, color=colors[2], alpha=0.5)

    x = tmp.trial.values
    y = tmp.auc_vel.values

    # Fit the sigmoid function to the data
    try:
        popt, pcov = curve_fit(sigmoid, x, y, p0=[max(y), np.median(x), 1, min(y)])
        y_fit = sigmoid(x, *popt)
        ax.plot(x, y_fit, color=colors[2], lw=2)
        
        print(f"id: {id}, params: {popt}")
        vel_x0.append(popt[1])  # x0 is the second parameter in popt
        vel_k.append(popt[2])  # k is the third parameter in popt
    except RuntimeError as e:
        print(f"Could not fit sigmoid for id {id}: {e}")
        vel_x0.append(np.nan)
        vel_k.append(np.nan)
vel_x0 = np.array(vel_x0)
vel_k = np.array(vel_k)


# %%
vel_x0

# %%
plt.scatter(snips_x0, vel_x0, alpha=0.5, color="black")

# %%
plt.scatter(snips_k, vel_k, alpha=0.5, color="black")

# %%
