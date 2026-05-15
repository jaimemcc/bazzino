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

savefigs = False

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
x_array

# %%
# order for data is rep_10, rep_45, dep_10, dep_45
 
colors = ["#67AFD2", "#016895", "#F4795B", "#C74632"]

f, ax = plt.subplots(ncols=len(x_array.cluster_photo.unique()), figsize=(1.8, 1.2),
                     gridspec_kw={'left': 0.05, 'right': 0.95}
                     )

all_clusters_observed_counts = []

for cluster in np.arange(0,2):
    tmp = x_array.query("cluster_photo == @cluster")# .groupby(["condition", "infusiontype"]).count().cluster
    total = len(tmp)
    
    observed_counts = [
        len(tmp.query("condition == 'replete' & infusiontype == '10NaCl'")),
        len(tmp.query("condition == 'replete' & infusiontype == '45NaCl'")),
        len(tmp.query("condition == 'deplete' & infusiontype == '10NaCl'")),
        len(tmp.query("condition == 'deplete' & infusiontype == '45NaCl'"))
    ]
    all_clusters_observed_counts.append(observed_counts) # Add to list for overall test
    
    pie_props = [count / total if total > 0 else 0 for count in observed_counts]
    
    if total > 0:
        expected_counts = [total / 4, total / 4, total / 4, total / 4]
        
        # Perform Chi-squared goodness-of-fit test
        chi2_stat, p_value = stats.chisquare(f_obs=observed_counts, f_exp=expected_counts)
        print(f"Cluster {cluster+1}: n={total}, Chi-squared p-value: {p_value:.4f}")
        
        # You can add an asterisk or other marker to the plot if p_value is significant
        significance_marker = "*" if p_value < 0.05 else ""
    else:
        print(f"Cluster {cluster+1}: n={total}, skipping Chi-squared test (empty cluster)")
        significance_marker = ""
        p_value = 1.0 # or np.nan
    
    ax[cluster].pie(pie_props,
              colors=colors,
              explode=(0.1, 0.1, 0.1, 0.1),
              #autopct='%1.1f%%'
              )
    
    ax[cluster].text(0,-1.7, f"n={total}", ha="center", va="center", fontsize=10, color="k")

if savefigs:
    f.savefig(FIGSFOLDER / "pies_clusters_all.pdf", dpi=600, transparent=True)


# %%
# Perform overall Chi-squared test for independence across all clusters
contingency_table = np.array(all_clusters_observed_counts)
np.savetxt(RESULTSFOLDER / "cluster_contingency_table.csv", contingency_table, delimiter=",", fmt='%d') # fmt='%d' for integers


# Check if the contingency table has valid data for the test
if contingency_table.sum() > 0 and not (np.any(contingency_table.sum(axis=1) == 0) or np.any(contingency_table.sum(axis=0) == 0)):
    chi2_overall, p_overall, dof_overall, expected_freq_overall = stats.chi2_contingency(contingency_table)
    print("\nOverall Chi-squared Test for Independence (Clusters vs. Conditions):")
    print(f"Chi2 Statistic: {chi2_overall:.4f}")
    print(f"P-value: {p_overall:.4f}")
    print(f"Degrees of Freedom: {dof_overall}")
    # print("Expected Frequencies Table:\n", expected_freq_overall)
else:
    print("\nOverall Chi-squared Test for Independence could not be performed due to zero sums in rows/columns or empty table.")

# %%
len(x_array)

# %%
from matplotlib.patches import Rectangle, Patch

f, ax = plt.subplots(figsize=(0.35, 0.45), gridspec_kw={'left': 0.01, 'right': 0.99}
                    )
# for ref: groups = ["10NaCl replete", "45NaCl replete", "10NaCl deplete", "45NaCl deplete"]

groups = ["", ""]

legend_patches = []
for color, text in zip(colors[:2], groups):
    patch = Patch(color=color, label=text)  # Create a patch for the legend
    legend_patches.append(patch)

ax.legend(handles=legend_patches, loc='center', fontsize=8, frameon=False,
          # bbox_to_anchor=(0.5, 0.5),
          # handletextpad=0.5,
          labelspacing=0.2,
          borderpad=0.05)
ax.axis('off')

if savefigs:
    f.savefig(FIGSFOLDER / "legend_pies_replete.pdf", dpi=600, transparent=True)

f, ax = plt.subplots(figsize=(0.35, 0.25), gridspec_kw={'left': 0.01, 'right': 0.99}                     )

legend_patches = []
for color, text in zip(colors[2:], groups):
    patch = Patch(color=color, label=text)  # Create a patch for the legend
    legend_patches.append(patch)

ax.legend(handles=legend_patches, loc='center', fontsize=8, frameon=False,
          # bbox_to_anchor=(0.5, 0.5),
          # handletextpad=0.5,
          labelspacing=0.2,
          borderpad=0.05)
ax.axis('off')

if savefigs:
    f.savefig(FIGSFOLDER / "legend_pies_deplete.pdf", dpi=600, transparent=True)

# %%
# plot heatmap of all trials sorted first by cluster, then by strength of response using
# sortresponse = np.argsort(np.mean(temp[:,sortwindow[0]:sortwindow[1]], axis=1))[::-1]
# or similar

# savefigs=True

list_of_clustered_snips = []
n_of_clusters = []

for cluster in np.arange(0,2):
    print("Cluster", cluster)
    tmp = snips_photo[x_array.cluster_photo == cluster, :]
    sort_order = np.argsort(np.mean(tmp[:,50:150], axis=1))[::-1]
    snips_cluster = tmp[sort_order, :]

    print(snips_cluster.shape)

    snips_cluster.shape
    list_of_clustered_snips.append(snips_cluster)
    n_of_clusters.append(snips_cluster.shape[0])
    
clustered_snips = np.vstack(list_of_clustered_snips)

layout = [["ax", "cbar_ax"],
          ["ax", "empty"]]

f, ax = plt.subplot_mosaic(layout, figsize=(2.3, 3.6), gridspec_kw={'left': 0.1, 'right': 0.9, 'width_ratios': [10, 1], 'wspace': 0.2})  # Adjust width_ratios for colorbar width

cbar_ax = ax["cbar_ax"]
ax["empty"].remove()

ax = ax["ax"]

sns.heatmap(clustered_snips, ax=ax, vmin=-1, vmax=1, cmap="coolwarm", cbar_ax=cbar_ax)
ax.set_yticks([])
ax.set_xticks([])

for c in np.cumsum(n_of_clusters)[:-1]:
    ax.axhline(c, color="white")


ax.plot([150, 198], [2000, 2000], color="k", linewidth=1.5, clip_on=False)
ax.text(174, 2050, "5 s", ha="center", va="top")

ax.plot([50, 150], [-30, -30], color="k", alpha=0.5, linewidth=1.5, clip_on=False)

cbar_ax.set_yticks([])

if savefigs:
    # f.savefig(FIGSFOLDER / "clusters_heatmap_all.pdf", transparent=True)
    f.savefig(FIGSFOLDER / "clusters_heatmap_all.tif",  dpi=300, transparent=True)
    

# %%
len(x_array)

# %%
df = x_array

def get_prop_of_cluster(df, cluster, condition, infusiontype):
    df_temp = (df
               .query("condition == @condition & infusiontype == @infusiontype")
               )
               
    prop_by_trial = []

    for trial in df_temp.trial.unique():
        n = len(df_temp.query("trial == @trial"))
        n_cluster = len(df_temp.query("trial == @trial & cluster_photo == @cluster"))
        prop_by_trial.append(n_cluster / n)
    
    return prop_by_trial

p = get_prop_of_cluster(df, 0, "replete", "45NaCl")

# %%
# make line plot with shaded error of cluster 1 and 2

colors = ["#6FA287", "blue"]
colors = ["black", "grey"]


f, [ax, empty] = plt.subplots(ncols=2, figsize=(2.3, 2.3), gridspec_kw={'left': 0.1, 'right': 0.9, 'width_ratios': [10, 1], 'wspace': 0.2})  # Adjust width_ratios for colorbar width

empty.remove()

for cluster, color in zip(np.arange(0,2), colors):
    print("Cluster", cluster)
    snips_cluster = snips_photo[x_array.cluster_photo == cluster, :]
    x = np.arange(snips_cluster.shape[1]) / 10
    mean = np.mean(snips_cluster, axis=0)
    sd = np.std(snips_cluster, axis=0)
    sem = sd / np.sqrt(snips_cluster.shape[0])
    ci = sem * 1.96
    
    ax.plot(x, np.mean(snips_cluster, axis=0), color=color, lw=1.5)
    ax.fill_between(x, mean-ci, mean+ci, alpha=0.1, color=color)

ax.set_xticks([])
ax.set_yticks([])
ax.set_xlim(0, 20)
ax.axvline(5, color="k", linestyle="--", alpha=0.3)
ax.axvline(15, color="k", linestyle="--", alpha=0.3)
sns.despine(ax=ax, top=True, right=True, left=True, bottom=True)
    

ax.plot([20, 20], [1, 2], color="k")
ax.text(21, 1.5, "1 Z", ha="left", va="center", fontsize=10)

ax.plot([15, 20], [-0.8, -0.8], color="k")
ax.text(17.5, -0.9, "5 s", ha="center", va="top", fontsize=10)

if savefigs:
    f.savefig(FIGSFOLDER / "snips_2clusters_all.pdf", dpi=600, transparent=True)



# %%
vlim = 2

f, ax = plt.subplots(ncols=2, figsize=(2, 1.4))

id = "PB39"
replete_data = snips_photo[(x_array.id == id) & (x_array.condition == "replete"), :] 
deplete_data = snips_photo[(x_array.id == id) & (x_array.condition == "deplete"), :]

sns.heatmap(data=replete_data, vmin=-vlim, vmax=vlim, ax=ax[0], cmap="coolwarm", cbar=False)
sns.heatmap(data=deplete_data, vmin=-vlim, vmax=vlim, ax=ax[1], cmap="coolwarm", cbar=False)

for axis in ax:
    axis.set_yticks([])
    axis.set_xticks([])
    axis.plot([50,150], [-2, -2], color="k", alpha=0.5, linewidth=1.5, clip_on=False)

if savefigs:
    f.savefig(FIGSFOLDER / "rep_heatmaps_10_all.pdf", transparent=True)




# %%
vlim = 2

f, ax = plt.subplots(ncols=2, figsize=(2, 1.4))
f2, cbar_ax = plt.subplots(figsize=(0.2, 1.4))

id = "PB73"
replete_data = snips_photo[(x_array.id == id) & (x_array.condition == "replete"), :] 
deplete_data = snips_photo[(x_array.id == id) & (x_array.condition == "deplete"), :]

sns.heatmap(data=replete_data, vmin=-vlim, vmax=vlim, ax=ax[0], cmap="coolwarm", cbar=True, cbar_ax=cbar_ax)
sns.heatmap(data=deplete_data, vmin=-vlim, vmax=vlim, ax=ax[1], cmap="coolwarm", cbar=False)

for axis in ax:
    axis.set_yticks([])
    axis.set_xticks([])
    axis.plot([50,150], [-2, -2], color="k", alpha=0.5, linewidth=1.5, clip_on=False)

cbar_ax.set_yticks([])

if savefigs:
    f.savefig(FIGSFOLDER / "rep_heatmaps_45_all.pdf", transparent=True)
    f2.savefig(FIGSFOLDER / "rep_heatmaps_45_cbar_all.pdf", transparent=True)

# %%
df = x_array

f, ax = plt.subplots(ncols=2, figsize=(5, 2), sharey=True,
                     gridspec_kw={'wspace': 0.3, "bottom": 0.2}
                     )

cluster=0
colors = ["#67AFD2", "#016895", "#F4795B", "#C74632"]

infusiontype = "10NaCl"
for condition, color in zip(["replete", "deplete"], [colors[0], colors[2]]):
        prop_by_trial = get_prop_of_cluster(df, cluster, condition, infusiontype)
        df_to_plot = pd.DataFrame({"prop": prop_by_trial})
        df_to_plot.rolling = df_to_plot.prop.rolling(window=3).mean()
        ax[0].plot(df_to_plot
                   .rolling
                   , label=condition,
                   marker="o",
                   markersize=5,
                   color=color)
        
infusiontype = "45NaCl"
for condition, color in zip(["replete", "deplete"], [colors[1], colors[3]]):
        prop_by_trial = get_prop_of_cluster(df, cluster, condition, infusiontype)
        df_to_plot = pd.DataFrame({"prop": prop_by_trial})
        df_to_plot.rolling = df_to_plot.prop.rolling(window=3).mean()
        ax[1].plot(df_to_plot
                   .rolling
                   , label=condition,
                   marker="o",
                   markersize=5,
                   color=color)


ax[0].set_ylabel("Proportion of trials in Cluster 1")

for axis in ax:
        axis.set_xlabel("Trial number")
        sns.despine(ax=axis)
        axis.set_xticks(np.arange(0, 51, 10),
                        # labels=["0", "8", "16", "24", "32", "40"]
                        )
        axis.set_yticks([0, 0.5, 1])
        axis.legend(frameon=False)
        
if savefigs:
        f.savefig(FIGSFOLDER / "cluster_1_prop_all.pdf", dpi=600, transparent=True)

# %%
# Question 1 from Mitch/Alex - probability of transitions between clusters

all_transitions = []
for rat in x_array.id.unique():
    for condition in ["replete", "deplete"]:
        x_array_r = x_array.query("id == @rat & condition == @condition")
        transition_matrix = np.zeros((2,2))
        for i in range(len(x_array_r)-1):
            current_cluster = x_array_r.iloc[i].cluster_photo
            next_cluster = x_array_r.iloc[i+1].cluster_photo
            transition_matrix[current_cluster, next_cluster] += 1
        
        all_transitions.append(transition_matrix)

# %%
all_transitions = np.array(all_transitions)
all_transitions_mean = np.mean(all_transitions, axis=0)

all_transitions_mean = all_transitions_mean / np.sum(all_transitions_mean, axis=1, keepdims=True)
all_transitions_mean = np.nan_to_num(all_transitions_mean, nan=0.0)
all_transitions_mean = np.clip(all_transitions_mean, 0, 1)
all_transitions_mean = all_transitions_mean.round(2)

print(all_transitions_mean)
f = plt.figure(figsize=(1.6,1.6))
gs = f.add_gridspec(1, 2, width_ratios=[15, 1],
                    left=0.25, bottom=0.3, right=0.8, top=0.8)  # Adjust width_ratios for colorbar width

ax = f.add_subplot(gs[0])
cbar_ax = f.add_subplot(gs[1])

sns.heatmap(all_transitions_mean,
            cmap="Greys",
            vmin=0, vmax=1,
            annot=True,
            ax=ax, cbar_ax=cbar_ax,
            cbar_kws={"orientation": "vertical"},
            # linewidths=0.5, linecolor="black", annot_kws={"size": 8}, fmt=".2f"
)
ax.set_ylabel("From cluster")
ax.set_xlabel("To cluster")

ax.set_yticks(np.arange(2)+0.5, labels=np.arange(1,3), rotation=0)
ax.set_xticks(np.arange(2)+0.5, labels=np.arange(1,3))

cbar_ax.set_yticks([0, 1])
cbar_ax.set_ylabel("Transition probability", rotation=-90, labelpad=10, fontsize=10)
cbar_ax.tick_params(axis='both', which='both', length=0)

if savefigs:
        f.savefig(FIGSFOLDER / "cluster_transitions_all.pdf", dpi=600, transparent=True)

# %%
all_transitions_sum = np.sum(all_transitions, axis=0)

# %%
all_transitions_sum


# %%
all_transitions_sum[0,1]

# %%
from statsmodels.stats.proportion import proportion_confint
all_transitions_sum = np.sum(all_transitions, axis=0)

# Example: Observed transitions
n_A_to_B = all_transitions_sum[0,1]
n_A_total = np.sum(all_transitions_sum, axis=1)[0]
n_B_to_A = all_transitions_sum[1,0]
n_B_total = np.sum(all_transitions_sum, axis=1)[1]

# Added: Observed transitions for staying in the same state
n_A_to_A = all_transitions_sum[0,0]
n_B_to_B = all_transitions_sum[1,1]

# Confidence intervals
ci_A_to_B = proportion_confint(n_A_to_B, n_A_total, alpha=0.05, method='normal')
ci_B_to_A = proportion_confint(n_B_to_A, n_B_total, alpha=0.05, method='normal')

# Added: Confidence intervals for staying in the same state
ci_A_to_A = proportion_confint(n_A_to_A, n_A_total, alpha=0.05, method='normal')
ci_B_to_B = proportion_confint(n_B_to_B, n_B_total, alpha=0.05, method='normal')

print(f"95% CI for P(A -> A): {ci_A_to_A}, {np.sum(ci_A_to_A) / 2:.4f} +/- {np.diff(ci_A_to_A)[0] / 2:.4f}")
print(f"95% CI for P(A -> B): {ci_A_to_B}, {np.sum(ci_A_to_B) / 2:.4f} +/- {np.diff(ci_A_to_B)[0] / 2:.4f}")
print(f"95% CI for P(B -> A): {ci_B_to_A}, {np.sum(ci_B_to_A) / 2:.4f} +/- {np.diff(ci_B_to_A)[0] / 2:.4f}")
print(f"95% CI for P(B -> B): {ci_B_to_B}, {np.sum(ci_B_to_B) / 2:.4f} +/- {np.diff(ci_B_to_B)[0] / 2:.4f}")


# %%
import statsmodels.api as sm
import statsmodels.formula.api as smf

model = smf.glm("cluster_photo ~ trial * C(condition) * C(infusiontype)",
                              data=x_array,
                              #groups="id", # Specifies the random effects grouping
                              family=sm.families.Binomial())
result = model.fit()
print(result.summary())

# %%
df = x_array.query("condition == 'replete'")

model = smf.glm("cluster_photo ~ trial * C(infusiontype)",
                              data=df,
                              #groups="id", # Specifies the random effects grouping
                              family=sm.families.Binomial())
result = model.fit()
print(result.summary())

# %%
df = x_array.query("condition == 'deplete'")

model = smf.glm("cluster_photo ~ trial * C(infusiontype)",
                              data=df,
                              #groups="id", # Specifies the random effects grouping
                              family=sm.families.Binomial())
result = model.fit()
print(result.summary())

# %% [markdown]
# ### Photo vs behav clusters
# Some stuff on comparison of photometry and behaviour clusters but generally behav clustering not working so well

# %%
import numpy as np
from scipy.stats import pearsonr


real_corr, pval = pearsonr(x_array.cluster_photo, x_array.cluster_vel)
print("Phi coefficient (Pearson r):", real_corr)
print("p-value:", pval)

# Shuffled distribution
corrs = []
for i in range(1000):
    # Shuffle the cluster_vel array
    shuffled_cluster_vel = np.random.permutation(x_array.cluster_vel)
    # Now you can correlate x_array.cluster_photo with shuffled_cluster_vel
    corr, pval = pearsonr(x_array.cluster_photo, shuffled_cluster_vel)
    corrs.append(corr)

# Convert corrs to a numpy array for easier analysis
corrs = np.array(corrs)

# Calculate the mean and standard deviation of the shuffled correlations
mean_corr = np.mean(corrs)
std_corr = np.std(corrs)
# Print the results
print("Mean of shuffled correlations:", mean_corr)
print("Standard deviation of shuffled correlations:", std_corr)

np.sum(corrs < real_corr) / len(corrs)

# %%
df2 = pd.DataFrame(data={"photo": x_array.cluster_photo,
                         "vel": x_array.cluster_vel
})

cluster_corr = np.zeros((2,2))
cluster_corr[0,0] = len(df2.query("photo == 0 & vel == 0"))
cluster_corr[1,0] = len(df2.query("photo == 1 & vel == 0"))
cluster_corr[0,1] = len(df2.query("photo == 0 & vel == 1"))
cluster_corr[1,1] = len(df2.query("photo == 1 & vel == 1"))

# %%
cluster_corr

# %%
from sklearn.metrics import confusion_matrix

cm = confusion_matrix(x_array.cluster_photo, x_array.cluster_vel)
print(cm)

# %%
sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', cbar=False,
             xticklabels=["Cluster 1", "Cluster 2"],
             yticklabels=["Cluster 1", "Cluster 2"])

# %%
