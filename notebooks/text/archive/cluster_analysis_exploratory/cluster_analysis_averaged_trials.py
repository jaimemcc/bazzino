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
from matplotlib.patches import Rectangle, Patch
import seaborn as sns
import trompy as tp

import dill

rcParams['font.family'] = 'Arial'
 
colors = ["#67AFD2", "#016895", "#F4795B", "#C74632"]

savefigs = False

DATAFOLDER = Path("..//data")
RESULTSFOLDER = Path("..//results")
FIGSFOLDER = Path("C:/Users/jmc010/Dropbox/Publications in Progress/Bazzino Roitman_sodium/figs")

# %%
with open(DATAFOLDER / "x_array_with_7clusters.pickle", "rb") as f:
    x_array = dill.load(f)

# %%
x_array

# %%
with open(DATAFOLDER / "snips_data_selected_conditions_reduced.pickle", "rb") as f:
    data = dill.load(f)

snips_10NaCl = data["snips_10NaCl_reduced"]
snips_45NaCl = data["snips_45NaCl_reduced"]

snips_all = np.vstack([snips_10NaCl, snips_45NaCl])
snips_all.shape

# %%
f, ax = plt.subplots(ncols=len(x_array.cluster.unique()), sharey=True, figsize=(4.5, 1.5),
                                          gridspec_kw={'left': 0.05, 'right': 0.95, 'top': 0.85}
                                          )

for cluster in np.arange(0,7):
    print("Cluster", cluster)
    snips_cluster = snips_all[x_array.cluster == cluster, :]
    x = np.arange(snips_cluster.shape[1]) / 10
    mean = np.mean(snips_cluster, axis=0)
    sem = np.std(snips_cluster, axis=0) / np.sqrt(snips_cluster.shape[0])
    sd = np.std(snips_cluster, axis=0)
    
    ax[cluster].plot(x, np.mean(snips_cluster, axis=0), color="#6FA287", lw=1)
    ax[cluster].fill_between(x, mean-sd, mean+sd, alpha=0.3, color="#6FA287")
    ax[cluster].set_title(f"{cluster+1}")

for axis in ax:
    axis.set_xticks([])
    axis.set_yticks([])
    axis.axvline(5, color="k", linestyle="--", alpha=0.3)
    axis.axvline(15, color="k", linestyle="--", alpha=0.3)
    sns.despine(ax=axis, top=True, right=True, left=True, bottom=True)
    

ax[0].plot([0, 0], [1, 2], color="k")
ax[0].text(-1, 1.5, "1 Z", ha="right", va="center", fontsize=10)

ax[6].plot([20, 25], [-1.5, -1.5], color="k")
ax[6].text(22.5, -1.7, "5 s", ha="center", va="top", fontsize=10)

if savefigs:
    f.savefig(FIGSFOLDER / "snips_clusters.pdf", dpi=600, transparent=True)

# %%
from scipy import stats
# order for data is rep_10, rep_45, dep_10, dep_45

f, ax = plt.subplots(ncols=len(x_array.cluster.unique()), figsize=(4.5, 0.8),
                     gridspec_kw={'left': 0.05, 'right': 0.95}
                     )

all_clusters_observed_counts = []

for cluster in np.arange(0,7):
    tmp = x_array.query("cluster == @cluster")# .groupby(["condition", "infusiontype"]).count().cluster
    total = len(tmp)
    
    # Observed counts for each category
    observed_counts = [
        len(tmp.query("condition == 'replete' & infusiontype == '10NaCl'")),
        len(tmp.query("condition == 'replete' & infusiontype == '45NaCl'")),
        len(tmp.query("condition == 'deplete' & infusiontype == '10NaCl'")),
        len(tmp.query("condition == 'deplete' & infusiontype == '45NaCl'"))
    ]
    all_clusters_observed_counts.append(observed_counts) # Add to list for overall test
    
    pie_props = [count / total if total > 0 else 0 for count in observed_counts]
    
    # Expected counts for each category (assuming random distribution, so total/4 for each)
    # Ensure total is not zero to avoid division by zero if a cluster is empty
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
    f.savefig(FIGSFOLDER / "pies_clusters.pdf", dpi=600, transparent=True)
    
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

# --- New code for saving in long format for JASP ---
# Define row and column headers for context.
# The number of clusters is derived from the shape of the contingency table.
cluster_labels = [f"Cluster {i+1}" for i in range(contingency_table.shape[0])]
# These condition_type_labels must match the order used when creating 'observed_counts':
# i.e., replete_10NaCl, replete_45NaCl, deplete_10NaCl, deplete_45NaCl
condition_type_labels = ["Replete_10NaCl", "Replete_45NaCl", "Deplete_10NaCl", "Deplete_45NaCl"]

# Convert the NumPy contingency_table to a pandas DataFrame (this is wide format)
# Assign cluster_labels as the index and condition_type_labels as columns.
contingency_df_wide = pd.DataFrame(data=contingency_table, index=cluster_labels, columns=condition_type_labels)

# Name the index so it becomes a column named 'Cluster' after reset_index()
contingency_df_wide.index.name = 'Cluster'
contingency_df_wide = contingency_df_wide.reset_index() # 'Cluster' is now a regular column

# Melt the DataFrame from wide to long format.
# 'Cluster' is the identifier variable.
# The other columns (condition_type_labels) are unpivoted.
contingency_df_long = contingency_df_wide.melt(
    id_vars=['Cluster'],                    # Column(s) to keep as identifiers
    value_vars=condition_type_labels,       # Columns to unpivot
    var_name='Condition_Type',              # Name for the new column holding the original column names
    value_name='Counts'                     # Name for the new column holding the values
)

# Save the long format DataFrame to a CSV file
# index=False prevents pandas from writing its DataFrame index as a column in the CSV.
long_format_csv_filename = RESULTSFOLDER / "contingency_table_long_format_for_jasp.csv"
contingency_df_long.to_csv(long_format_csv_filename, index=False)
print(f"\nContingency table in long format saved to {long_format_csv_filename}")
# --- End of new code ---



# %%
f, ax = plt.subplots(figsize=(0.35, 0.25), gridspec_kw={'left': 0.01, 'right': 0.99}
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

list_of_clustered_snips = []
n_of_clusters = []

for cluster in np.arange(0,7):
    print("Cluster", cluster)
    tmp = snips_all[x_array.cluster == cluster, :]
    sort_order = np.argsort(np.mean(tmp[:,50:150], axis=1))[::-1]
    snips_cluster = tmp[sort_order, :]

    print(snips_cluster.shape)

    snips_cluster.shape
    list_of_clustered_snips.append(snips_cluster)
    n_of_clusters.append(snips_cluster.shape[0])
    
clustered_snips = np.vstack(list_of_clustered_snips)

f = plt.figure(figsize=(1.8, 2.8))
gs = f.add_gridspec(2, 2, width_ratios=[10, 1])  # Adjust width_ratios for colorbar width

# Add heatmap and colorbar axes
ax = f.add_subplot(gs[:, 0])
cbar_ax = f.add_subplot(gs[0, 1])

sns.heatmap(clustered_snips, ax=ax, vmin=-1, vmax=1, cmap="coolwarm", cbar_ax=cbar_ax)
ax.set_yticks([])
ax.set_xticks([])

for c in np.cumsum(n_of_clusters)[:-1]:
    ax.axhline(c, color="white")


ax.plot([150, 198], [485, 485], color="k", linewidth=1.5, clip_on=False)
ax.text(174, 490, "5 s", ha="center", va="top")

cbar_ax.set_yticks([])

if savefigs:
    f.savefig(FIGSFOLDER / "clusters_heatmap.pdf", dpi=600, transparent=True)

# %%
with open(DATAFOLDER / "x_array_with_2clusters.pickle", "rb") as f:
    x_array2 = dill.load(f)

## alternative - construct x_array2 by combining clusters 1,2,3,4

def combine_clusters(df):

    supercluster1 = [1,2,3,4]

    superclusters = []
    for c in df.cluster:
        if c in supercluster1:
            superclusters.append(1)
        else:
            superclusters.append(2)

    return (df
            .assign(cluster=superclusters))

# x_array2 = combine_clusters(x_array)


# %%
df = x_array2

def get_prop_of_cluster(df, cluster, condition, infusiontype):
    df_temp = (df
               .query("condition == @condition & infusiontype == @infusiontype")
               )
               
    prop_by_trial = []

    for trial in df_temp.trial.unique():
        n = len(df_temp.query("trial == @trial"))
        n_cluster = len(df_temp.query("trial == @trial & cluster == @cluster"))
        prop_by_trial.append(n_cluster / n)
    
    return prop_by_trial

p = get_prop_of_cluster(df, 0, "replete", "45NaCl")

# %%

list_of_clustered_snips = []
n_of_clusters = []

for cluster in np.arange(0,2):
    print("Cluster", cluster)
    tmp = snips_all[x_array2.cluster == cluster, :]
    sort_order = np.argsort(np.mean(tmp[:,50:150], axis=1))[::-1]
    snips_cluster = tmp[sort_order, :]

    print(snips_cluster.shape)

    snips_cluster.shape
    list_of_clustered_snips.append(snips_cluster)
    n_of_clusters.append(snips_cluster.shape[0])
    
clustered_snips = np.vstack(list_of_clustered_snips)

f = plt.figure(figsize=(1.8, 2.8))
gs = f.add_gridspec(2, 2, width_ratios=[10, 1])  # Adjust width_ratios for colorbar width

# Add heatmap and colorbar axes
ax = f.add_subplot(gs[:, 0])
cbar_ax = f.add_subplot(gs[0, 1])

sns.heatmap(clustered_snips, ax=ax, vmin=-1, vmax=1, cmap="coolwarm", cbar_ax=cbar_ax)
ax.set_yticks([])
ax.set_xticks([])

for c in np.cumsum(n_of_clusters)[:-1]:
    ax.axhline(c, color="white")


ax.plot([150, 198], [485, 485], color="k", linewidth=1.5, clip_on=False)
ax.text(174, 490, "5 s", ha="center", va="top")

cbar_ax.set_yticks([])

if savefigs:
    f.savefig(FIGSFOLDER / "2clusters_heatmap.pdf", dpi=600, transparent=True)

# %%
# make line plot with shaded error of cluster 1 and 2

snips_all.shape
x_array2.shape

colors = ["#6FA287", "blue"]

f = plt.figure(figsize=(1.8, 1.8))
gs = f.add_gridspec(2, 2, width_ratios=[10, 1])  # Adjust width_ratios for colorbar width

ax = f.add_subplot(gs[:, 0])

for cluster, color in zip(np.arange(0,2), colors):
    print("Cluster", cluster)
    snips_cluster = snips_all[x_array2.cluster == cluster, :]
    x = np.arange(snips_cluster.shape[1]) / 10
    mean = np.mean(snips_cluster, axis=0)
    sem = np.std(snips_cluster, axis=0) / np.sqrt(snips_cluster.shape[0])
    
    ax.plot(x, np.mean(snips_cluster, axis=0), color=color, lw=1)
    ax.fill_between(x, mean-sd, mean+sd, alpha=0.3, color=color)

ax.set_xticks([])
ax.set_yticks([])
ax.axvline(5, color="k", linestyle="--", alpha=0.3)
ax.axvline(15, color="k", linestyle="--", alpha=0.3)
sns.despine(ax=ax, top=True, right=True, left=True, bottom=True)
    

ax.plot([20, 20], [1, 2], color="k")
ax.text(21, 1.5, "1 Z", ha="left", va="center", fontsize=10)

ax.plot([15, 20], [-0.8, -0.8], color="k")
ax.text(17.5, -0.9, "5 s", ha="center", va="top", fontsize=10)

if savefigs:
    f.savefig(FIGSFOLDER / "snips_2clusters.pdf", dpi=600, transparent=True)



# %%
# make rep heatplots of two different rats over time - make for all rats and see who looks best
# try with original data and reduced data

for id in x_array2.id.unique():
    
    print("ID", id)
    for condition in ["replete", "deplete"]:
        tmp = (x_array2
                .query("id == @id & condition == @condition")
                )
        
        snips_to_plot = snips_all[(x_array2.id == id) & (x_array2.condition == condition), :]
        
        f, ax = plt.subplots(figsize=(1.8, 1.8))
        sns.heatmap(snips_to_plot, vmin=-1, vmax=1, cmap="coolwarm", ax=ax)
        plt.title(f"{id} {condition}")

# %%
# possible representatives
# PB73 - 45
# PB39 - 10

# %%
reps = ["PB73, PB39"]
x_array2.query("id == 'PB39'").infusiontype.unique()

# %%
vlim = 2

f, ax = plt.subplots(ncols=2, figsize=(2.5, 1.4))

id = "PB39"
replete_data = snips_all[(x_array2.id == id) & (x_array2.condition == "replete"), :] 
deplete_data = snips_all[(x_array2.id == id) & (x_array2.condition == "deplete"), :]

sns.heatmap(data=replete_data, vmin=-vlim, vmax=vlim, ax=ax[0], cmap="coolwarm", cbar=False)
sns.heatmap(data=deplete_data, vmin=-vlim, vmax=vlim, ax=ax[1], cmap="coolwarm", cbar=False)

for axis in ax:
    axis.set_yticks([])
    axis.set_xticks([])

if savefigs:
    f.savefig(FIGSFOLDER / "rep_heatmaps_10.pdf", transparent=True)




# %%
vlim = 2

f, ax = plt.subplots(ncols=2, figsize=(2.5, 1.4))
f2, cbar_ax = plt.subplots(figsize=(0.2, 1.4))

id = "PB73"
replete_data = snips_all[(x_array2.id == id) & (x_array2.condition == "replete"), :] 
deplete_data = snips_all[(x_array2.id == id) & (x_array2.condition == "deplete"), :]

sns.heatmap(data=replete_data, vmin=-vlim, vmax=vlim, ax=ax[0], cmap="coolwarm", cbar=True, cbar_ax=cbar_ax)
sns.heatmap(data=deplete_data, vmin=-vlim, vmax=vlim, ax=ax[1], cmap="coolwarm", cbar=False)

for axis in ax:
    axis.set_yticks([])
    axis.set_xticks([])

cbar_ax.set_yticks([])

if savefigs:
    # f.savefig(FIGSFOLDER / "rep_heatmaps_45.pdf", transparent=True)
    f2.savefig(FIGSFOLDER / "rep_heatmaps_45_cbar.pdf", transparent=True)

# %%
df = x_array2

f, ax = plt.subplots(ncols=2, figsize=(5, 2.5), sharey=True)

cluster=0
colors = ["#67AFD2", "#016895", "#F4795B", "#C74632"]

infusiontype = "10NaCl"
for condition, color in zip(["replete", "deplete"], [colors[0], colors[2]]):
        prop_by_trial = get_prop_of_cluster(df, cluster, condition, infusiontype)
        df_to_plot = pd.DataFrame({"prop": prop_by_trial})
        df_to_plot.rolling = df_to_plot.prop.rolling(window=3).mean()
        ax[0].plot(df_to_plot
                   #.rolling
                   , label=condition,
                   marker="o", color=color)
        
infusiontype = "45NaCl"
for condition, color in zip(["replete", "deplete"], [colors[1], colors[3]]):
        prop_by_trial = get_prop_of_cluster(df, cluster, condition, infusiontype)
        df_to_plot = pd.DataFrame({"prop": prop_by_trial})
        df_to_plot.rolling = df_to_plot.prop.rolling(window=3).mean()
        ax[1].plot(df_to_plot
                   #.rolling
                   , label=condition,
                   marker="o", color=color)


ax[0].set_ylabel("Proportion of trials in Cluster 1")

for axis in ax:
        axis.set_xlabel("Trial number")
        sns.despine(ax=axis)
        axis.set_xticks(np.arange(0, 12, 2),
                        labels=["0", "8", "16", "24", "32", "40"]
                        )
        axis.legend(frameon=False)
        #axis.set_ylim(0, 1)
        
if savefigs:
        f.savefig(FIGSFOLDER / "cluster_1_prop.pdf", dpi=600, transparent=True)

# %%
# stats for prop cluster 1 over time
x_array2

# %%
import statsmodels.api as sm
import statsmodels.formula.api as smf

model = smf.glm("cluster ~ trial * C(condition) * C(infusiontype)",
                              data=x_array2,
                              #groups="id", # Specifies the random effects grouping
                              family=sm.families.Binomial())
result = model.fit()
print(result.summary())

# %%
# Question 1 from Mitch/Alex - probability of transitions between clusters

all_transitions = []
for rat in x_array.id.unique():
    for condition in ["replete", "deplete"]:
        x_array_r = x_array.query("id == @rat & condition == @condition")
        transition_matrix = np.zeros((7,7))
        for i in range(len(x_array_r)-1):
            current_cluster = x_array_r.iloc[i].cluster
            next_cluster = x_array_r.iloc[i+1].cluster
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
f = plt.figure(figsize=(5,4))
gs = f.add_gridspec(2, 2, width_ratios=[15, 1])  # Adjust width_ratios for colorbar width

# f, [ax, cbar_ax] = plt.subplots(ncols=2, figsize=(5,4),
#                                 gridspec_kw={'left': 0.1, 'right': 0.9, 'top': 0.9, 'bottom': 0.1,
#                                              'width_ratios': [10, 1]})

# Add heatmap and colorbar axes
ax = f.add_subplot(gs[:, 0])
cbar_ax = f.add_subplot(gs[0, 1])

sns.heatmap(all_transitions_mean,
            cmap="Greys",
            vmin=0, vmax=0.7,
            annot=True,
            ax=ax, cbar_ax=cbar_ax,
            cbar_kws={"orientation": "vertical"},
            # linewidths=0.5, linecolor="black", annot_kws={"size": 8}, fmt=".2f"
)
ax.set_ylabel("From cluster")
ax.set_xlabel("To cluster")

ax.set_yticks(np.arange(7)+0.5, labels=np.arange(1,8), rotation=0)
ax.set_xticks(np.arange(7)+0.5, labels=np.arange(1,8))

cbar_ax.set_yticks([0, 0.7])
cbar_ax.set_ylabel("Transition probability", rotation=-90, labelpad=10, fontsize=10)
cbar_ax.tick_params(axis='both', which='both', length=0)

if savefigs:
        f.savefig(FIGSFOLDER / "cluster_transitions.pdf", dpi=600, transparent=True)

# %%
n = all_transitions_mean.shape[0]
sum_above_diagonal = np.sum(all_transitions_mean[np.triu_indices(n, k=1)])
sum_below_diagonal = np.sum(all_transitions_mean[np.tril_indices(n, k=-1)])


# %%
sum_above_diagonal
sum_below_diagonal

# %%
# Question 2 from Mitch/Alex - percent of rats in each cluster

prop_rats_per_cluster = []
for cluster in np.arange(0,7):
    print("Cluster", cluster)
    tmp = x_array.groupby(["cluster", "id", "condition"]).count().reset_index().query("cluster == @cluster")
    
    rep = len(tmp.query("condition == 'replete'"))
    dep = len(tmp.query("condition == 'deplete'"))
    prop_rats_per_cluster.append([rep, dep])

# %%
f, ax = plt.subplots(figsize=(5, 3),
                     gridspec_kw={'left': 0.15, 'right': 0.9, 'top': 0.85, 'bottom': 0.2}
                     )

spacing = 0.2
width = 0.35
for i in np.arange(0,7):

    ax.bar(i-spacing, prop_rats_per_cluster[i][0]/20, width=width, color=colors[0])
    ax.bar(i+spacing, prop_rats_per_cluster[i][1]/20, width=width, color=colors[2])

for i in np.arange(0,7):
    
    prop1 = prop_rats_per_cluster[i][0] / 20
    prop2 = prop_rats_per_cluster[i][1] / 20
    
    if prop1 == 0:
        ax.text(i-spacing, 0.02, f"{int(prop1*20)}", ha="center", va="bottom", fontsize=8, color=colors[0])
    else:
        ax.text(i-spacing, 0.02, f"{int(prop1*20)}", ha="center", va="bottom", fontsize=8, color="white")
        
    if prop2 == 0:
        ax.text(i+spacing, 0.02, f"{int(prop2*20)}", ha="center", va="bottom", fontsize=8, color=colors[2])
    else:
        ax.text(i+spacing, 0.02, f"{int(prop2*20)}", ha="center", va="bottom", fontsize=8, color="white")

    
ax.set_xticks(np.arange(0,7), labels=np.arange(1,8))
ax.set_ylabel("Proportion of rats")
ax.set_xlabel("Cluster")

ax.set_yticks([0,0.5,1])


sns.despine(ax=ax, top=True, right=True)

if savefigs:
    f.savefig(FIGSFOLDER / "nrats_per_cluster.pdf", dpi=600, transparent=True)

# %%
# code to work out how well the 2 clusters map on to the 7 clusters
x_array[x_array2.cluster == 0].groupby("cluster").count()

# %%
x_array[x_array2.cluster == 1].groupby("cluster").count()
