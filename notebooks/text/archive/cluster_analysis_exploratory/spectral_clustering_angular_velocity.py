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
import os
os.environ["OMP_NUM_THREADS"] = "2"  # To avoide warning in clustering, can try removing

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import trompy as tp

import dill

# %run "..//src//plots_for_spectral"

DATAFOLDER = Path("..//data")
RESULTSFOLDER = Path("..//results")

# %%
# for clustering all trials with no averaging
with open(DATAFOLDER / "snips_angvel_selected_conditions.pickle", "rb") as f:
    data = dill.load(f)

snips_10NaCl = data["snips_10NaCl_selected"]
snips_45NaCl = data["snips_45NaCl_selected"]
x_10NaCl = data["x_10NaCl_selected"]
x_45NaCl = data["x_45NaCl_selected"]

# %%
with open(DATAFOLDER / "snips_angvel.pickle", "rb") as f:
    snips_vel = dill.load(f)
    
with open(DATAFOLDER / "x_array_forDLC.pickle", "rb") as f:
    x_vel = dill.load(f)

snips_vel.shape

## pre-processing of the DLC angvel snips

# interpolation
df_snips_vel = pd.DataFrame(snips_vel)


df_interpolated = df_snips_vel.interpolate(method='linear', axis=1)
df_filled = df_interpolated.ffill(axis=1).bfill(axis=1)

snips_vel_processed = df_filled.to_numpy()

# find absolute
snips_vel_processed = np.abs(snips_vel_processed)

# how about instead of zscoring, we just adjust to baseline subtraction
baseline = np.nanmean(snips_vel_processed[:, :50], axis=1)
snips_vel_processed = snips_vel_processed - baseline[:, None]

print(np.sum(np.isnan(snips_vel_processed)))
snips_vel = np.array(snips_vel_processed)

# removing NaNs
rows_with_nans_mask = np.isnan(snips_vel).any(axis=1)
snips_vel = snips_vel[~rows_with_nans_mask]
x_vel = x_vel[~rows_with_nans_mask].reset_index(drop=True)
print(snips_vel.shape)
x_vel.shape


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
framerate = 10
pre_window_size = 5 * framerate # frames before infusion
window_size = 200 # Total number of frames plotted around an infusion
frames_for_infusion = 10 * framerate # 10 seconds until end of infusion

populationdata = snips_vel
sortwindow = [pre_window_size, pre_window_size + frames_for_infusion] # Sort responses between CS+ onset and reward.
sortresponse = np.argsort(np.mean(populationdata[:,sortwindow[0]:sortwindow[1]], axis=1))[::-1]
# sortresponse corresponds to an ordering of the neurons based on their average response in the sortwindow

cmax = 0.1 # Maximum colormap value. 

trial_types = ['trials']

# %%

fig, axs = plt.subplots(2,2,figsize=(3*2,3*2), sharex='all', sharey='row')
cbar_ax = fig.add_axes([0.86, .3, .01, .4])
cbar_ax.tick_params(width=0.5) 

for t in range(len(trial_types)):
    axs[0,t].set_title(trial_types[t])
    ax = axs[0,t]
    sns.heatmap(populationdata[sortresponse, t*window_size: (t+1)*window_size],
                ax=ax,
                cmap=plt.get_cmap('coolwarm'),
                vmin=-cmax,
                vmax=cmax,
                cbar=(t==0),
                cbar_ax=cbar_ax if (t==0) else None,
                cbar_kws={'label': 'Normalized fluorescence'})
    ax.grid(False)
    ax.tick_params(width=0.5)   
    ax.set_xticks([0, pre_window_size, pre_window_size + frames_for_infusion, window_size]) 
    ax.set_xticklabels([str(int((a-pre_window_size+0.0)/framerate))
                                     for a in [0, pre_window_size,
                                               pre_window_size + frames_for_infusion, window_size]])
    ax.set_yticks([])
    ax.axvline(pre_window_size, linestyle='--', color='k', linewidth=0.5)    
    ax.axvline(pre_window_size + frames_for_infusion, linestyle='--', color='k', linewidth=0.5)    
    ax.set_xlabel('Time from inf (s)')
    
        
    ax = axs[1,t]
    tsplot(ax, populationdata[sortresponse, t*window_size:(t+1)*window_size],
              )
    ax.axvline(pre_window_size, linestyle='--', color='k', linewidth=0.5)    
    ax.axvline(pre_window_size + frames_for_infusion, linestyle='--', color='k', linewidth=0.5)    
    ax.set_xlabel('Time from inf (s)')
    ax.set_xticks([0, pre_window_size, pre_window_size + frames_for_infusion, window_size]) 
    ax.set_xticklabels([str(int((a-pre_window_size+0.0)/framerate))
                                     for a in [0, pre_window_size,
                                               pre_window_size + frames_for_infusion, window_size]])
    
axs[0,0].set_ylabel('Trials')
axs[1,0].set_ylabel('Mean norm. fluor.')

fig.tight_layout()
fig.subplots_adjust(right=0.82)

# %%
load_savedpca_or_dopca = 'dopca'
# Select 'dopca' for doing PCA on the data. Select 'savedpca' for loading my previous results

if load_savedpca_or_dopca == 'dopca':
    pca = PCA(n_components=populationdata.shape[1], whiten=True)
    pca.fit(populationdata) 
    with open(RESULTSFOLDER / 'pcaresults.pickle', 'wb') as f:
        pickle.dump(pca, f)
elif load_savedpca_or_dopca == 'savedpca':
    with open(RESULTSFOLDER / 'pcaresults.pickle', 'rb') as f:
        pca = pickle.load(f)
    
transformed_data = pca.transform(populationdata)
with open(RESULTSFOLDER / 'transformed_data.pickle', 'wb') as f:
        pickle.dump(transformed_data, f)

pca_vectors = pca.components_
print('Number of PCs = %d'%(pca_vectors.shape[0]))

x = 100*pca.explained_variance_ratio_
xprime = x - (x[0] + (x[-1]-x[0])/(x.size-1)*np.arange(x.size))
num_retained_pcs = np.argmin(xprime)
# Number of PCs to be kept is defined as the number at which the 
# scree plot bends. This is done by simply bending the scree plot
# around the line joining (1, variance explained by first PC) and
# (num of PCs, variance explained by the last PC) and finding the 
# number of components just below the minimum of this rotated plot
print('Number of PCs to keep = %d'%(num_retained_pcs))

fig, ax = plt.subplots(figsize=(2,2))
ax.plot(np.arange(pca.explained_variance_ratio_.shape[0]).astype(int)+1, x, 'k')
ax.set_ylabel('Percentage of\nvariance explained')
ax.set_xlabel('PC number')
ax.axvline(num_retained_pcs, linestyle='--', color='k', linewidth=0.5)
ax.set_title('Scree plot')
# ax.set_xlim([0,50])
#[i.set_linewidth(0.5) for i in ax.spines.itervalues()]
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

fig.subplots_adjust(left=0.3)
fig.subplots_adjust(right=0.98)
fig.subplots_adjust(bottom=0.25)
fig.subplots_adjust(top=0.9)

colors_for_key = {}
colors_for_key['trials'] = (0,0.5,1)

numcols = 3.0
fig, axs = plt.subplots(int(np.ceil(num_retained_pcs/numcols)), int(numcols), sharey='all',
                        figsize=(2*numcols, 2*int(np.ceil(num_retained_pcs/numcols))))
for pc in range(num_retained_pcs):
    ax = axs.flat[pc]
    for k, tempkey in enumerate(trial_types):
        print(tempkey)
        ax.plot(pca_vectors[pc, k*window_size:(k+1)*window_size], color=colors_for_key[tempkey],
                label='PC %d: %s'%(pc+1, tempkey))
    ax.axvline(pre_window_size, linestyle='--', color='k', linewidth=1)
    # ax.annotate(s='PC %d'%(pc+1), xy=(0.45, 0.06), xytext=(0.45, 0.06), xycoords='axes fraction',
    #            textcoords='axes fraction', multialignment='center', size='large')
    if pc >= num_retained_pcs-numcols:
        ax.set_xticks([0, pre_window_size,
                       pre_window_size + frames_for_infusion, window_size])
        ax.set_xticklabels([str(int((a-pre_window_size+0.0)/framerate))
                             for a in [0, pre_window_size,
                                       pre_window_size + frames_for_infusion, window_size]])
    else:
        ax.set_xticks([])
        ax.xaxis.set_ticks_position('none')
    if pc%numcols:
        ax.yaxis.set_ticks_position('none')
    # [i.set_linewidth(0.5) for i in ax.spines.itervalues()]
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)


fig.text(0.5, 0.05, 'Time from inf (s)', horizontalalignment='center', rotation='horizontal')
fig.text(0.02, 0.6, 'PCA weights', verticalalignment='center', rotation='vertical')
fig.tight_layout()
for ax in axs.flat[num_retained_pcs:]:
    ax.set_visible(False)

fig.subplots_adjust(wspace=0.08, hspace=0.08)
fig.subplots_adjust(bottom=0.13)

# %%
transformed_data = transformed_data[:, :num_retained_pcs]

# %%
np.min(transformed_data)

# %%
### Best silhouette score was with spectral clustering using a sigmoid and discretize
# when all trials are included
# 
max_n_clusters = 9 # Maximum number of clusters expected. I already ran this with up to 20 clusters and know
# that the optimal number is 9. So, I am leaving this at 11. In your data, might be worth increasing this, but
# it will take more time to run.

possible_n_clusters = np.arange(2, max_n_clusters+1) #This requires a minimum of 2 clusters.
# When the data contain no clusters at all, it will be quite visible when inspecting the two obtained clusters, 
# as the responses of the clusters will be quite similar. This will also be visible when plotting the data in
# the reduced dimensionality PC space (done below).

silhouette_scores = np.nan*np.ones((possible_n_clusters.size,
                                    ))

possible_n_nearest_neighbors = np.array([400])

for n_clustersidx, n_clusters in enumerate(possible_n_clusters):
        model = SpectralClustering(n_clusters=n_clusters,
                                   affinity='sigmoid',
                                   assign_labels='discretize',
                                   random_state=0)
        model.fit(transformed_data[:,:num_retained_pcs])
        silhouette_scores[n_clustersidx] = silhouette_score(transformed_data[:,:num_retained_pcs],
                                                                   model.labels_,
                                                                   metric='cosine')
        print('Done with numclusters = %d: score = %.3f'%(n_clusters, silhouette_scores[n_clustersidx,                                                                           
                                                                                                ]))

print('Done with model fitting')

# %%
temp = {}
temp['possible_n_clusters'] = possible_n_clusters
temp['silhouette_scores'] = silhouette_scores
temp['shape'] = 'cluster_nn'
with open(RESULTSFOLDER / 'silhouette_scores.pickle', 'wb') as f:
    pickle.dump(temp, f)

# %%
with open(RESULTSFOLDER / 'silhouette_scores.pickle', 'rb') as f:
    silhouette_scores = pickle.load(f)
    
# Identify optimal parameters from the above parameter space
temp = np.where(silhouette_scores['silhouette_scores']==np.nanmax(silhouette_scores['silhouette_scores']))
n_clusters = silhouette_scores['possible_n_clusters'][temp[0][0]]

print (n_clusters)
    
# transformed_data = np.load(os.path.join(ROOT_FOLDER, 'transformed_data.npy'))

# Redo clustering with these optimal parameters
model = SpectralClustering(n_clusters=2,
                           affinity='sigmoid',
                           assign_labels='discretize',)

model.fit(transformed_data[:,:num_retained_pcs])

temp = silhouette_score(transformed_data[:,:num_retained_pcs], model.labels_, metric='cosine')

print('Number of clusters = %d, average silhouette = %.3f'%(len(set(model.labels_)), temp))

# Save this optimal clustering model.
# with open(os.path.join(basedir, 'clusteringmodel.pickle'), 'wb') as f:
#     pickle.dump(model, f)

          
# Since the clustering labels are arbitrary, I rename the clusters so that the first cluster will have the most
# positive response and the last cluster will have the most negative response.
def reorder_clusters(rawlabels):
    uniquelabels = list(set(rawlabels))
    responses = np.nan*np.ones((len(uniquelabels),))
    for l, label in enumerate(uniquelabels):
        responses[l] = np.mean(populationdata[rawlabels==label, pre_window_size:2*pre_window_size])
    temp = np.argsort(responses).astype(int)[::-1]
    temp = np.array([np.where(temp==a)[0][0] for a in uniquelabels])
    outputlabels = np.array([temp[a] for a in list(np.digitize(rawlabels, uniquelabels)-1)])
    return outputlabels
newlabels = reorder_clusters(model.labels_)

# Create a new variable containing all unique cluster labels
uniquelabels = list(set(newlabels))

np.save(RESULTSFOLDER / 'clusterlabels.npy', newlabels)

colors_for_cluster = [[0.933, 0.250, 0.211],
                      "green",
                      [0.941, 0.352, 0.156],
                      [0.964, 0.572, 0.117],
                      [0.980, 0.686, 0.250],
                      [0.545, 0.772, 0.247],
                      [0.215, 0.701, 0.290],
                      [0, 0.576, 0.270],
                      [0, 0.650, 0.611],
                      [0.145, 0.662, 0.878]]


# %%
def make_silhouette_plot(X, cluster_labels):
    
    n_clusters = len(set(cluster_labels))
    
    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(4, 4)

    # The 1st subplot is the silhouette plot
    # The silhouette coefficient can range from -1, 1 but in this example all
    # lie within [-0.1, 1]
    ax.set_xlim([-0.4, 1])
    # The (n_clusters+1)*10 is for inserting blank space between silhouette
    # plots of individual clusters, to demarcate them clearly.
    ax.set_ylim([0, len(X) + (n_clusters + 1) * 10])
    silhouette_avg = silhouette_score(X, cluster_labels, metric='cosine')

    # Compute the silhouette scores for each sample
    sample_silhouette_values = silhouette_samples(X, cluster_labels, metric='cosine')

    y_lower = 10
    for i in range(n_clusters):
        # Aggregate the silhouette scores for samples belonging to
        # cluster i, and sort them
        ith_cluster_silhouette_values = \
            sample_silhouette_values[cluster_labels == i]

        ith_cluster_silhouette_values.sort()

        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i

        color = colors_for_cluster[i]
        ax.fill_betweenx(np.arange(y_lower, y_upper),
                          0, ith_cluster_silhouette_values,
                          facecolor=color, edgecolor=color, alpha=0.9)

        # Label the silhouette plots with their cluster numbers at the middle
        ax.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i+1))

        # Compute the new y_lower for next plot
        y_lower = y_upper + 10  # 10 for the 0 samples

    ax.set_title("The silhouette plot for the various clusters.")
    ax.set_xlabel("The silhouette coefficient values")
    ax.set_ylabel("Cluster label")

    # The vertical line for average silhouette score of all the values
    ax.axvline(x=silhouette_avg, color="red", linestyle="--")

    ax.set_yticks([])  # Clear the yaxis labels / ticks
    ax.set_xticks([-0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1])


# %%
make_silhouette_plot(transformed_data[:,:num_retained_pcs], model.labels_)

# %%
cmax = 0.1
sortwindow = [50, 100]

fig, axs = plt.subplots(len(trial_types),len(uniquelabels),
                        figsize=(2*len(uniquelabels),2*len(trial_types)),
                        gridspec_kw={"top": 0.8})

cbar_ax = fig.add_axes([.94, .3, .01, .4])
cbar_ax.tick_params(width=0.5) 

numroisincluster = np.nan*np.ones((len(uniquelabels),))
for c, cluster in enumerate(uniquelabels):
    k=0
    tempkey='cas'
    temp = populationdata[np.where(newlabels==cluster)[0], k*window_size:(k+1)*window_size]
    numroisincluster[c] = temp.shape[0]
    sortresponse = np.argsort(np.mean(temp[:,sortwindow[0]:sortwindow[1]], axis=1))[::-1]
    sns.heatmap(temp[sortresponse],
                ax=axs[cluster],
                cmap=plt.get_cmap('coolwarm'),
                vmin=-cmax,
                vmax=cmax,
                cbar=(cluster==0),
                cbar_ax=cbar_ax if (cluster==0) else None,
                cbar_kws={'label': 'Normalized fluorescence'})
    axs[cluster].grid(False)
    
    if k==len(trial_types)-1:
        axs[cluster].set_xticks([0, pre_window_size,
                                    pre_window_size + frames_for_infusion, window_size])
    else:
        axs[cluster].set_xticks([])
    axs[cluster].tick_params(width=0.5)    
    axs[cluster].set_xticklabels([str(int((a-pre_window_size+0.0)/framerate))
                                     for a in [0, pre_window_size,
                                               pre_window_size + frames_for_infusion, window_size]])
    axs[cluster].set_yticks([])
    axs[cluster].axvline(pre_window_size, linestyle='--', color='k', linewidth=0.5)
    axs[cluster].axvline(pre_window_size + frames_for_infusion, linestyle='--', color='k', linewidth=0.5)
    if cluster==0:
        axs[0].set_ylabel('%s\nNeurons'%(tempkey))
    axs[cluster].set_title('Cluster %d\n(n=%d)'%(cluster+1, numroisincluster[c]))

fig.text(0.5, 0.05, 'Time from lick (s)', fontsize=12,
         horizontalalignment='center', verticalalignment='center', rotation='horizontal')
fig.tight_layout()

fig.subplots_adjust(wspace=0.1, hspace=0.1)
fig.subplots_adjust(left=0.03)
fig.subplots_adjust(right=0.93)
fig.subplots_adjust(bottom=0.2)
fig.subplots_adjust(top=0.83)

fig.savefig(RESULTSFOLDER / "heatplots_clusters.png")

# %%
cmax = 0.1
sortwindow = [15, 100]

fig, axs = plt.subplots(len(trial_types),len(uniquelabels),
                        figsize=(2*len(uniquelabels),2*len(trial_types)),
                        gridspec_kw={"top": 0.8})

numroisincluster = np.nan*np.ones((len(uniquelabels),))
for c, cluster in enumerate(uniquelabels):
    k=0
    tempkey='cas'
    temp = populationdata[np.where(newlabels==cluster)[0], k*window_size:(k+1)*window_size]
    numroisincluster[c] = temp.shape[0]
    sortresponse = np.argsort(np.mean(temp[:,sortwindow[0]:sortwindow[1]], axis=1))[::-1]
    
    tp.shadedError(axs[cluster], temp, linecolor=colors_for_cluster[c])
    print(temp.shape)
    
    axs[cluster].grid(False)
    
    if k==len(trial_types)-1:
        axs[cluster].set_xticks([0, pre_window_size,
                                    pre_window_size + frames_for_infusion, window_size])
    else:
        axs[cluster].set_xticks([])
    axs[cluster].tick_params(width=0.5)    
    axs[cluster].set_xticklabels([str(int((a-pre_window_size+0.0)/framerate))
                                     for a in [0, pre_window_size,
                                               pre_window_size + frames_for_infusion, window_size]])
    
    axs[cluster].set_yticks([])
    axs[cluster].axvline(pre_window_size, linestyle='--', color='k', linewidth=0.5)
    axs[cluster].axvline(pre_window_size + frames_for_infusion, linestyle='--', color='k', linewidth=0.5)
    if cluster==0:
        axs[0].set_ylabel('%s\nNormalized fluorescence'%(tempkey))
    axs[cluster].set_title('Cluster %d\n(n=%d)'%(cluster+1, numroisincluster[c]))
    
    fig.text(0.5, 0.05, 'Time from lick (s)', fontsize=12,
         horizontalalignment='center', verticalalignment='center', rotation='horizontal')
    fig.tight_layout()

    fig.subplots_adjust(wspace=0.1, hspace=0.1)
    fig.subplots_adjust(left=0.03)
    fig.subplots_adjust(right=0.93)
    fig.subplots_adjust(bottom=0.2)
    fig.subplots_adjust(top=0.83)
    
fig.savefig(RESULTSFOLDER / "shadederror_clusters.png")

# %%
num_clusterpairs = len(uniquelabels)*(len(uniquelabels)-1)/2

numrows = int(np.ceil(num_clusterpairs**0.5))
numcols = int(np.ceil(num_clusterpairs/np.ceil(num_clusterpairs**0.5)))
fig, axs = plt.subplots(numrows, numcols, figsize=(3*numrows, 3*numcols))

tempsum = 0
for c1, cluster1 in enumerate(uniquelabels):
    for c2, cluster2 in enumerate(uniquelabels):
        if cluster1>=cluster2:
            continue
        temp1 = transformed_data[np.where(newlabels==cluster1)[0], :num_retained_pcs]
        temp2 = transformed_data[np.where(newlabels==cluster2)[0], :num_retained_pcs]
        X = np.concatenate((temp1, temp2), axis=0)
        tsne = TSNE(n_components=2, init='random',
                    random_state=0, perplexity=100)
        Y = tsne.fit_transform(X)
        #ax = axs[int(tempsum/numcols),
        #         int(tempsum - tempsum/numcols*numcols)]
        
        try:
            ax = axs.flat[tempsum]
        except:
            ax = axs
        ax.scatter(Y[:np.sum(newlabels==cluster1),0],
                   Y[:np.sum(newlabels==cluster1),1],
                   color=colors_for_cluster[cluster1], label='Cluster %d'%(cluster1+1), alpha=0.2)
        ax.scatter(Y[np.sum(newlabels==cluster1):,0],
                   Y[np.sum(newlabels==cluster1):,1],
                   color=colors_for_cluster[cluster2], label='Cluster %d'%(cluster2+1), alpha=0.2)

        ax.set_xlabel('tsne dimension 1')
        ax.set_ylabel('tsne dimension 2')
        ax.legend()
        tempsum += 1
fig.tight_layout()

# %%
x_10NaCl = x_10NaCl.assign(infusiontype="10NaCl")
x_45NaCl = x_45NaCl.assign(infusiontype="45NaCl")

x_array = (pd.concat([x_10NaCl, x_45NaCl], axis=0, ignore_index=True)
           .assign(cluster=newlabels)
)

with open(DATAFOLDER / "x_array_2clusters_alltrials.pickle", "wb") as f:
    dill.dump(x_array, f)

# %%
x_array
