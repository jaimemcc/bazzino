
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import subprocess
import time
import pandas
import pickle
import math
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC, SVR, LinearSVC
from sklearn.metrics import accuracy_score, silhouette_score, adjusted_rand_score, silhouette_samples
from sklearn.cluster import AgglomerativeClustering, SpectralClustering, KMeans
from sklearn.model_selection import KFold, LeaveOneOut, train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.kernel_ridge import KernelRidge
from sklearn import linear_model
from sklearn.manifold import TSNE
import scipy.stats as stats
import statsmodels.api as sm
import statsmodels.formula.api as smf
from patsy import (ModelDesc, EvalEnvironment, Term, EvalFactor, LookupFactor, dmatrices, INTERCEPT)
from statsmodels.distributions.empirical_distribution import ECDF
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.colorbar as colorbar
import sys
import re

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.colorbar as colorbar

import seaborn as sns
sns.set_style('ticks')

mpl.rcParams['axes.titlesize'] = 16
mpl.rcParams['axes.labelsize'] = 12
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12
mpl.rcParams['legend.fontsize'] = 12
mpl.rcParams['legend.labelspacing'] = 0.2
mpl.rcParams['axes.labelpad'] = 2
mpl.rcParams['xtick.major.size'] = 2
mpl.rcParams['xtick.major.width'] = 0.5
mpl.rcParams['xtick.major.pad'] = 1
mpl.rcParams['ytick.major.size'] = 2
mpl.rcParams['ytick.major.width'] = 0.5
mpl.rcParams['ytick.major.pad'] = 1
mpl.rcParams['lines.scale_dashes'] = False
mpl.rcParams['lines.dashed_pattern'] = (2, 1)
mpl.rcParams['font.sans-serif'] = ['Arial']
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['text.color'] = 'k'

def standardize_plot_graphics(ax):
    [i.set_linewidth(0.5) for i in ax.spines.itervalues()]
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    return ax

def fit_regression(x, y):
    lm = sm.OLS(y, sm.add_constant(x)).fit()
    x_range = sm.add_constant(np.array([x.min(), x.max()]))
    x_range_pred = lm.predict(x_range)
    return lm.pvalues[1], lm.params[1], x_range[:,1], x_range_pred, lm.rsquared

def CDFplot(x, ax, **kwargs):
    x = np.array(x)
    ix=np.argsort(x)
    ax.plot(x[ix], ECDF(x)(x)[ix], **kwargs)
    return ax

def fit_regression_and_plot(x, y, ax, plot_label='', color='k', linecolor='r', markersize=3,
                            show_pval=True):
    #linetype is a string like 'bo'
    pvalue, slope, temp, temppred, R2 = fit_regression(x, y)   
    if show_pval:
        plot_label = '%s p=%.2e\nr=%.3f'% (plot_label, pvalue, np.sign(slope)*np.sqrt(R2))
    else:
        plot_label = '%s r=%.3f'% (plot_label, np.sign(slope)*np.sqrt(R2))
    ax.scatter(x, y, color=color, label=plot_label, s=markersize)
    ax.plot(temp, temppred, color=linecolor)
    return ax, slope, pvalue, R2

def tsplot(ax, data,**kw):
    x = np.arange(data.shape[1])
    est = np.mean(data, axis=0)
    sd = np.std(data, axis=0)
    cis = (est - sd, est + sd)
    ax.fill_between(x,cis[0],cis[1],alpha=0.2, **kw)
    ax.plot(x,est,**kw)
    ax.margins(x=0)

def make_tsplot(populationdata,
           sortresponse,
           framerate,
           pre_window_size,
           window_size,
           frames_for_infusion,
           cmax,
           trial_types,
           **kwargs):

    f, axs = plt.subplots(nrows=2, figsize=(3,6), sharex='all', sharey='row')
    cbar_ax = f.add_axes([0.86, .3, .01, .4])
    cbar_ax.tick_params(width=0.5)

    for t in range(len(trial_types)):
        axs[0].set_title(trial_types[t])
        ax = axs[0]
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
        
            
        ax = axs[1]
        tsplot(ax, populationdata[sortresponse, t*window_size:(t+1)*window_size],
                )
        ax.axvline(pre_window_size, linestyle='--', color='k', linewidth=0.5)    
        ax.axvline(pre_window_size + frames_for_infusion, linestyle='--', color='k', linewidth=0.5)    
        ax.set_xlabel('Time from inf (s)')
        ax.set_xticks([0, pre_window_size, pre_window_size + frames_for_infusion, window_size]) 
        ax.set_xticklabels([str(int((a-pre_window_size+0.0)/framerate))
                                        for a in [0, pre_window_size,
                                                pre_window_size + frames_for_infusion, window_size]])
        
    axs[0].set_ylabel('Trials')
    axs[1].set_ylabel('Mean norm. fluor.')

    f.tight_layout()
    f.subplots_adjust(right=0.82)

    return f, axs

def plot_scree(x,
               pca,
               num_retained_pcs):
        
    f, ax = plt.subplots(figsize=(2,2))
    ax.plot(np.arange(pca.explained_variance_ratio_.shape[0]).astype(int)+1, x, 'k')
    ax.set_ylabel('Percentage of\nvariance explained')
    ax.set_xlabel('PC number')
    ax.axvline(num_retained_pcs, linestyle='--', color='k', linewidth=0.5)
    ax.set_title('Scree plot')
    # ax.set_xlim([0,50])
    #[i.set_linewidth(0.5) for i in ax.spines.itervalues()]
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    f.subplots_adjust(left=0.3)
    f.subplots_adjust(right=0.98)
    f.subplots_adjust(bottom=0.25)
    f.subplots_adjust(top=0.9)

    return f, ax

def plot_pca_vectors(pca_vectors,
                     num_retained_pcs,
                     framerate,
                     pre_window_size,
                     window_size,
                     frames_for_infusion,):

    numcols = 3.0
    f, axs = plt.subplots(int(np.ceil(num_retained_pcs/numcols)), int(numcols), sharey='all',
                            figsize=(2*numcols, 2*int(np.ceil(num_retained_pcs/numcols))))
    for pc in range(num_retained_pcs):
        ax = axs.flat[pc]

        ax.plot(pca_vectors[pc, :], color="k",
                label='PC %d: %s'%(pc+1, "trial"))
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


    f.text(0.5, 0.05, 'Time from inf (s)', horizontalalignment='center', rotation='horizontal')
    f.text(0.02, 0.6, 'PCA weights', verticalalignment='center', rotation='vertical')
    f.tight_layout()
    for ax in axs.flat[num_retained_pcs:]:
        ax.set_visible(False)

    f.subplots_adjust(wspace=0.08, hspace=0.08)
    f.subplots_adjust(bottom=0.13)

    return f, axs

def make_silhouette_plot(X,
                         cluster_labels,
                         colors_for_cluster=["red", "green", "blue", "orange", "purple", "brown"]
                         ):
    
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

def show_clusters(populationdata,
                  newlabels,
                  uniquelabels,
                  pre_window_size,
                  frames_for_infusion,
                  window_size,
                  framerate,
                  cmax,
                  colors_for_cluster=["red", "green", "blue", "orange", "purple", "brown", "pink", "gray", "cyan"]
                  ):

    f, ax = plt.subplots(ncols=2, nrows=2, figsize=(5,5))
    
    # make heatmaps first and then line plots
    for c, cluster in enumerate(uniquelabels):
        k=0
        temp = populationdata[np.where(newlabels==cluster)[0], k*window_size:(k+1)*window_size]
        sortresponse = np.argsort(np.mean(temp[:,pre_window_size:pre_window_size + frames_for_infusion], axis=1))[::-1]
        
        axis = ax[0, c]
        axis.set_title('Cluster %d\n(n=%d)'%(cluster+1, temp.shape[0]))
        sns.heatmap(temp[sortresponse],
                    ax=axis,
                    cmap=plt.get_cmap('coolwarm'),
                    vmin=-cmax,
                    vmax=cmax,
                    cbar=(c==1),
                    cbar_kws={'label': 'Normalized fluorescence'})
        
        axis = ax[1, c]
        axis.plot(np.arange(0, window_size)/framerate - pre_window_size/framerate,
                      np.mean(temp, axis=0),
                      color=colors_for_cluster[c],)
        axis.fill_between(np.arange(0, window_size)/framerate - pre_window_size/framerate,
                          np.mean(temp, axis=0) - np.std(temp, axis=0),
                          np.mean(temp, axis=0) + np.std(temp, axis=0),
                          color=colors_for_cluster[c], alpha=0.2)
        
        axis.axvline(-pre_window_size/framerate, linestyle='--', color='k', linewidth=0.5)
        axis.axvline(-pre_window_size/framerate + frames_for_infusion/framerate, linestyle='--', color='k', linewidth=0.5)
        
        sns.despine(ax=axis, left=True, bottom=True)
    
    ax[1,1].sharey(ax[1,0])
     