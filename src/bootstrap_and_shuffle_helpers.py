import random
import pandas as pd
import numpy as np
from scipy import stats

from model_fit_helpers import fit_sigmoid


def get_bootstrapped_distribution(df, trial_col, value_col, n_bootstraps=1000, random_seed=42):

    random.seed(random_seed)
    rats = df.id.unique()
    
    random_lists = [random.choices(list(rats), k=len(rats)) for _ in range(n_bootstraps)]
    bootstrapped_ks = []
    
    for i, random_list in enumerate(random_lists, 1):
        subset_bootstrapped = []
        for rat in random_list:
            rat_trials = df[df.id == rat]
            subset_bootstrapped.append(rat_trials)
            
        boot_df = pd.concat(subset_bootstrapped)
        boot_df = filter_complete_trials(boot_df, trial_col, nrats=len(rats))
        popt_boot = fit_sigmoid(boot_df, trial_col, value_col)
        k_boot = popt_boot[2] if np.all(np.isfinite(popt_boot)) else np.nan
        bootstrapped_ks.append(k_boot)
        
    print(f"Number of NaNs in bootstrapped distribution: {np.isnan(bootstrapped_ks).sum()}")
    
    return np.array(bootstrapped_ks)[np.isfinite(bootstrapped_ks)]  # also removes inf

def get_bootstrapped_distribution_using_slopes(df, trial_col, value_col, n_bootstraps=1000, random_seed=42):

    random.seed(random_seed)
    rats = df.id.unique()
    
    random_lists = [random.choices(list(rats), k=len(rats)) for _ in range(n_bootstraps)]
    bootstrapped_slopes = []
    
    for i, random_list in enumerate(random_lists, 1):
        subset_bootstrapped = []
        for rat in random_list:
            rat_trials = df[df.id == rat]
            subset_bootstrapped.append(rat_trials)
            
        boot_df = pd.concat(subset_bootstrapped)
        
        xvals = boot_df[trial_col]
        yvals = boot_df[value_col]
        slope, intercept, r, p, std_err = stats.linregress(xvals, yvals)
        
        slope = slope if np.all(np.isfinite(slope)) else np.nan
        
        
        bootstrapped_slopes.append(slope)

    print(f"Number of NaNs in bootstrapped distribution: {np.isnan(bootstrapped_slopes).sum()}")
    
    return np.array(bootstrapped_slopes)[np.isfinite(bootstrapped_slopes)]  # also removes inf

def filter_complete_trials(df, trial_col, nrats=None):
    trial_counts = df.groupby([trial_col]).size().reset_index(name='count')
    if nrats is None:
        nrats = len(df.id.unique())
    complete_trials = trial_counts[trial_counts['count'] == nrats]
    return df.merge(complete_trials[[trial_col]], on=trial_col, how='inner')