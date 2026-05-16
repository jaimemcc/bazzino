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
# # Compare Angular Velocity vs Movement
#
# Comparison of head rotation (angular velocity) vs body movement metrics across individual rats and conditions. Explore different threshold values for both metrics.

# %%
import os
import sys
import pickle
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path.cwd().parent / 'src'))

# Configuration
DATAFOLDER = Path.cwd().parent / 'data'

# Plotting settings
sns.set_style('whitegrid')
plt.rcParams['figure.figsize'] = (14, 8)
plt.rcParams['font.size'] = 10

# %%
# Load assembled data
pkl_file = DATAFOLDER / 'assembled_data.pickle'
with open(pkl_file, 'rb') as f:
    data_dict = pickle.load(f)

# Extract key components
x_array = data_dict['x_array']
snips_movement = data_dict['snips_movement']
snips_angvel = data_dict['snips_angvel']
snips_photo = data_dict['snips_photo']

print(f"x_array shape: {x_array.shape}")
print(f"x_array columns: {list(x_array.columns)}")
print(f"snips_movement shape: {snips_movement.shape}")
print(f"snips_angvel shape: {snips_angvel.shape}")
print(f"\nFirst few rows of x_array:")
print(x_array.head())

# %%
# Filter to deplete/45NaCl condition
x_filtered = x_array.query("condition == 'deplete' & infusiontype == '45NaCl'")
print(f"Filtered to: {len(x_filtered)} trials")
print(f"Number of unique rats: {x_filtered.id.nunique()}")
print(f"Rats: {sorted(x_filtered.id.unique())}")
print(f"\nAUC movement - mean: {x_filtered['auc_movement'].mean():.3f}, std: {x_filtered['auc_movement'].std():.3f}")
print(f"AUC angular velocity - mean: {x_filtered['auc_angvel'].mean():.3f}, std: {x_filtered['auc_angvel'].std():.3f}")

# %%
# Mean snips comparison
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Movement snips
trial_indices = x_filtered.index.values
movement_snips_filtered = snips_movement[trial_indices]
mean_movement = movement_snips_filtered.mean(axis=0)
std_movement = movement_snips_filtered.std(axis=0)

axes[0].plot(mean_movement, label='Mean', linewidth=2, color='blue')
axes[0].fill_between(range(len(mean_movement)), 
                      mean_movement - std_movement, 
                      mean_movement + std_movement, 
                      alpha=0.3, color='blue', label='±1 std')
axes[0].set_title(f'Movement Snips', fontsize=12, fontweight='bold')
axes[0].set_xlabel('Time bins')
axes[0].set_ylabel('Movement (normalized)')
axes[0].legend()
axes[0].grid(True, alpha=0.3)

# Angular velocity snips
angvel_snips_filtered = snips_angvel[trial_indices]
mean_angvel = angvel_snips_filtered.mean(axis=0)
std_angvel = angvel_snips_filtered.std(axis=0)

axes[1].plot(mean_angvel, label='Mean', linewidth=2, color='red')
axes[1].fill_between(range(len(mean_angvel)), 
                      mean_angvel - std_angvel, 
                      mean_angvel + std_angvel, 
                      alpha=0.3, color='red', label='±1 std')
axes[1].set_title(f'Angular Velocity Snips ({condition})', fontsize=12, fontweight='bold')
axes[1].set_xlabel('Time bins')
axes[1].set_ylabel('Angular velocity (deg/frame)')
axes[1].legend()
axes[1].grid(True, alpha=0.3)

print(f"Movement snips: n_trials={len(movement_snips_filtered)}, n_bins={len(mean_movement)}")
print(f"Angular velocity snips: n_trials={len(angvel_snips_filtered)}, n_bins={len(mean_angvel)}")

# %%
# Individual rat heatmaps
rats = sorted(x_filtered.id.unique())
n_rats = len(rats)

fig, axes = plt.subplots(n_rats, 2, figsize=(12, 3*n_rats))
if n_rats == 1:
    axes = axes[np.newaxis, :]

for i, rat in enumerate(rats):
    rat_trials = x_filtered[x_filtered.id == rat].index.values
    
    # Movement heatmap
    movement_rat = snips_movement[rat_trials].values
    sns.heatmap(movement_rat, ax=axes[i, 0], cmap='YlOrRd', cbar_kws={'label': 'Movement'})
    axes[i, 0].set_title(f'Rat {rat} - Movement ({len(rat_trials)} trials)', fontweight='bold')
    axes[i, 0].set_xlabel('Time bins')
    axes[i, 0].set_ylabel('Trial')
    
    # Angular velocity heatmap
    angvel_rat = snips_angvel[rat_trials].values
    sns.heatmap(angvel_rat, ax=axes[i, 1], cmap='RdYlBu_r', cbar_kws={'label': 'Ang Vel (deg/frame)'})
    axes[i, 1].set_title(f'Rat {rat} - Angular Velocity ({len(rat_trials)} trials)', fontweight='bold')
    axes[i, 1].set_xlabel('Time bins')
    axes[i, 1].set_ylabel('Trial')

plt.tight_layout()
plt.savefig(RESULTS_FOLDER / 'individual_rat_heatmaps.png', dpi=150, bbox_inches='tight')
plt.show()


# %%
# Threshold exploration functions
def calculate_time_above_threshold(snips_df, threshold):
    """Calculate proportion of bins above threshold for each trial"""
    return (snips_df > threshold).sum(axis=1) / len(snips_df.columns)

def calculate_threshold_metrics(x_data, snips_movement_df, snips_angvel_df, 
                               movement_threshold, angvel_threshold):
    """Calculate metrics at specific thresholds"""
    trial_indices = x_data.index.values
    
    movement_thresh = calculate_time_above_threshold(
        snips_movement_df.loc[trial_indices], movement_threshold)
    angvel_thresh = calculate_time_above_threshold(
        snips_angvel_df.loc[trial_indices], angvel_threshold)
    
    results = x_data.copy()
    results['prop_above_movement_thresh'] = movement_thresh
    results['prop_above_angvel_thresh'] = angvel_thresh
    
    return results

# Test thresholds
mov_thresholds = np.array([0.01, 0.02, 0.05, 0.1])
angvel_thresholds = np.array([0.5, 1.0, 1.5, 2.0])

print("Movement thresholds:", mov_thresholds)
print("Angular velocity thresholds (deg/frame):", angvel_thresholds)

# %%
# Threshold exploration - effect on trial distribution
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Movement threshold exploration
for thresh in mov_thresholds:
    mov_metric = calculate_time_above_threshold(movement_snips_filtered, thresh)
    axes[0, 0].hist(mov_metric, bins=20, alpha=0.5, label=f'thresh={thresh}')
axes[0, 0].set_xlabel('Proportion of bins above threshold')
axes[0, 0].set_ylabel('Number of trials')
axes[0, 0].set_title('Movement Threshold Exploration')
axes[0, 0].legend()
axes[0, 0].grid(True, alpha=0.3)

# Angular velocity threshold exploration
for thresh in angvel_thresholds:
    angvel_metric = calculate_time_above_threshold(angvel_snips_filtered, thresh)
    axes[0, 1].hist(angvel_metric, bins=20, alpha=0.5, label=f'thresh={thresh}')
axes[0, 1].set_xlabel('Proportion of bins above threshold')
axes[0, 1].set_ylabel('Number of trials')
axes[0, 1].set_title('Angular Velocity Threshold Exploration')
axes[0, 1].legend()
axes[0, 1].grid(True, alpha=0.3)

# Mean values at different thresholds
mean_mov = [calculate_time_above_threshold(movement_snips_filtered, t).mean() for t in mov_thresholds]
axes[1, 0].plot(mov_thresholds, mean_mov, 'o-', linewidth=2, markersize=8, color='blue')
axes[1, 0].set_xlabel('Threshold')
axes[1, 0].set_ylabel('Mean proportion above threshold')
axes[1, 0].set_title('Movement: Mean Activity vs Threshold')
axes[1, 0].grid(True, alpha=0.3)

mean_angvel = [calculate_time_above_threshold(angvel_snips_filtered, t).mean() for t in angvel_thresholds]
axes[1, 1].plot(angvel_thresholds, mean_angvel, 'o-', linewidth=2, markersize=8, color='red')
axes[1, 1].set_xlabel('Threshold (deg/frame)')
axes[1, 1].set_ylabel('Mean proportion above threshold')
axes[1, 1].set_title('Angular Velocity: Mean Activity vs Threshold')
axes[1, 1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(RESULTS_FOLDER / 'threshold_exploration.png', dpi=150, bbox_inches='tight')
plt.show()

# %%
# Correlation analysis between movement and angular velocity
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# Scatter: AUC movement vs AUC angular velocity
axes[0].scatter(x_filtered['auc_movement'], x_filtered['auc_angvel'], 
               s=100, alpha=0.6, c=pd.Categorical(x_filtered['rat_id']).codes, cmap='tab10')
corr = x_filtered['auc_movement'].corr(x_filtered['auc_angvel'])
axes[0].set_xlabel('AUC Movement')
axes[0].set_ylabel('AUC Angular Velocity')
axes[0].set_title(f'AUC Correlation (r={corr:.3f})')
axes[0].grid(True, alpha=0.3)

# Per-rat comparison of mean metrics
rat_comparison = []
for rat in rats:
    rat_mask = x_filtered['rat_id'] == rat
    rat_comparison.append({
        'rat_id': rat,
        'mean_auc_movement': x_filtered[rat_mask]['auc_movement'].mean(),
        'mean_auc_angvel': x_filtered[rat_mask]['auc_angvel'].mean(),
        'n_trials': rat_mask.sum()
    })
rat_df = pd.DataFrame(rat_comparison)

x_pos = np.arange(len(rats))
width = 0.35
ax2 = axes[1]
ax2_twin = ax2.twinx()

bars1 = ax2.bar(x_pos - width/2, rat_df['mean_auc_movement'], width, 
                label='Movement', alpha=0.8, color='blue')
bars2 = ax2_twin.bar(x_pos + width/2, rat_df['mean_auc_angvel'], width, 
                     label='Angular Velocity', alpha=0.8, color='red')

ax2.set_xlabel('Rat ID')
ax2.set_ylabel('Mean AUC Movement', color='blue')
ax2_twin.set_ylabel('Mean AUC Angular Velocity', color='red')
ax2.set_title('Mean Metrics by Rat')
ax2.set_xticks(x_pos)
ax2.set_xticklabels(rats)
ax2.tick_params(axis='y', labelcolor='blue')
ax2_twin.tick_params(axis='y', labelcolor='red')
ax2.grid(True, alpha=0.3, axis='y')

# Ratio of activity metrics
rat_df['movement_to_angvel_ratio'] = rat_df['mean_auc_movement'] / rat_df['mean_auc_angvel']
axes[2].bar(rats, rat_df['movement_to_angvel_ratio'], alpha=0.8, color='purple')
axes[2].set_xlabel('Rat ID')
axes[2].set_ylabel('Movement / Angular Velocity Ratio')
axes[2].set_title('Relative Activity: Movement vs Rotation')
axes[2].grid(True, alpha=0.3, axis='y')
axes[2].axhline(y=rat_df['movement_to_angvel_ratio'].mean(), 
                color='k', linestyle='--', label=f'Mean={rat_df["movement_to_angvel_ratio"].mean():.2f}')
axes[2].legend()

plt.tight_layout()
plt.savefig(RESULTS_FOLDER / 'correlation_analysis.png', dpi=150, bbox_inches='tight')
plt.show()

print("\nPer-rat comparison:")
print(rat_df.to_string(index=False))

# %%
# Summary statistics table
summary_stats = []

for mov_thresh in mov_thresholds:
    for angvel_thresh in angvel_thresholds:
        results = calculate_threshold_metrics(
            x_filtered, movement_snips_filtered, angvel_snips_filtered,
            mov_thresh, angvel_thresh)
        
        summary_stats.append({
            'movement_threshold': mov_thresh,
            'angvel_threshold': angvel_thresh,
            'n_high_movement': (results['prop_above_movement_thresh'] > 0.5).sum(),
            'n_high_angvel': (results['prop_above_angvel_thresh'] > 0.5).sum(),
            'mean_movement_activity': results['prop_above_movement_thresh'].mean(),
            'mean_angvel_activity': results['prop_above_angvel_thresh'].mean(),
        })

summary_df = pd.DataFrame(summary_stats)

print("Summary: Trial counts with high activity (>50% bins above threshold)")
print("\nMovement threshold vs Angular velocity threshold:")
pivot_table = summary_df.pivot_table(
    index='movement_threshold', 
    columns='angvel_threshold',
    values='n_high_movement',
    aggfunc='first'
)
print("\nNumber of trials with high MOVEMENT (>50% bins above threshold):")
print(pivot_table)

pivot_angvel = summary_df.pivot_table(
    index='movement_threshold', 
    columns='angvel_threshold',
    values='n_high_angvel',
    aggfunc='first'
)
print("\nNumber of trials with high ANGULAR VELOCITY (>50% bins above threshold):")
print(pivot_angvel)

# %%
# Interactive threshold testing
# Try different thresholds here
test_movement_threshold = 0.02
test_angvel_threshold = 1.0

print(f"\n--- Testing Movement Threshold: {test_movement_threshold} ---")
movement_above = calculate_time_above_threshold(movement_snips_filtered, test_movement_threshold)
print(f"Trials above threshold: {(movement_above > 0.5).sum()} / {len(movement_above)}")
print(f"Mean proportion active: {movement_above.mean():.3f}")
print(f"Std deviation: {movement_above.std():.3f}")

print(f"\n--- Testing Angular Velocity Threshold: {test_angvel_threshold} deg/frame ---")
angvel_above = calculate_time_above_threshold(angvel_snips_filtered, test_angvel_threshold)
print(f"Trials above threshold: {(angvel_above > 0.5).sum()} / {len(angvel_above)}")
print(f"Mean proportion active: {angvel_above.mean():.3f}")
print(f"Std deviation: {angvel_above.std():.3f}")

# Comparison by rat
print(f"\n--- Activity levels by rat (movement_thresh={test_movement_threshold}, angvel_thresh={test_angvel_threshold}) ---")
for rat in rats:
    rat_mask = x_filtered['rat_id'] == rat
    rat_trials = x_filtered[rat_mask].index.values
    
    rat_movement = movement_snips_filtered.loc[rat_trials]
    mov_perc = (calculate_time_above_threshold(rat_movement, test_movement_threshold) > 0.5).sum() / len(rat_trials)
    
    rat_angvel = angvel_snips_filtered.loc[rat_trials]
    angvel_perc = (calculate_time_above_threshold(rat_angvel, test_angvel_threshold) > 0.5).sum() / len(rat_trials)
    
    print(f"Rat {rat}: Movement active {mov_perc:.1%}, Angular vel active {angvel_perc:.1%}")
