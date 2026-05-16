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
# # Troubleshooting Angular Velocity Calculations
#
# This notebook helps debug and visualize angular velocity calculations using both the ears-only method and the new head_base-referenced method.
#
# Provide a DLC file stub (e.g., '2024-01-15-162145') to analyze a single rat's behavior.

# %%
# Imports
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Add src to path
sys.path.insert(0, str(Path.cwd().parent / 'src'))
from extract_behav_parameters import read_DLC_csv, interpolate_low_likehood, calc_angular_velocity, calc_angvel_earsonly

# Set plotting style
sns.set_style('whitegrid')
plt.rcParams['figure.figsize'] = (14, 10)

# %% [markdown]
# ## Configuration
#
# Edit the stub and folder paths below to specify which DLC data to analyze.

# %%
# Configuration - EDIT THIS SECTION


dlc_folder = Path("D:/TestData/bazzino/output_csv_shuffle4") #office computer, Dionne's tracking
dlc_folder = Path("C:/Users/jmc010/Data/bazzino/Output DLC shuffle 4 csv files") # laptop

stub = "PB71-221123-113609"  # Change this to the stub of the DLC file you want to analyze

# Print configuration
print(f"Analyzing stub: {stub}")
print(f"DLC folder: {dlc_folder}")

# %% [markdown]
# ## Load and Preprocess Data

# %%
# Load DLC data
print(f"Loading DLC data for stub: {stub}")
df = read_DLC_csv(stub, dlc_folder)

if df is None:
    print("Error: Could not load DLC data. Check the stub and folder path.")
else:
    print(f"Loaded DLC data with shape: {df.shape}")
    print(f"\nColumns: {df.columns.tolist()}")
    print(f"\nFirst few rows:")
    print(df.head())

# %%
# Interpolate low-likelihood detections
print("Interpolating low-likelihood detections...")
df_interp = interpolate_low_likehood(df, threshold=0.6)
print("Done.")

# %% [markdown]
# ## Calculate Angular Velocity - Both Methods

# %%
# Calculate angular velocity using both methods
print("Calculating angular velocity (ears-only method)...")
df_earsonly = calc_angvel_earsonly(df_interp, rightear="r_ear", leftear="l_ear")
print("Done.")

print("Calculating angular velocity (head_base reference method)...")
df_headbase = calc_angular_velocity(df_interp, rightear="r_ear", leftear="l_ear")
print("Done.")

# Check for required columns
print(f"\nEars-only method columns added: {[c for c in df_earsonly.columns if c not in df_interp.columns]}")
print(f"Head_base method columns added: {[c for c in df_headbase.columns if c not in df_interp.columns]}")

# %% [markdown]
# ## Inspect Geometry Data

# %%
# Check if required bodyparts are present
required_parts = ['r_ear', 'l_ear', 'head_base']

print("Checking for required bodyparts...\n")
for part in required_parts:
    x_col = f"{part}_x"
    y_col = f"{part}_y"
    
    if x_col in df_interp.columns and y_col in df_interp.columns:
        x_vals = pd.to_numeric(df_interp[x_col], errors='coerce')
        y_vals = pd.to_numeric(df_interp[y_col], errors='coerce')
        print(f"{part}:")
        print(f"  X range: {x_vals.min():.1f} - {x_vals.max():.1f}")
        print(f"  Y range: {y_vals.min():.1f} - {y_vals.max():.1f}")
        print(f"  NaN values: {x_vals.isna().sum()} / {len(x_vals)}")
        print()
    else:
        print(f"{part}: NOT FOUND in data")
        print()

# %%
# Check ear distance distribution
ear_dist = pd.to_numeric(df_headbase['ear_distance'], errors='coerce')

print("Ear Distance Statistics:")
print(ear_dist.describe())
print(f"\nFrames with ear_distance >= 90: {(ear_dist >= 90).sum()} / {len(ear_dist)}")

# Visualize ear distance
fig, ax = plt.subplots(1, 1, figsize=(12, 4))
ax.plot(ear_dist, linewidth=0.5, alpha=0.7)
ax.axhline(y=90, color='r', linestyle='--', label='Threshold (90 px)')
ax.set_xlabel('Frame')
ax.set_ylabel('Ear Distance (pixels)')
ax.set_title(f'Ear Distance Over Time - {stub}')
ax.legend()
plt.tight_layout()
plt.show()

# Histogram
fig, ax = plt.subplots(1, 1, figsize=(8, 4))
ax.hist(ear_dist.dropna(), bins=50, edgecolor='black', alpha=0.7)
ax.axvline(x=90, color='r', linestyle='--', linewidth=2, label='Threshold (90 px)')
ax.set_xlabel('Ear Distance (pixels)')
ax.set_ylabel('Frequency')
ax.set_title(f'Ear Distance Distribution')
ax.legend()
plt.tight_layout()
plt.show()

# %% [markdown]
# ## Compare Angular Velocity Outputs

# %%
# Get angular velocity columns
angvel_earsonly = pd.to_numeric(df_earsonly['d_angle_deg'], errors='coerce')
angvel_headbase = pd.to_numeric(df_headbase['d_angle_deg'], errors='coerce')

# Print statistics
print("=" * 60)
print("EARS-ONLY METHOD (d_angle_deg)")
print("=" * 60)
print(angvel_earsonly.describe())
print(f"NaN values: {angvel_earsonly.isna().sum()} / {len(angvel_earsonly)}")

print("\n" + "=" * 60)
print("HEAD_BASE REFERENCE METHOD (d_angle_deg)")
print("=" * 60)
print(angvel_headbase.describe())
print(f"NaN values: {angvel_headbase.isna().sum()} / {len(angvel_headbase)}")

# %% [markdown]
# ## Visualization: Comparison Plot

# %%
# Create comparison plot
fig, axes = plt.subplots(3, 1, figsize=(14, 10))

# Plot 1: Time series overlay
ax = axes[0]
ax.plot(angvel_earsonly, label='Ears-only method', linewidth=0.8, alpha=0.7)
ax.plot(angvel_headbase, label='Head_base reference', linewidth=0.8, alpha=0.7)
ax.set_xlabel('Frame')
ax.set_ylabel('Angular Velocity (deg/frame)')
ax.set_title(f'Angular Velocity Comparison - {stub}')
ax.legend()
ax.grid(True, alpha=0.3)

# Plot 2: Distribution comparison
ax = axes[1]
ax.hist(angvel_earsonly.dropna(), bins=50, alpha=0.6, label='Ears-only', edgecolor='black')
ax.hist(angvel_headbase.dropna(), bins=50, alpha=0.6, label='Head_base', edgecolor='black')
ax.set_xlabel('Angular Velocity (deg/frame)')
ax.set_ylabel('Frequency')
ax.set_title('Distribution Comparison')
ax.legend()
ax.grid(True, alpha=0.3)

# Plot 3: Scatter comparison
ax = axes[2]
# Remove NaN values for scatter plot
valid_idx = ~(angvel_earsonly.isna() | angvel_headbase.isna())
ax.scatter(angvel_earsonly[valid_idx], angvel_headbase[valid_idx], alpha=0.3, s=5)
# Add diagonal line
lim = [min(angvel_earsonly.min(), angvel_headbase.min()), 
       max(angvel_earsonly.max(), angvel_headbase.max())]
ax.plot(lim, lim, 'r--', linewidth=2, label='Perfect agreement')
ax.set_xlabel('Ears-only method (deg/frame)')
ax.set_ylabel('Head_base reference (deg/frame)')
ax.set_title('Scatter: Method Comparison')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

# %% [markdown]
# ## Zoom into a Time Window

# %%
# Select a time window to examine in detail
start_frame = 1000  # Change these values to zoom into different parts
end_frame = 2000

window_earsonly = angvel_earsonly.iloc[start_frame:end_frame]
window_headbase = angvel_headbase.iloc[start_frame:end_frame]

fig, ax = plt.subplots(figsize=(14, 5))
ax.plot(range(start_frame, end_frame), window_earsonly, label='Ears-only', marker='o', markersize=2, linewidth=1)
ax.plot(range(start_frame, end_frame), window_headbase, label='Head_base', marker='s', markersize=2, linewidth=1)
ax.set_xlabel('Frame')
ax.set_ylabel('Angular Velocity (deg/frame)')
ax.set_title(f'Zoomed View: Frames {start_frame} - {end_frame}')
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

print(f"\nFrames {start_frame}-{end_frame} Statistics:")
print(f"Ears-only: mean={window_earsonly.mean():.2f}, std={window_earsonly.std():.2f}")
print(f"Head_base: mean={window_headbase.mean():.2f}, std={window_headbase.std():.2f}")

# %% [markdown]
# ## Directional Analysis

# %% [markdown]
# ## Smoothing Analysis: Position Data Before vs After Angular Velocity
#
# **Question 1: Should you smooth before calculating angular velocity?**
#
# Generally **yes** - smoothing position data (x, y coordinates) *before* angle calculation is more principled than smoothing after. Here's why:
# - DLC tracking has inherent jitter in the detected positions
# - Jitter in positions → noisier angle calculations
# - Smoothing positions first reduces this jitter, leading to cleaner angles
# - This is especially important for derivatives (velocity/acceleration) which amplify noise
#
# Let's compare approaches:

# %%
from scipy.ndimage import gaussian_filter1d

# Create smoothed version of position data
df_smoothed = df_interp.copy()

# Smooth the x,y position columns using Gaussian filter
smooth_sigma = 2.0  # Adjust this to control smoothing strength (higher = more smoothing)

bodyparts_to_smooth = ['r_ear', 'l_ear', 'head_base']

for bp in bodyparts_to_smooth:
    x_col = f"{bp}_x"
    y_col = f"{bp}_y"
    
    if x_col in df_smoothed.columns and y_col in df_smoothed.columns:
        x_vals = pd.to_numeric(df_smoothed[x_col], errors='coerce')
        y_vals = pd.to_numeric(df_smoothed[y_col], errors='coerce')
        
        # Only smooth non-NaN values
        if x_vals.notna().any():
            # Create a mask for NaN values
            x_mask = x_vals.notna()
            y_mask = y_vals.notna()
            
            # Smooth non-NaN values
            x_smooth = x_vals.copy()
            y_smooth = y_vals.copy()
            
            x_smooth[x_mask] = gaussian_filter1d(x_vals[x_mask], sigma=smooth_sigma)
            y_smooth[y_mask] = gaussian_filter1d(y_vals[y_mask], sigma=smooth_sigma)
            
            df_smoothed[x_col] = x_smooth
            df_smoothed[y_col] = y_smooth

print(f"Smoothed position data with Gaussian filter (sigma={smooth_sigma})")")

# %%
# Calculate angular velocity from smoothed position data
print("Calculating angular velocity from smoothed position data...")
df_smoothed_angvel = calc_angular_velocity(df_smoothed, rightear="r_ear", leftear="l_ear")
angvel_smoothed = pd.to_numeric(df_smoothed_angvel['d_angle_deg'], errors='coerce')
print("Done.")

# Visualization: Effect of smoothing on angular velocity
fig, axes = plt.subplots(3, 1, figsize=(14, 10))

# Plot 1: Original vs smoothed angular velocity (head_base method)
ax = axes[0]
ax.plot(angvel_headbase, label='Angular velocity (raw positions)', linewidth=0.8, alpha=0.7)
ax.plot(angvel_smoothed, label='Angular velocity (smoothed positions)', linewidth=0.8, alpha=0.7)
ax.set_xlabel('Frame')
ax.set_ylabel('Angular Velocity (deg/frame)')
ax.set_title(f'Effect of Smoothing Position Data Before Angle Calculation - {stub}')
ax.legend()
ax.grid(True, alpha=0.3)

# Plot 2: Zoomed comparison
ax = axes[1]
start, end = 1000, 2000
ax.plot(range(start, end), angvel_headbase.iloc[start:end], label='Raw positions', marker='o', markersize=2, linewidth=1)
ax.plot(range(start, end), angvel_smoothed.iloc[start:end], label='Smoothed positions', marker='s', markersize=2, linewidth=1)
ax.set_xlabel('Frame')
ax.set_ylabel('Angular Velocity (deg/frame)')
ax.set_title(f'Zoomed View: Frames {start}-{end}')
ax.legend()
ax.grid(True, alpha=0.3)

# Plot 3: Distribution comparison
ax = axes[2]
ax.hist(angvel_headbase.dropna(), bins=50, alpha=0.6, label='Raw positions', edgecolor='black')
ax.hist(angvel_smoothed.dropna(), bins=50, alpha=0.6, label='Smoothed positions', edgecolor='black')
ax.set_xlabel('Angular Velocity (deg/frame)')
ax.set_ylabel('Frequency')
ax.set_title('Distribution Comparison')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

# Statistics
print("\n" + "="*60)
print("EFFECT OF SMOOTHING POSITION DATA")
print("="*60)
print(f"\nRaw position data - Angular velocity stats:")
print(f"  Mean: {angvel_headbase.mean():.3f}, Std: {angvel_headbase.std():.3f}")
print(f"  Min: {angvel_headbase.min():.3f}, Max: {angvel_headbase.max():.3f}")

print(f"\nSmoothed position data - Angular velocity stats:")
print(f"  Mean: {angvel_smoothed.mean():.3f}, Std: {angvel_smoothed.std():.3f}")
print(f"  Min: {angvel_smoothed.min():.3f}, Max: {angvel_smoothed.max():.3f}")

print(f"\nSmoothing reduced angular velocity std dev by {(1 - angvel_smoothed.std()/angvel_headbase.std())*100:.1f}%")

# %% [markdown]
# ## Absolute Value Analysis: Direction vs Magnitude
#
# **Question 2: Should you take the absolute value of angular velocity?**
#
# Use absolute values when you care about **how much** the head moves (magnitude), not **which direction** it rotates. Key considerations:
# - **Signed angular velocity**: Preserves direction information (positive=clockwise, negative=counterclockwise)
# - **Absolute angular velocity**: Only tells you speed of rotation, useful for measuring overall motor behavior
# - **Use cases for absolute value**: Relating to dopamine activity, motor engagement, exploratory behavior
# - **Use cases for signed**: Understanding turning preferences, directional bias, head checking patterns

# %%
# Calculate absolute values
angvel_headbase_abs = np.abs(angvel_headbase)
angvel_smoothed_abs = np.abs(angvel_smoothed)

# Visualization: Signed vs Absolute angular velocity
fig, axes = plt.subplots(3, 2, figsize=(15, 10))

# Row 1: Signed angular velocity (with direction)
ax = axes[0, 0]
colors = np.sign(angvel_headbase.fillna(0))
scatter = ax.scatter(range(len(angvel_headbase)), angvel_headbase, c=colors, cmap='RdBu', s=3, alpha=0.6)
ax.axhline(y=0, color='k', linestyle='-', linewidth=0.5)
ax.set_ylabel('Angular Velocity (deg/frame)')
ax.set_title('Signed Angular Velocity (headbase, raw data)\nred=clockwise, blue=counterclockwise')
ax.grid(True, alpha=0.3)

ax = axes[0, 1]
colors = np.sign(angvel_smoothed.fillna(0))
scatter = ax.scatter(range(len(angvel_smoothed)), angvel_smoothed, c=colors, cmap='RdBu', s=3, alpha=0.6)
ax.axhline(y=0, color='k', linestyle='-', linewidth=0.5)
ax.set_ylabel('Angular Velocity (deg/frame)')
ax.set_title('Signed Angular Velocity (headbase, smoothed data)')
ax.grid(True, alpha=0.3)

# Row 2: Absolute angular velocity (magnitude only)
ax = axes[1, 0]
ax.plot(angvel_headbase_abs, linewidth=0.8, alpha=0.7, color='green')
ax.fill_between(range(len(angvel_headbase_abs)), 0, angvel_headbase_abs, alpha=0.3, color='green')
ax.set_ylabel('|Angular Velocity| (deg/frame)')
ax.set_title('Absolute Angular Velocity (headbase, raw data)\nMagnitude of rotation (no direction)')
ax.grid(True, alpha=0.3)

ax = axes[1, 1]
ax.plot(angvel_smoothed_abs, linewidth=0.8, alpha=0.7, color='green')
ax.fill_between(range(len(angvel_smoothed_abs)), 0, angvel_smoothed_abs, alpha=0.3, color='green')
ax.set_ylabel('|Angular Velocity| (deg/frame)')
ax.set_title('Absolute Angular Velocity (headbase, smoothed data)')
ax.grid(True, alpha=0.3)

# Row 3: Distribution comparison
ax = axes[2, 0]
ax.hist(angvel_headbase.dropna(), bins=50, alpha=0.6, label='Signed', edgecolor='black', color='steelblue')
ax.hist(angvel_headbase_abs.dropna(), bins=50, alpha=0.6, label='Absolute', edgecolor='black', color='green')
ax.set_xlabel('Angular Velocity (deg/frame)')
ax.set_ylabel('Frequency')
ax.set_title('Distribution: Signed vs Absolute (raw data)')
ax.legend()
ax.grid(True, alpha=0.3)

ax = axes[2, 1]
ax.hist(angvel_smoothed.dropna(), bins=50, alpha=0.6, label='Signed', edgecolor='black', color='steelblue')
ax.hist(angvel_smoothed_abs.dropna(), bins=50, alpha=0.6, label='Absolute', edgecolor='black', color='green')
ax.set_xlabel('Angular Velocity (deg/frame)')
ax.set_ylabel('Frequency')
ax.set_title('Distribution: Signed vs Absolute (smoothed data)')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

# Statistics comparison
print("\n" + "="*70)
print("SIGNED vs ABSOLUTE ANGULAR VELOCITY")
print("="*70)

print("\n--- RAW POSITION DATA ---")
print(f"\nSigned (includes direction):")
print(f"  Mean: {angvel_headbase.mean():.3f}, Std: {angvel_headbase.std():.3f}")
print(f"  Min: {angvel_headbase.min():.3f}, Max: {angvel_headbase.max():.3f}")
print(f"  Frames moving right (positive): {(angvel_headbase > 0).sum()}")
print(f"  Frames moving left (negative): {(angvel_headbase < 0).sum()}")

print(f"\nAbsolute (magnitude only):")
print(f"  Mean: {angvel_headbase_abs.mean():.3f}, Std: {angvel_headbase_abs.std():.3f}")
print(f"  Min: {angvel_headbase_abs.min():.3f}, Max: {angvel_headbase_abs.max():.3f}")

print("\n--- SMOOTHED POSITION DATA ---")
print(f"\nSigned (includes direction):")
print(f"  Mean: {angvel_smoothed.mean():.3f}, Std: {angvel_smoothed.std():.3f}")
print(f"  Min: {angvel_smoothed.min():.3f}, Max: {angvel_smoothed.max():.3f}")
print(f"  Frames moving right (positive): {(angvel_smoothed > 0).sum()}")
print(f"  Frames moving left (negative): {(angvel_smoothed < 0).sum()}")

print(f"\nAbsolute (magnitude only):")
print(f"  Mean: {angvel_smoothed_abs.mean():.3f}, Std: {angvel_smoothed_abs.std():.3f}")
print(f"  Min: {angvel_smoothed_abs.min():.3f}, Max: {angvel_smoothed_abs.max():.3f}")

# %% [markdown]
# ## Summary & Recommendations
#
# Use this table to decide which version to use for your analysis:
#
# | Metric | When to Use | Good For |
# |--------|------------|----------|
# | **Signed + Raw** | Want raw data without any filtering | Debugging, understanding jitter |
# | **Signed + Smoothed** | Want to preserve direction with cleaner signal | Directional bias analysis, turning patterns |
# | **Absolute + Raw** | Want magnitude only, but can tolerate noise | Quick exploratory analysis |
# | **Absolute + Smoothed** (✓ Recommended) | Want clean magnitude signal | Relating to dopamine, motor engagement, activity levels |
#
# **My recommendation for relating to dopamine/neural activity**: Use **absolute value of angular velocity from smoothed position data** - you get a clean signal of head movement magnitude without directional artifacts.
#
# Consider adjusting `smooth_sigma` (currently set to 2.0) to control smoothing strength:
# - Lower values (0.5-1.0): More detail, some jitter remains
# - Medium values (1.5-2.5): Good balance of smoothing and detail
# - Higher values (3.0+): Very smooth, but may lose fast movements

# %%
# Analyze angular displacement (absolute value) vs direction
fig, axes = plt.subplots(2, 2, figsize=(14, 8))

# Ears-only: magnitude and direction
ax = axes[0, 0]
ax.scatter(range(len(angvel_earsonly)), angvel_earsonly, c=np.sign(angvel_earsonly.fillna(0)), 
          cmap='RdBu', s=3, alpha=0.6)
ax.set_xlabel('Frame')
ax.set_ylabel('Angular Velocity (deg/frame)')
ax.set_title('Ears-only: Direction (red=right, blue=left)')
ax.axhline(y=0, color='k', linestyle='-', linewidth=0.5)
ax.grid(True, alpha=0.3)

# Head_base: magnitude and direction
ax = axes[0, 1]
ax.scatter(range(len(angvel_headbase)), angvel_headbase, c=np.sign(angvel_headbase.fillna(0)), 
          cmap='RdBu', s=3, alpha=0.6)
ax.set_xlabel('Frame')
ax.set_ylabel('Angular Velocity (deg/frame)')
ax.set_title('Head_base: Direction (red=right, blue=left)')
ax.axhline(y=0, color='k', linestyle='-', linewidth=0.5)
ax.grid(True, alpha=0.3)

# Cumulative angle over time
ax = axes[1, 0]
cumsum_earsonly = angvel_earsonly.fillna(0).cumsum()
cumsum_headbase = angvel_headbase.fillna(0).cumsum()
ax.plot(cumsum_earsonly, label='Ears-only', linewidth=1)
ax.plot(cumsum_headbase, label='Head_base', linewidth=1)
ax.set_xlabel('Frame')
ax.set_ylabel('Cumulative Angular Displacement (degrees)')
ax.set_title('Cumulative Rotation Over Time')
ax.legend()
ax.grid(True, alpha=0.3)

# Raw angles
ax = axes[1, 1]
angle_earsonly = pd.to_numeric(df_earsonly['angle_rad'], errors='coerce')
angle_headbase = pd.to_numeric(df_headbase['angle_rad'], errors='coerce')
ax.plot(np.rad2deg(angle_earsonly), label='Ears-only', linewidth=0.5, alpha=0.7)
ax.plot(np.rad2deg(angle_headbase), label='Head_base', linewidth=0.5, alpha=0.7)
ax.set_xlabel('Frame')
ax.set_ylabel('Angle (degrees)')
ax.set_title('Absolute Head Angle Over Time')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

# %% [markdown]
# ## Export Data for Further Analysis

# %%
# Option to export the processed data
export_folder = Path.cwd() / 'troubleshooting_outputs'
export_folder.mkdir(exist_ok=True)

# Export both dataframes with angular velocity calculations
output_earsonly = export_folder / f'{stub}_angvel_earsonly.csv'
output_headbase = export_folder / f'{stub}_angvel_headbase.csv'

# Select relevant columns to export
export_cols_earsonly = [
    'r_ear_x', 'r_ear_y', 'l_ear_x', 'l_ear_y',
    'ear_distance', 'angle_rad', 'd_angle_deg'
]

export_cols_headbase = [
    'r_ear_x', 'r_ear_y', 'l_ear_x', 'l_ear_y', 'head_base_x', 'head_base_y',
    'ear_distance', 'ear_midpoint_x', 'ear_midpoint_y', 'rel_head_x', 'rel_head_y',
    'angle_rad', 'd_angle_deg'
]

# Check which columns exist and filter
export_cols_earsonly = [c for c in export_cols_earsonly if c in df_earsonly.columns]
export_cols_headbase = [c for c in export_cols_headbase if c in df_headbase.columns]

df_earsonly[export_cols_earsonly].to_csv(output_earsonly, index=False)
df_headbase[export_cols_headbase].to_csv(output_headbase, index=False)

print(f"Exported ears-only results to: {output_earsonly}")
print(f"Exported head_base results to: {output_headbase}")
