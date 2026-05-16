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
# # Troubleshoot SIMBA CI Threshold
#
# Diagnostic plots to understand `simba_pct_time_above_95ci`:
#
# 1. **Strip/box plot** of pct_time values by the four groups (from assembled pickle)
# 2. **Session-level: snips vs CI band** — for a chosen session, re-run the shuffle and overlay
#    individual real trial snips on the per-bin CI band
# 3. **Histogram**: distribution of shuffled bin values vs real bin values in the infusion window
# 4. **Per-animal summary** of pct_time values by condition

# %%
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import dill

sys.path.insert(0, str(Path('../src').resolve()))
from pickle_compat import enable_dill_pathlib_compat
enable_dill_pathlib_compat()

from extract_simba import (
    read_simba_probability,
    get_cam1_onsets_for_stub,
    make_simba_snips,
    get_shifted_snip_means,
    baseline_simba_snips,
    get_time_above_simba_ci,
)
from extract_behav_parameters import get_ttls

from figure_config import DATAFOLDER, COLORS

sns.set_style('whitegrid')
plt.rcParams.update({'figure.dpi': 100, 'font.size': 10})
colors = COLORS

# %% [markdown]
# ## Load assembled data

# %%
assembled_data_path = DATAFOLDER / 'assembled_data.pickle'
with open(assembled_data_path, 'rb') as f:
    data = dill.load(f)

x_array    = data['x_array']
snips_simba = data['snips_simba']
params     = data.get('params', {})
metadata   = data.get('metadata', {})

print('Loaded:', assembled_data_path)
print('x_array shape:', x_array.shape)
print('snips_simba shape:', snips_simba.shape)
print('simba_ci_percentile in metadata:', metadata.get('simba_ci_percentile', 'NOT FOUND'))
print()
print('pct_time_above_95ci summary by group:')
print(
    x_array.groupby(['condition', 'infusiontype'])['simba_pct_time_above_95ci']
    .describe().round(2)
)

# %% [markdown]
# ## Plot 1: Strip plot of pct_time by group
#
# Per-animal mean (x_array already has one row per trial; we average within animal first).

# %%
group_order = [
    ('replete', '10NaCl'),
    ('replete', '45NaCl'),
    ('deplete', '10NaCl'),
    ('deplete', '45NaCl'),
]
group_labels = ['Rep 10', 'Rep 45', 'Dep 10', 'Dep 45']
group_colors = colors[:4]

# Per-animal averages
per_animal = (
    x_array
    .groupby(['id', 'condition', 'infusiontype'])
    ['simba_pct_time_above_95ci']
    .mean()
    .reset_index()
)

fig, axes = plt.subplots(1, 2, figsize=(9, 4), sharey=False)

# Left: per-trial values (all trials)
ax = axes[0]
for i, ((cond, inf), label, col) in enumerate(zip(group_order, group_labels, group_colors)):
    vals = x_array.query('condition==@cond and infusiontype==@inf')['simba_pct_time_above_95ci'].dropna()
    ax.scatter(np.full(len(vals), i) + np.random.uniform(-0.15, 0.15, len(vals)),
               vals, alpha=0.3, s=8, color=col)
    ax.errorbar(i, vals.mean(), yerr=vals.sem(), fmt='o', color=col,
                markersize=7, capsize=4, linewidth=1.5, zorder=5)

ax.axhline(0, color='k', lw=0.8, ls='--', alpha=0.5)
ax.set_xticks(range(4)); ax.set_xticklabels(group_labels)
ax.set_ylabel('% time above CI')
ax.set_title('Per-trial values')

# Right: per-animal means
ax = axes[1]
for i, ((cond, inf), label, col) in enumerate(zip(group_order, group_labels, group_colors)):
    vals = per_animal.query('condition==@cond and infusiontype==@inf')['simba_pct_time_above_95ci'].dropna()
    ax.scatter(np.full(len(vals), i) + np.random.uniform(-0.15, 0.15, len(vals)),
               vals, alpha=0.7, s=30, color=col, zorder=4)
    ax.errorbar(i, vals.mean(), yerr=vals.sem(), fmt='o', color=col,
                markersize=9, capsize=5, linewidth=2, zorder=5)

ax.axhline(0, color='k', lw=0.8, ls='--', alpha=0.5)
ax.set_xticks(range(4)); ax.set_xticklabels(group_labels)
ax.set_ylabel('% time above CI (animal mean)')
ax.set_title('Per-animal means')

fig.tight_layout()
plt.show()

# Print per-animal table
print('Per-animal means:')
display(
    per_animal
    .pivot_table(index='id', columns=['condition','infusiontype'],
                 values='simba_pct_time_above_95ci', aggfunc='mean')
    .round(1)
)

# %% [markdown]
# ## Plot 2: Example session — individual snips vs per-bin CI band
#
# Re-runs the shuffle for one session and plots:
# - All real trial snips (thin translucent lines)
# - Mean real snip (thick line)
# - Per-bin CI band from the shuffled null distribution (5th–95th percentile shaded)
# - The per-bin threshold used for pct_time (red dashed line = `simba_ci_percentile`)
#
# **Edit `STUB` and `CONDITION` to choose the session to inspect.**

# %%
# ─── EDIT HERE ────────────────────────────────────────────────────────
# x_array['id'] contains subject IDs (e.g. "PB23").
# Simba files use the full session stub (e.g. "PB23-220608-131619") from the FileKey.
# Run this cell to see the mapping before choosing a subject below.
data_folder_temp = Path('../data')
file_key = pd.concat([
    pd.read_csv(data_folder_temp / '10NaCl_FileKey.csv'),
    pd.read_csv(data_folder_temp / '45NaCl_FileKey.csv'),
]).drop_duplicates(subset='Subject')

print('Subject → Session stub mapping:')
display(file_key[['Subject', 'Folder', 'Physiological state', 'TreatNum']].sort_values('Subject').reset_index(drop=True))

# %%
# Set SUBJECT_ID to any value shown in the cell above, e.g. "PB23"
SUBJECT_ID = sorted(x_array['id'].unique())[0]   # ← change this
SUBJECT_ID = "PB24"
# Look up the full session stub from the FileKey
stub_lookup = file_key.set_index('Subject')['Folder'].to_dict()
if SUBJECT_ID not in stub_lookup:
    raise ValueError(f"'{SUBJECT_ID}' not in FileKey. Available: {sorted(stub_lookup.keys())}")
STUB = stub_lookup[SUBJECT_ID]

print(f'Subject ID  : {SUBJECT_ID}')
print(f'Session stub: {STUB}')

# Path params: params dict in the pickle may contain None or relative paths.
# _resolve_param_path falls back to the hardcoded default if the pickled value
# is None or points to a path that does not exist on this machine.
def _resolve_param_path(params, key, fallback):
    v = params.get(key)
    if v is None:
        return Path(fallback)
    p = Path(v)
    return p if p.exists() else Path(fallback)

simba_folder  = _resolve_param_path(params, 'simba_folder', 'C:/Users/jmc010/Data/bazzino/simba')
tank_folder   = _resolve_param_path(params, 'tank_folder',  'C:/Users/jmc010/Data/bazzino/tanks')
data_folder   = _resolve_param_path(params, 'data_folder',  '../data')

prob_col      = params.get('simba_probability_column', 'Probability_Appetitive')
fps           = params.get('simba_fps', 10)
pre_bins      = params.get('simba_pre_bins', 50)
post_bins     = params.get('simba_post_bins', 150)
shift_frames  = params.get('simba_shift_frames', 300)
n_shuffles    = params.get('simba_n_shuffles', 100)
ci_percentile = params.get('simba_ci_percentile', 95)
auc_start     = params.get('auc_start_bin', 50)
auc_end       = params.get('auc_end_bin', 150)

print(f'CI percentile: {ci_percentile}')
print(f'simba_folder : {simba_folder}')
print(f'data_folder  : {data_folder}')

# Use global TTL timing stats to pick a safer shuffle shift in seconds.
INFUSION_DUR_S = 10.0
ttls_table = pd.read_csv(data_folder / 'ttls.csv')
if ttls_table.columns[0].startswith('Unnamed') or ttls_table.columns[0] == '':
    ttls_table = ttls_table.drop(columns=ttls_table.columns[0])

ttls_by_stub = {
    col: np.sort(pd.to_numeric(ttls_table[col], errors='coerce').dropna().to_numpy(dtype=float))
    for col in ttls_table.columns
}
all_isis = np.concatenate([np.diff(v) for v in ttls_by_stub.values() if len(v) > 1])
print('TTL ISI summary (s): min={:.2f}, p5={:.2f}, median={:.2f}, p95={:.2f}, max={:.2f}'.format(
    np.min(all_isis), np.percentile(all_isis, 5), np.median(all_isis), np.percentile(all_isis, 95), np.max(all_isis)
))

def _contamination_rate(shift_s, duration_s=INFUSION_DUR_S):
    total = 0
    contaminated = 0
    for vals in ttls_by_stub.values():
        if len(vals) == 0:
            continue
        total += len(vals)
        for t in vals:
            a, b = t + shift_s, t + shift_s + duration_s
            overlap = np.any((a < vals + duration_s) & (b > vals))
            contaminated += int(overlap)
    return contaminated / total

def _global_margin(shift_s, duration_s=INFUSION_DUR_S):
    margins = []
    for vals in ttls_by_stub.values():
        for t in vals:
            a, b = t + shift_s, t + shift_s + duration_s
            # Positive margin means disjoint; negative means overlap amount.
            gaps = np.maximum(vals - b, a - (vals + duration_s))
            margins.append(np.min(gaps))
    return np.min(margins), np.percentile(margins, 5), np.median(margins)

candidate_shifts_s = np.arange(10, 121, dtype=float)
contam = np.array([_contamination_rate(s) for s in candidate_shifts_s])
safe_shifts_s = candidate_shifts_s[contam == 0]

if len(safe_shifts_s) > 0:
    safe_m5 = np.array([_global_margin(s)[1] for s in safe_shifts_s])
    recommended_shift_s = float(safe_shifts_s[np.argmax(safe_m5)])
else:
    # Conservative fallback if no zero-overlap shifts are found.
    recommended_shift_s = float(max(10.0, np.percentile(all_isis, 5) - INFUSION_DUR_S - 1.0))

current_shift_s = shift_frames / fps
print(f'Current shift from params: {current_shift_s:.1f}s ({shift_frames} frames @ {fps} fps)')
print(f'Contamination at current shift: {_contamination_rate(current_shift_s)*100:.2f}%')
print(f'Recommended shift from ttls.csv: {recommended_shift_s:.1f}s')
print(f'Contamination at recommended shift: {_contamination_rate(recommended_shift_s)*100:.2f}%')
m_min, m_p5, m_med = _global_margin(recommended_shift_s)
print(f'Margin to nearest infusion at recommended shift (min/p5/median): {m_min:.2f}/{m_p5:.2f}/{m_med:.2f}s')

# Apply recommendation for this diagnostic run.
shift_seconds = recommended_shift_s
shift_frames = int(round(shift_seconds * fps))
print(f'Using shift = {shift_seconds:.1f}s ({shift_frames} frames) for shuffle null')

# Load probability vector using session stub
prob_vec, src_file = read_simba_probability(STUB, simba_folder, probability_column=prob_col)
print(f'Simba file: {src_file.name}, length: {len(prob_vec)} frames')

# Load TTLs and cam1 onsets
solenoid_ts = get_ttls(STUB, data_folder)
try:
    cam1_onsets = get_cam1_onsets_for_stub(STUB, tank_folder)
    alignment = 'cam1'
except Exception as e:
    cam1_onsets = None
    alignment = 'fps'
    print(f'Cam1 unavailable ({e}), using fps alignment')
print(f'Alignment: {alignment}, n_trials (solenoid TTLs): {len(solenoid_ts)-1}')

# Build real snips
real_snips = make_simba_snips(prob_vec, solenoid_ts, fps=fps,
                               pre_bins=pre_bins, post_bins=post_bins,
                               cam1_onsets=cam1_onsets)
print(f'Real snips shape: {real_snips.shape}')

# Build shuffled null distribution (all individual snips stacked)
print(f'Running {n_shuffles} shuffles...')
shuffled_snips = get_shifted_snip_means(prob_vec, solenoid_ts, fps=fps,
                                         pre_bins=pre_bins, post_bins=post_bins,
                                         shift_frames=shift_frames,
                                         n_shuffles=n_shuffles,
                                         cam1_onsets=cam1_onsets)
print(f'Shuffled snips stacked shape: {shuffled_snips.shape}')

# Compute pct_time for this session
pct_time, mean_threshold = get_time_above_simba_ci(
    real_snips, shuffled_snips,
    percentile=ci_percentile,
    start_bin=auc_start, end_bin=auc_end
)
print(f'pct_time (per trial): {np.round(pct_time * 100, 1)}')
print(f'Mean CI threshold across infusion window: {mean_threshold:.4f}')

# %%
n_bins = pre_bins + post_bins
time_axis = np.arange(n_bins) / fps - pre_bins / fps  # seconds relative to infusion

# Per-bin CI from shuffled snips (full time axis)
ci_lo  = np.nanpercentile(shuffled_snips, 100 - ci_percentile, axis=0)
ci_hi  = np.nanpercentile(shuffled_snips, ci_percentile, axis=0)
shuf_median = np.nanmedian(shuffled_snips, axis=0)

# Also get the infusion-window-only threshold line for context
threshold_per_bin_window = np.nanpercentile(
    shuffled_snips[:, auc_start:auc_end], ci_percentile, axis=0
)  # shape: (auc_end - auc_start,)
threshold_line = np.full(n_bins, np.nan)
threshold_line[auc_start:auc_end] = threshold_per_bin_window

fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)

for ax, show_baseline in zip(axes, [False, True]):
    if show_baseline:
        baseline_val = float(np.nanmean(shuffled_snips))
        plot_snips   = real_snips - baseline_val
        plot_ci_lo   = ci_lo - baseline_val
        plot_ci_hi   = ci_hi - baseline_val
        plot_shuf    = shuf_median - baseline_val
        plot_thresh  = threshold_line - baseline_val
        title_suffix = '(baseline-subtracted)'
    else:
        plot_snips  = real_snips
        plot_ci_lo  = ci_lo
        plot_ci_hi  = ci_hi
        plot_shuf   = shuf_median
        plot_thresh = threshold_line
        title_suffix = '(raw probability)'

    # Individual real trial snips
    for snip in plot_snips:
        ax.plot(time_axis, snip, color='steelblue', alpha=0.2, lw=0.8)

    # Mean real snip
    ax.plot(time_axis, np.nanmean(plot_snips, axis=0),
            color='steelblue', lw=2.5, label='Mean real')

    # Shuffled CI band
    ax.fill_between(time_axis, plot_ci_lo, plot_ci_hi,
                    color='grey', alpha=0.25, label=f'Shuffle {100-ci_percentile:.0f}–{ci_percentile:.0f}% CI')
    ax.plot(time_axis, plot_shuf, color='grey', lw=1.2, ls='--', label='Shuffle median')

    # Per-bin CI threshold (infusion window only)
    ax.plot(time_axis, plot_thresh, color='red', lw=1.5, ls=':', label=f'{ci_percentile}th pct threshold')

    ax.axvspan(0, 10, color='k', alpha=0.05, zorder=-10)  # infusion window
    ax.axvline(0, color='k', lw=0.8, ls='--')
    ax.axhline(0, color='k', lw=0.5, ls=':')
    ax.set_xlabel('Time from infusion (s)')
    ax.set_ylabel('Appetitive probability')
    ax.set_title(f'{STUB}  {title_suffix}')
    ax.legend(fontsize=8)

fig.suptitle(f'SIMBA CI diagnostics — {STUB}', fontsize=12)
fig.tight_layout()
plt.show()

print(f'\nPer-trial % time above {ci_percentile}th pct CI (infusion window):')
for i, p in enumerate(pct_time):
    print(f'  Trial {i:2d}: {p*100:.1f}%')

# %% [markdown]
# ## Plot 2b: Does smoothing improve real-vs-shuffled separation?
#
# Runs the same TTL-aligned analysis on smoothed versions of the SIMBA probability vector.
#
# - Uses centered moving-average smoothing (window in frames)
# - Recomputes real and shuffled snips for each smoothing level
# - Compares separation metrics in the infusion window

# %%
# Smoothing windows in frames; 1 means no smoothing (raw).
smooth_windows = [1, 3, 5, 9]  # at 10 fps: 0.1s, 0.3s, 0.5s, 0.9s

def smooth_centered_moving_average(x, window):
    x = np.asarray(x, dtype=float)
    if window <= 1:
        return x.copy()
    k = np.ones(int(window), dtype=float) / float(window)
    return np.convolve(x, k, mode='same')

def evaluate_signal_separation(probability_vector, smooth_window):
    p = smooth_centered_moving_average(probability_vector, smooth_window)

    rs = make_simba_snips(
        p, solenoid_ts, fps=fps, pre_bins=pre_bins, post_bins=post_bins, cam1_onsets=cam1_onsets
    )
    ss = get_shifted_snip_means(
        p, solenoid_ts, fps=fps, pre_bins=pre_bins, post_bins=post_bins,
        shift_frames=shift_frames, n_shuffles=n_shuffles, cam1_onsets=cam1_onsets
    )

    pct, _ = get_time_above_simba_ci(
        rs, ss, percentile=ci_percentile, start_bin=auc_start, end_bin=auc_end
    )

    rw = rs[:, auc_start:auc_end].reshape(-1)
    sw = ss[:, auc_start:auc_end].reshape(-1)
    rw = rw[~np.isnan(rw)]
    sw = sw[~np.isnan(sw)]

    threshold = float(np.nanpercentile(sw, ci_percentile))
    frac_real_above = float(np.mean(rw > threshold)) if len(rw) else np.nan
    effect_z = (float(np.mean(rw)) - float(np.mean(sw))) / float(np.std(sw) + 1e-9)

    ci_lo_w = np.nanpercentile(ss[:, auc_start:auc_end], 100 - ci_percentile, axis=0)
    ci_hi_w = np.nanpercentile(ss[:, auc_start:auc_end], ci_percentile, axis=0)
    ci_band_width = float(np.nanmedian(ci_hi_w - ci_lo_w))

    return {
        'smooth_window_frames': int(smooth_window),
        'smooth_window_seconds': float(smooth_window / fps),
        'pct_time_trial_mean': float(np.nanmean(pct)),
        'pct_time_trial_median': float(np.nanmedian(pct)),
        'real_above_threshold_frac': frac_real_above,
        'effect_z_vs_shuffle': effect_z,
        'median_ci_band_width': ci_band_width,
        'real_snips': rs,
        'shuffled_snips': ss,
    }

results = [evaluate_signal_separation(prob_vec, w) for w in smooth_windows]
results_df = pd.DataFrame([{k: v for k, v in r.items() if k not in ['real_snips', 'shuffled_snips']} for r in results])

display(
    results_df.assign(
        pct_time_trial_mean=lambda d: d['pct_time_trial_mean'] * 100,
        pct_time_trial_median=lambda d: d['pct_time_trial_median'] * 100,
        real_above_threshold_frac=lambda d: d['real_above_threshold_frac'] * 100,
    ).rename(columns={
        'pct_time_trial_mean': 'pct_time_mean_%',
        'pct_time_trial_median': 'pct_time_median_%',
        'real_above_threshold_frac': f'real_bins_above_{ci_percentile}pct_%',
    }).round(3)
 )

# Visual comparison: raw vs best smoothing by effect size
best_idx = int(np.nanargmax(results_df['effect_z_vs_shuffle'].values))
best = results[best_idx]
raw = results[0]

fig, axes = plt.subplots(1, 2, figsize=(13, 4), sharey=True)
time_axis = np.arange(pre_bins + post_bins) / fps - pre_bins / fps

for ax, label, bundle in [
    (axes[0], f'Raw (window={raw["smooth_window_frames"]} frame)', raw),
    (axes[1], f'Best smooth (window={best["smooth_window_frames"]} frames)', best),
]:
    rs = bundle['real_snips']
    ss = bundle['shuffled_snips']
    ci_lo = np.nanpercentile(ss, 100 - ci_percentile, axis=0)
    ci_hi = np.nanpercentile(ss, ci_percentile, axis=0)

    ax.plot(time_axis, np.nanmean(rs, axis=0), color='steelblue', lw=2.2, label='Mean real')
    ax.fill_between(time_axis, ci_lo, ci_hi, color='grey', alpha=0.25, label='Shuffle CI')
    ax.plot(time_axis, np.nanmedian(ss, axis=0), color='grey', ls='--', lw=1.2, label='Shuffle median')

    ax.axvspan(0, 10, color='k', alpha=0.06, zorder=-10)
    ax.axvline(0, color='k', lw=0.8, ls='--')
    ax.set_title(label)
    ax.set_xlabel('Time from infusion (s)')
    ax.set_ylabel('Appetitive probability')

axes[0].legend(fontsize=8)
fig.suptitle(f'Smoothing impact on real-vs-shuffled separation — {STUB}', fontsize=12)
fig.tight_layout()
plt.show()

print('Interpretation helper:')
print('- Larger effect_z_vs_shuffle and real_bins_above_threshold suggest better separation.')
print('- Smaller median_ci_band_width means a tighter null CI band.')

# %%
# Re-test with more shuffles and compare against the currently displayed results_df
baseline_n_shuffles = int(n_shuffles)
baseline_results_df = results_df.copy()

high_n_shuffles = 1000  # increase for a stabler null estimate
n_shuffles = high_n_shuffles

results_hi = [evaluate_signal_separation(prob_vec, w) for w in smooth_windows]
results_hi_df = pd.DataFrame([{k: v for k, v in r.items() if k not in ['real_snips', 'shuffled_snips']} for r in results_hi])

compare = baseline_results_df.merge(
    results_hi_df, on=['smooth_window_frames', 'smooth_window_seconds'], suffixes=('_base', '_hi')
)

for metric in ['pct_time_trial_mean', 'real_above_threshold_frac', 'effect_z_vs_shuffle', 'median_ci_band_width']:
    compare[f'delta_{metric}'] = compare[f'{metric}_hi'] - compare[f'{metric}_base']

print(f'Baseline n_shuffles: {baseline_n_shuffles}')
print(f'High n_shuffles: {high_n_shuffles}')

display(
    compare[[
        'smooth_window_frames',
        'pct_time_trial_mean_base', 'pct_time_trial_mean_hi', 'delta_pct_time_trial_mean',
        'real_above_threshold_frac_base', 'real_above_threshold_frac_hi', 'delta_real_above_threshold_frac',
        'effect_z_vs_shuffle_base', 'effect_z_vs_shuffle_hi', 'delta_effect_z_vs_shuffle',
        'median_ci_band_width_base', 'median_ci_band_width_hi', 'delta_median_ci_band_width',
    ]].round(4)
 )

# Keep notebook state on the high-shuffle setting for subsequent cells
results_df = results_hi_df.copy()
results = results_hi
print('Notebook variables updated to high-shuffle results (results_df, results, n_shuffles).')

# %% [markdown]
# ## Plot 2c: Null model excluding infusion windows from shifted source
#
# Tests whether the shifted null is inflated by sampling from true infusion periods.
#
# - Masks infusion windows in the source probability trace before circular shifts
# - Recomputes shuffled snips and CI threshold
# - Compares metrics vs the standard shifted null

# %%
# Compare standard shifted-null vs shifted-null with infusion windows excluded from source trace
INFUSION_DUR_S = 10.0
EXCLUSION_GUARD_S = 2.0  # extra padding on each side of infusion window

def _trial_ttls_from_solenoid(solenoid_ts):
    ts = np.asarray(solenoid_ts, dtype=float)
    if ts.size < 2:
        return np.asarray([], dtype=float)
    ts = ts[:-1]
    return ts[np.isfinite(ts)]

def _time_to_frame_idx(times_s, fps, n_frames, cam1_onsets=None):
    t = np.asarray(times_s, dtype=float)
    if cam1_onsets is None:
        idx = np.round(t * fps).astype(int)
    else:
        cam = np.asarray(cam1_onsets, dtype=float)
        idx = np.searchsorted(cam, t, side='left')
        idx = np.clip(idx, 1, cam.size - 1)
        left = idx - 1
        choose_left = np.abs(t - cam[left]) <= np.abs(cam[idx] - t)
        idx = np.where(choose_left, left, idx).astype(int)
    return np.clip(idx, 0, max(n_frames - 1, 0))

def _make_infusion_exclusion_mask(n_frames, trial_ttls, fps, infusion_dur_s=10.0, guard_s=0.0, cam1_onsets=None):
    mask = np.zeros(int(n_frames), dtype=bool)
    if n_frames <= 0 or len(trial_ttls) == 0:
        return mask

    start_times = np.asarray(trial_ttls, dtype=float) - float(guard_s)
    end_times = np.asarray(trial_ttls, dtype=float) + float(infusion_dur_s) + float(guard_s)

    start_idx = _time_to_frame_idx(start_times, fps=fps, n_frames=n_frames, cam1_onsets=cam1_onsets)
    end_idx = _time_to_frame_idx(end_times, fps=fps, n_frames=n_frames, cam1_onsets=cam1_onsets)

    for s, e in zip(start_idx, end_idx):
        lo, hi = (int(s), int(e)) if s <= e else (int(e), int(s))
        mask[lo:hi + 1] = True
    return mask

def get_shifted_snip_means_excluding_infusions(
    probability_vector, solenoid_ts, fps=10, pre_bins=50, post_bins=150,
    shift_frames=300, n_shuffles=1000, cam1_onsets=None, infusion_dur_s=10.0, guard_s=0.0,
 ):
    p = np.asarray(probability_vector, dtype=float).copy()
    trial_ttls = _trial_ttls_from_solenoid(solenoid_ts)
    exclusion_mask = _make_infusion_exclusion_mask(
        n_frames=len(p), trial_ttls=trial_ttls, fps=fps,
        infusion_dur_s=infusion_dur_s, guard_s=guard_s, cam1_onsets=cam1_onsets,
    )
    p[exclusion_mask] = np.nan

    shifted = p.copy()
    all_snips = []
    for _ in range(int(n_shuffles)):
        shifted = np.roll(shifted, int(shift_frames))
        snips = make_simba_snips(
            shifted, solenoid_ts, fps=fps, pre_bins=pre_bins, post_bins=post_bins, cam1_onsets=cam1_onsets
        )
        all_snips.append(snips)

    return np.concatenate(all_snips, axis=0), exclusion_mask

# Build standard null with the current n_shuffles setting (currently high-shuffle)
shuffled_snips_standard = get_shifted_snip_means(
    prob_vec, solenoid_ts, fps=fps, pre_bins=pre_bins, post_bins=post_bins,
    shift_frames=shift_frames, n_shuffles=n_shuffles, cam1_onsets=cam1_onsets,
 )

# Build exclusion null
shuffled_snips_excl, exclusion_mask = get_shifted_snip_means_excluding_infusions(
    prob_vec, solenoid_ts, fps=fps, pre_bins=pre_bins, post_bins=post_bins,
    shift_frames=shift_frames, n_shuffles=n_shuffles, cam1_onsets=cam1_onsets,
    infusion_dur_s=INFUSION_DUR_S, guard_s=EXCLUSION_GUARD_S,
 )

# Compare using the same real snips
pct_std, thr_std = get_time_above_simba_ci(
    real_snips, shuffled_snips_standard, percentile=ci_percentile, start_bin=auc_start, end_bin=auc_end
)
pct_excl, thr_excl = get_time_above_simba_ci(
    real_snips, shuffled_snips_excl, percentile=ci_percentile, start_bin=auc_start, end_bin=auc_end
)

ci_std_lo = np.nanpercentile(shuffled_snips_standard[:, auc_start:auc_end], 100 - ci_percentile, axis=0)
ci_std_hi = np.nanpercentile(shuffled_snips_standard[:, auc_start:auc_end], ci_percentile, axis=0)
ci_ex_lo = np.nanpercentile(shuffled_snips_excl[:, auc_start:auc_end], 100 - ci_percentile, axis=0)
ci_ex_hi = np.nanpercentile(shuffled_snips_excl[:, auc_start:auc_end], ci_percentile, axis=0)

summary = pd.DataFrame([
    {
        'null_model': 'standard_shifted',
        'mean_pct_time_above_ci': float(np.nanmean(pct_std)),
        'median_pct_time_above_ci': float(np.nanmedian(pct_std)),
        'mean_ci_threshold': float(thr_std),
        'median_ci_band_width': float(np.nanmedian(ci_std_hi - ci_std_lo)),
        'finite_fraction_in_null_window': float(np.mean(np.isfinite(shuffled_snips_standard[:, auc_start:auc_end]))),
    },
    {
        'null_model': f'exclude_infusions_guard_{EXCLUSION_GUARD_S:.1f}s',
        'mean_pct_time_above_ci': float(np.nanmean(pct_excl)),
        'median_pct_time_above_ci': float(np.nanmedian(pct_excl)),
        'mean_ci_threshold': float(thr_excl),
        'median_ci_band_width': float(np.nanmedian(ci_ex_hi - ci_ex_lo)),
        'finite_fraction_in_null_window': float(np.mean(np.isfinite(shuffled_snips_excl[:, auc_start:auc_end]))),
    },
])

print(f'n_shuffles used: {n_shuffles}')
print(f'Excluded source frames: {np.mean(exclusion_mask) * 100:.2f}% of trace')
display(summary.round(4))

# Optional quick visual: threshold profiles across infusion window
w_t = (np.arange(auc_start, auc_end) - pre_bins) / fps
fig, ax = plt.subplots(figsize=(7, 3.5))
ax.plot(w_t, ci_std_hi, label='Standard null threshold', color='grey', lw=1.8)
ax.plot(w_t, ci_ex_hi, label='Exclusion null threshold', color='darkorange', lw=1.8)
ax.set_xlabel('Time from infusion (s)')
ax.set_ylabel(f'{ci_percentile}th percentile threshold')
ax.set_title(f'Infusion-window thresholds: standard vs exclusion null ({STUB})')
ax.legend(fontsize=8)
fig.tight_layout()
plt.show()

# %% [markdown]
# ## Plot 2d: More data in the null (more shuffles + multi-shift pool)
#
# Expands null sampling in two ways:
#
# - More shuffles (`5000`)
# - Sampling from multiple non-overlapping shifts (instead of one fixed shift)
#
# This checks whether inference changes when the null has substantially more coverage.

# %%
# More null data: compare 1000 single-shift baseline vs 5000 single-shift vs 5000 multi-shift
BASELINE_SHUFFLES = int(n_shuffles)
MORE_SHUFFLES = 5000

def get_shifted_snip_means_multi_shift(
    probability_vector, solenoid_ts, fps=10, pre_bins=50, post_bins=150,
    shift_frames_pool=None, n_shuffles=1000, cam1_onsets=None, random_state=0,
 ):
    if shift_frames_pool is None or len(shift_frames_pool) == 0:
        raise ValueError('shift_frames_pool must contain at least one shift in frames.')

    p = np.asarray(probability_vector, dtype=float)
    rng = np.random.default_rng(random_state)
    pool = np.asarray(shift_frames_pool, dtype=int)
    all_snips = []

    for _ in range(int(n_shuffles)):
        sf = int(rng.choice(pool))
        shifted = np.roll(p, sf)
        snips = make_simba_snips(
            shifted, solenoid_ts, fps=fps, pre_bins=pre_bins, post_bins=post_bins, cam1_onsets=cam1_onsets
        )
        all_snips.append(snips)

    return np.concatenate(all_snips, axis=0)

def summarize_null_model(real_snips_local, shuffled_snips_local, label):
    pct, thr = get_time_above_simba_ci(
        real_snips_local, shuffled_snips_local, percentile=ci_percentile, start_bin=auc_start, end_bin=auc_end
    )
    ci_lo_local = np.nanpercentile(shuffled_snips_local[:, auc_start:auc_end], 100 - ci_percentile, axis=0)
    ci_hi_local = np.nanpercentile(shuffled_snips_local[:, auc_start:auc_end], ci_percentile, axis=0)
    return {
        'model': label,
        'mean_pct_time_above_ci': float(np.nanmean(pct)),
        'median_pct_time_above_ci': float(np.nanmedian(pct)),
        'mean_ci_threshold': float(thr),
        'median_ci_band_width': float(np.nanmedian(ci_hi_local - ci_lo_local)),
    }

# Baseline reference (already in memory from Plot 2c)
baseline_shuffled = shuffled_snips_standard.copy()

# Same shift, more shuffles
shuffled_more_single = get_shifted_snip_means(
    prob_vec, solenoid_ts, fps=fps, pre_bins=pre_bins, post_bins=post_bins,
    shift_frames=shift_frames, n_shuffles=MORE_SHUFFLES, cam1_onsets=cam1_onsets,
 )

# Multi-shift pool around safe shifts (seconds -> frames)
shift_pool_seconds = np.array([18, 20, 23, 25, 28, 30], dtype=float)
shift_pool_frames = np.unique(np.round(shift_pool_seconds * fps).astype(int))

shuffled_more_multi = get_shifted_snip_means_multi_shift(
    prob_vec, solenoid_ts, fps=fps, pre_bins=pre_bins, post_bins=post_bins,
    shift_frames_pool=shift_pool_frames, n_shuffles=MORE_SHUFFLES,
    cam1_onsets=cam1_onsets, random_state=42,
 )

comparison = pd.DataFrame([
    summarize_null_model(real_snips, baseline_shuffled, f'baseline_single_shift_{BASELINE_SHUFFLES}'),
    summarize_null_model(real_snips, shuffled_more_single, f'single_shift_{MORE_SHUFFLES}'),
    summarize_null_model(real_snips, shuffled_more_multi, f'multi_shift_{MORE_SHUFFLES}'),
])

display(comparison.round(4))
print(f'Shift pool (s): {shift_pool_seconds.tolist()}')
print(f'Shift pool (frames @ {fps} fps): {shift_pool_frames.tolist()}')

# Keep the larger single-shift null as default for downstream cells unless you prefer multi-shift
n_shuffles = MORE_SHUFFLES
shuffled_snips = shuffled_more_single
pct_time, mean_threshold = get_time_above_simba_ci(
    real_snips, shuffled_snips, percentile=ci_percentile, start_bin=auc_start, end_bin=auc_end
)
print(f'Updated defaults -> n_shuffles={n_shuffles}, shuffled_snips=single_shift_{MORE_SHUFFLES}')

# %% [markdown]
# ## Plot 2e: Repeat null-model tests on 3 additional deplete-10 subjects
#
# Runs the same comparison used above on three additional subjects from:
#
# - `condition == deplete`
# - `infusiontype == 10NaCl`
#
# using the current high-shuffle settings and shift pool.

# %%
# Batch test on 3 additional deplete-10NaCl subjects
candidate_ids = sorted(
    x_array.query("condition == 'deplete' and infusiontype == '10NaCl'")['id'].dropna().unique()
)
test_ids = [sid for sid in candidate_ids if sid != SUBJECT_ID][:3]

if len(test_ids) < 3:
    print(f'Only found {len(test_ids)} additional deplete-10 subjects: {test_ids}')
else:
    print(f'Testing subjects: {test_ids}')

batch_rows = []
for sid in test_ids:
    stub = stub_lookup.get(sid)
    if stub is None:
        print(f'Skipping {sid}: no stub in FileKey')
        continue

    # Load session-specific data
    pvec, _ = read_simba_probability(stub, simba_folder, probability_column=prob_col)
    s_ttls = get_ttls(stub, data_folder)
    try:
        s_cam1 = get_cam1_onsets_for_stub(stub, tank_folder)
    except Exception:
        s_cam1 = None

    s_real = make_simba_snips(
        pvec, s_ttls, fps=fps, pre_bins=pre_bins, post_bins=post_bins, cam1_onsets=s_cam1
    )

    # Standard single-shift null (same shift_frames, same n_shuffles)
    s_std = get_shifted_snip_means(
        pvec, s_ttls, fps=fps, pre_bins=pre_bins, post_bins=post_bins,
        shift_frames=shift_frames, n_shuffles=n_shuffles, cam1_onsets=s_cam1,
    )

    # Multi-shift null
    s_multi = get_shifted_snip_means_multi_shift(
        pvec, s_ttls, fps=fps, pre_bins=pre_bins, post_bins=post_bins,
        shift_frames_pool=shift_pool_frames, n_shuffles=n_shuffles, cam1_onsets=s_cam1, random_state=42,
    )

    # Exclusion null
    s_excl, s_ex_mask = get_shifted_snip_means_excluding_infusions(
        pvec, s_ttls, fps=fps, pre_bins=pre_bins, post_bins=post_bins,
        shift_frames=shift_frames, n_shuffles=n_shuffles, cam1_onsets=s_cam1,
        infusion_dur_s=INFUSION_DUR_S, guard_s=EXCLUSION_GUARD_S,
    )

    row_std = summarize_null_model(s_real, s_std, 'single_shift')
    row_multi = summarize_null_model(s_real, s_multi, 'multi_shift')
    row_excl = summarize_null_model(s_real, s_excl, f'exclusion_guard_{EXCLUSION_GUARD_S:.1f}s')

    for row in [row_std, row_multi, row_excl]:
        row['subject_id'] = sid
        row['stub'] = stub
        row['n_shuffles'] = int(n_shuffles)
        row['excluded_frame_frac'] = float(np.mean(s_ex_mask)) if 'exclusion' in row['model'] else 0.0
        batch_rows.append(row)

batch_df = pd.DataFrame(batch_rows)
if not batch_df.empty:
    cols = [
        'subject_id', 'model',
        'mean_pct_time_above_ci', 'median_pct_time_above_ci',
        'mean_ci_threshold', 'median_ci_band_width',
        'excluded_frame_frac', 'n_shuffles'
    ]
    display(batch_df[cols].sort_values(['subject_id', 'model']).round(4))

    # Compact across-subject summary by model
    model_summary = (
        batch_df.groupby('model')[['mean_pct_time_above_ci', 'mean_ci_threshold', 'median_ci_band_width']].agg(['mean', 'std'])
    )
    print('\nAcross-subject summary (3 additional deplete-10 subjects):')
    display(model_summary.round(4))

# %% [markdown]
# ## Plot 2f: Trial-by-trial pct_time for all 10 deplete-10 subjects
#
# Computes `% time above CI` by trial number for each subject under each null option, then plots:
#
# - One figure per subject
# - One across-subject mean figure

# %%
# Trial-by-trial pct_time curves for all deplete-10NaCl subjects under each null option
all_deplete10_ids = sorted(
    x_array.query("condition == 'deplete' and infusiontype == '10NaCl'")['id'].dropna().unique()
)
print(f'Running trial-by-trial plots for {len(all_deplete10_ids)} deplete-10 subjects: {all_deplete10_ids}')

trial_rows = []
for sid in all_deplete10_ids:
    stub = stub_lookup.get(sid)
    if stub is None:
        print(f'Skipping {sid}: no stub in FileKey')
        continue

    pvec, _ = read_simba_probability(stub, simba_folder, probability_column=prob_col)
    s_ttls = get_ttls(stub, data_folder)
    try:
        s_cam1 = get_cam1_onsets_for_stub(stub, tank_folder)
    except Exception:
        s_cam1 = None

    s_real = make_simba_snips(
        pvec, s_ttls, fps=fps, pre_bins=pre_bins, post_bins=post_bins, cam1_onsets=s_cam1
    )

    nulls = {
        'single_shift': get_shifted_snip_means(
            pvec, s_ttls, fps=fps, pre_bins=pre_bins, post_bins=post_bins,
            shift_frames=shift_frames, n_shuffles=n_shuffles, cam1_onsets=s_cam1,
        ),
        'multi_shift': get_shifted_snip_means_multi_shift(
            pvec, s_ttls, fps=fps, pre_bins=pre_bins, post_bins=post_bins,
            shift_frames_pool=shift_pool_frames, n_shuffles=n_shuffles, cam1_onsets=s_cam1, random_state=42,
        ),
        f'exclusion_guard_{EXCLUSION_GUARD_S:.1f}s': get_shifted_snip_means_excluding_infusions(
            pvec, s_ttls, fps=fps, pre_bins=pre_bins, post_bins=post_bins,
            shift_frames=shift_frames, n_shuffles=n_shuffles, cam1_onsets=s_cam1,
            infusion_dur_s=INFUSION_DUR_S, guard_s=EXCLUSION_GUARD_S,
        )[0],
    }

    for model_name, s_null in nulls.items():
        pct_local, _ = get_time_above_simba_ci(
            s_real, s_null, percentile=ci_percentile, start_bin=auc_start, end_bin=auc_end
        )
        for trial_idx, pct_val in enumerate(pct_local, start=1):
            trial_rows.append({
                'subject_id': sid,
                'stub': stub,
                'trial_number': int(trial_idx),
                'model': model_name,
                'pct_time_above_ci': float(pct_val),
            })

trial_df = pd.DataFrame(trial_rows)
display(trial_df.head())

model_colors = {
    'single_shift': 'grey',
    'multi_shift': 'darkorange',
    f'exclusion_guard_{EXCLUSION_GUARD_S:.1f}s': 'steelblue',
}

# One figure per subject
for sid in all_deplete10_ids:
    subj = trial_df.query('subject_id == @sid').copy()
    if subj.empty:
        continue

    fig, ax = plt.subplots(figsize=(7.5, 4))
    for model_name, grp in subj.groupby('model'):
        grp = grp.sort_values('trial_number')
        ax.plot(
            grp['trial_number'], grp['pct_time_above_ci'] * 100,
            marker='o', ms=3, lw=1.8, label=model_name,
            color=model_colors.get(model_name, None), alpha=0.9,
        )

    ax.axhline(0, color='k', lw=0.8, ls='--', alpha=0.5)
    ax.set_title(f'{sid} ({stub_lookup[sid]})')
    ax.set_xlabel('Trial number')
    ax.set_ylabel('% time above CI')
    ax.legend(fontsize=8)
    fig.tight_layout()
    plt.show()

# Across-subject mean +/- SEM by trial number
mean_df = (
    trial_df.groupby(['model', 'trial_number'])['pct_time_above_ci']
    .agg(['mean', 'std', 'count'])
    .reset_index()
)
mean_df['sem'] = mean_df['std'] / np.sqrt(mean_df['count'].clip(lower=1))

fig, ax = plt.subplots(figsize=(8, 4.5))
for model_name, grp in mean_df.groupby('model'):
    grp = grp.sort_values('trial_number')
    y = grp['mean'] * 100
    sem = grp['sem'] * 100
    x = grp['trial_number']
    color = model_colors.get(model_name, None)
    ax.plot(x, y, lw=2.2, label=model_name, color=color)
    ax.fill_between(x, y - sem, y + sem, alpha=0.2, color=color)

ax.axhline(0, color='k', lw=0.8, ls='--', alpha=0.5)
ax.set_title('All 10 deplete-10 subjects: mean % time above CI by trial')
ax.set_xlabel('Trial number')
ax.set_ylabel('% time above CI (mean +/- SEM)')
ax.legend(fontsize=8)
fig.tight_layout()
plt.show()

print(f'Used n_shuffles={n_shuffles} and shift pool {shift_pool_seconds.tolist()} s for multi-shift.')

# %% [markdown]
# ## Plot 2g: Sweep CI thresholds across groups
#
# Tests whether a less strict CI percentile gives a cleaner group-level pattern using the multi-shift null.
#
# Goal:
#
# - Deplete-10NaCl: relatively high and consistent `% time above CI`
# - Replete-10NaCl and Replete-45NaCl: low to very low `% time above CI`
#
# To keep runtime reasonable, this sweep uses a smaller shuffle count than the per-subject plotting section.

# %%
# Sweep CI percentiles and score how well they separate deplete-10 from the replete groups
THRESHOLD_SWEEP_SHUFFLES = 1000
threshold_percentiles = [80, 85, 90, 92.5, 95, 97.5]

full_file_key = pd.concat([
    pd.read_csv(data_folder / '10NaCl_FileKey.csv'),
    pd.read_csv(data_folder / '45NaCl_FileKey.csv'),
]).reset_index(drop=True)

group_sessions = {
    'deplete_10': full_file_key.query("TreatNum == 1 and `Physiological state` == 'Sodium Depleted'")
        [['Subject', 'Folder']].drop_duplicates().rename(columns={'Subject': 'subject_id', 'Folder': 'stub'}),
    'replete_10': full_file_key.query("TreatNum == 1 and `Physiological state` == 'Sodium Replete'")
        [['Subject', 'Folder']].drop_duplicates().rename(columns={'Subject': 'subject_id', 'Folder': 'stub'}),
    'replete_45': full_file_key.query("TreatNum == 45 and `Physiological state` == 'Sodium Replete'")
        [['Subject', 'Folder']].drop_duplicates().rename(columns={'Subject': 'subject_id', 'Folder': 'stub'}),
}

print('Threshold sweep session counts:')
for group_name, session_df in group_sessions.items():
    print(f'  {group_name}: {len(session_df)} sessions')

# Cache real snips and a multi-shift null per session so we can sweep thresholds cheaply.
session_cache = {}
for group_name, session_df in group_sessions.items():
    for row in session_df.itertuples(index=False):
        stub = row.stub
        if stub in session_cache:
            continue

        pvec, _ = read_simba_probability(stub, simba_folder, probability_column=prob_col)
        s_ttls = get_ttls(stub, data_folder)
        try:
            s_cam1 = get_cam1_onsets_for_stub(stub, tank_folder)
        except Exception:
            s_cam1 = None

        s_real = make_simba_snips(
            pvec, s_ttls, fps=fps, pre_bins=pre_bins, post_bins=post_bins, cam1_onsets=s_cam1
        )
        s_multi = get_shifted_snip_means_multi_shift(
            pvec, s_ttls, fps=fps, pre_bins=pre_bins, post_bins=post_bins,
            shift_frames_pool=shift_pool_frames, n_shuffles=THRESHOLD_SWEEP_SHUFFLES,
            cam1_onsets=s_cam1, random_state=42,
        )
        session_cache[stub] = {'real_snips': s_real, 'null_snips': s_multi}

sweep_trial_rows = []
for percentile in threshold_percentiles:
    for group_name, session_df in group_sessions.items():
        for row in session_df.itertuples(index=False):
            cached = session_cache[row.stub]
            pct_vals, _ = get_time_above_simba_ci(
                cached['real_snips'], cached['null_snips'],
                percentile=percentile, start_bin=auc_start, end_bin=auc_end,
            )
            for trial_idx, pct_val in enumerate(pct_vals, start=1):
                sweep_trial_rows.append({
                    'percentile': float(percentile),
                    'group': group_name,
                    'subject_id': row.subject_id,
                    'stub': row.stub,
                    'trial_number': int(trial_idx),
                    'pct_time_above_ci': float(pct_val),
                })

sweep_trial_df = pd.DataFrame(sweep_trial_rows)

summary_rows = []
for percentile in threshold_percentiles:
    subset = sweep_trial_df.query('percentile == @percentile')
    dep = subset.query("group == 'deplete_10'")['pct_time_above_ci'].to_numpy()
    rep10 = subset.query("group == 'replete_10'")['pct_time_above_ci'].to_numpy()
    rep45 = subset.query("group == 'replete_45'")['pct_time_above_ci'].to_numpy()

    dep_q25 = float(np.nanpercentile(dep, 25))
    dep_median = float(np.nanmedian(dep))
    dep_mean = float(np.nanmean(dep))
    dep_std = float(np.nanstd(dep))

    rep10_mean = float(np.nanmean(rep10))
    rep45_mean = float(np.nanmean(rep45))
    rep10_q75 = float(np.nanpercentile(rep10, 75))
    rep45_q75 = float(np.nanpercentile(rep45, 75))

    separation_score = dep_q25 - max(rep10_q75, rep45_q75)
    contrast_score = dep_mean - 0.5 * (rep10_mean + rep45_mean)

    summary_rows.append({
        'percentile': float(percentile),
        'deplete10_mean': dep_mean,
        'deplete10_median': dep_median,
        'deplete10_q25': dep_q25,
        'deplete10_std': dep_std,
        'replete10_mean': rep10_mean,
        'replete45_mean': rep45_mean,
        'replete_max_q75': max(rep10_q75, rep45_q75),
        'separation_score_q25_minus_replete_q75': separation_score,
        'contrast_score_mean_gap': contrast_score,
    })

threshold_summary_df = pd.DataFrame(summary_rows).sort_values(
    ['separation_score_q25_minus_replete_q75', 'contrast_score_mean_gap'], ascending=False
).reset_index(drop=True)

display((threshold_summary_df * pd.Series({
    'percentile': 1,
    'deplete10_mean': 100,
    'deplete10_median': 100,
    'deplete10_q25': 100,
    'deplete10_std': 100,
    'replete10_mean': 100,
    'replete45_mean': 100,
    'replete_max_q75': 100,
    'separation_score_q25_minus_replete_q75': 100,
    'contrast_score_mean_gap': 100,
})).round(2))

best_threshold = float(threshold_summary_df.iloc[0]['percentile'])
print(f'Best threshold by separation score: {best_threshold}th percentile')

# Plot 1: group means by threshold
group_mean_df = (
    sweep_trial_df.groupby(['percentile', 'group'])['pct_time_above_ci']
    .agg(['mean', 'std', 'count'])
    .reset_index()
)
group_mean_df['sem'] = group_mean_df['std'] / np.sqrt(group_mean_df['count'].clip(lower=1))

group_colors = {
    'deplete_10': 'crimson',
    'replete_10': 'steelblue',
    'replete_45': 'seagreen',
}

fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))

ax = axes[0]
for group_name, grp in group_mean_df.groupby('group'):
    grp = grp.sort_values('percentile')
    x = grp['percentile']
    y = grp['mean'] * 100
    sem = grp['sem'] * 100
    color = group_colors[group_name]
    ax.plot(x, y, marker='o', lw=2, color=color, label=group_name)
    ax.fill_between(x, y - sem, y + sem, color=color, alpha=0.18)
ax.set_xlabel('CI percentile threshold')
ax.set_ylabel('% time above CI')
ax.set_title('Group means vs threshold')
ax.legend(fontsize=8)

ax = axes[1]
plot_summary = threshold_summary_df.sort_values('percentile')
ax.plot(plot_summary['percentile'], plot_summary['separation_score_q25_minus_replete_q75'] * 100,
        marker='o', lw=2, color='black', label='Dep q25 - max rep q75')
ax.plot(plot_summary['percentile'], plot_summary['contrast_score_mean_gap'] * 100,
        marker='s', lw=2, color='darkorange', label='Dep mean - mean replete')
ax.axhline(0, color='k', lw=0.8, ls='--', alpha=0.5)
ax.set_xlabel('CI percentile threshold')
ax.set_ylabel('Separation score (percentage points)')
ax.set_title('Threshold ranking metrics')
ax.legend(fontsize=8)

fig.tight_layout()
plt.show()

# Plot 2: trial-by-trial mean curves for the best threshold
best_trial_df = sweep_trial_df.query('percentile == @best_threshold').copy()
best_mean_df = (
    best_trial_df.groupby(['group', 'trial_number'])['pct_time_above_ci']
    .agg(['mean', 'std', 'count'])
    .reset_index()
)
best_mean_df['sem'] = best_mean_df['std'] / np.sqrt(best_mean_df['count'].clip(lower=1))

fig, ax = plt.subplots(figsize=(8, 4.5))
for group_name, grp in best_mean_df.groupby('group'):
    grp = grp.sort_values('trial_number')
    x = grp['trial_number']
    y = grp['mean'] * 100
    sem = grp['sem'] * 100
    color = group_colors[group_name]
    ax.plot(x, y, lw=2.2, color=color, label=group_name)
    ax.fill_between(x, y - sem, y + sem, color=color, alpha=0.18)
ax.axhline(0, color='k', lw=0.8, ls='--', alpha=0.5)
ax.set_xlabel('Trial number')
ax.set_ylabel('% time above CI (mean +/- SEM)')
ax.set_title(f'Best threshold = {best_threshold}th percentile')
ax.legend(fontsize=8)
fig.tight_layout()
plt.show()

print(f'Used multi-shift null with {THRESHOLD_SWEEP_SHUFFLES} shuffles per session for threshold sweep.')

# %% [markdown]
# ## Plot 3: Histogram — shuffled vs real bin values in infusion window
#
# Shows whether real trial probability values (in the infusion window) are separated
# from the null distribution of shuffled bin values.

# %%
# Flatten infusion-window bins from shuffled and real snips
shuf_window_vals = shuffled_snips[:, auc_start:auc_end].flatten()
shuf_window_vals = shuf_window_vals[~np.isnan(shuf_window_vals)]

real_window_vals = real_snips[:, auc_start:auc_end].flatten()
real_window_vals = real_window_vals[~np.isnan(real_window_vals)]

threshold_scalar = float(np.nanpercentile(shuf_window_vals, ci_percentile))

fig, ax = plt.subplots(figsize=(7, 4))

bins = np.linspace(
    min(shuf_window_vals.min(), real_window_vals.min()),
    max(shuf_window_vals.max(), real_window_vals.max()),
    60
)

ax.hist(shuf_window_vals, bins=bins, density=True, alpha=0.5,
        color='grey', label=f'Shuffled (n={len(shuf_window_vals):,} bins)')
ax.hist(real_window_vals, bins=bins, density=True, alpha=0.6,
        color='steelblue', label=f'Real trials (n={len(real_window_vals):,} bins)')

ax.axvline(threshold_scalar, color='red', lw=1.5, ls='--',
           label=f'{ci_percentile}th pct of shuffled = {threshold_scalar:.3f}')

pct_real_above = np.mean(real_window_vals > threshold_scalar) * 100
ax.set_title(f'{STUB} — infusion window bins\n'
             f'{pct_real_above:.1f}% of real bins above threshold')
ax.set_xlabel('Appetitive probability')
ax.set_ylabel('Density')
ax.legend(fontsize=9)
fig.tight_layout()
plt.show()

# %% [markdown]
# ## Plot 4: Per-trial pct_time across trials (by trial number)
#
# Checks whether there is a trial-order trend in the pct_time metric.

# %%
# Compare trial-number trend from assembled x_array for this subject
subject_rows = x_array.query('id == @SUBJECT_ID').copy()

fig, axes = plt.subplots(1, 2, figsize=(11, 4))

# Left: pct_time by trial number from assembled data
for (cond, inf), grp in subject_rows.groupby(['condition', 'infusiontype']):
    axes[0].scatter(grp['trial'], grp['simba_pct_time_above_95ci'],
                    label=f'{cond} {inf}', alpha=0.7, s=30)
axes[0].axhline(0, color='k', lw=0.8, ls='--')
axes[0].set_xlabel('Trial number')
axes[0].set_ylabel('% time above CI')
axes[0].set_title(f'{SUBJECT_ID} ({STUB}) — pct_time by trial (from pickle)')
axes[0].legend(fontsize=8)

# Right: pct_time from the freshly recomputed shuffle
axes[1].scatter(np.arange(len(pct_time)), pct_time * 100,
                color='steelblue', alpha=0.8, s=30)
axes[1].axhline(0, color='k', lw=0.8, ls='--')
axes[1].set_xlabel('Trial number')
axes[1].set_ylabel('% time above CI')
axes[1].set_title(f'{SUBJECT_ID} ({STUB}) — pct_time from fresh recompute')

fig.tight_layout()
plt.show()


# %% [markdown]
# ## Plot 5: Mean snips ± CI band, all four groups
#
# Uses the already-assembled `snips_simba` and overlays the CI band computed from a
# pooled shuffle across all infusion-window bins.

# %%
from figure_plotting import smooth_array

snips_smooth = np.array(smooth_array(snips_simba))

group_order = [
    ('replete', '10NaCl'),
    ('replete', '45NaCl'),
    ('deplete', '10NaCl'),
    ('deplete', '45NaCl'),
]
group_labels = ['Replete 10NaCl', 'Replete 45NaCl', 'Deplete 10NaCl', 'Deplete 45NaCl']

fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True)

fps_plot = params.get('simba_fps', 10)
pre_plot = params.get('simba_pre_bins', 50)
post_plot = params.get('simba_post_bins', 150)
n_bins_plot = pre_plot + post_plot
time_ax = np.arange(n_bins_plot) / fps_plot - pre_plot / fps_plot

for ax, (cond, inf), label, col in zip(axes.flat, group_order, group_labels, colors):
    mask = ((x_array['condition'] == cond) & (x_array['infusiontype'] == inf)).values
    grp_snips = snips_smooth[mask]

    # Animal-averaged snips
    animal_ids = x_array.loc[mask, 'id'].values
    unique_ids = np.unique(animal_ids)
    animal_means = np.array([
        np.nanmean(grp_snips[animal_ids == aid], axis=0)
        for aid in unique_ids
    ])

    grand_mean = np.nanmean(animal_means, axis=0)
    sem_curve  = np.nanstd(animal_means, axis=0) / np.sqrt(len(unique_ids))

    # Individual trial snips (translucent)
    for snip in grp_snips:
        ax.plot(time_ax, snip, color=col, alpha=0.07, lw=0.6)

    # Mean ± SEM
    ax.fill_between(time_ax, grand_mean - sem_curve, grand_mean + sem_curve,
                    color=col, alpha=0.3)
    ax.plot(time_ax, grand_mean, color=col, lw=2)

    ax.axvspan(0, 10, color='k', alpha=0.05, zorder=-10)
    ax.axvline(0, color='k', lw=0.8, ls='--')
    ax.axhline(0, color='k', lw=0.5, ls=':')
    ax.set_title(f'{label} (n={len(unique_ids)} animals, {mask.sum()} trials)')
    ax.set_xlabel('Time from infusion (s)')
    ax.set_ylabel('Appetitive prob. (baseline-sub)')

fig.suptitle('All groups: mean ± SEM snips (assembled data)', fontsize=12)
fig.tight_layout()
plt.show()

# %%
# Debug mapping for replete/deplete groups
full_file_key = pd.concat([
    pd.read_csv(data_folder / '10NaCl_FileKey.csv'),
    pd.read_csv(data_folder / '45NaCl_FileKey.csv'),
]).reset_index(drop=True)
display(
    full_file_key.groupby(['TreatNum', 'Physiological state']).size().reset_index(name='n_rows')
        .sort_values(['TreatNum', 'Physiological state']).reset_index(drop=True)
 )
print('x_array group counts by condition/infusiontype:')
display(x_array.groupby(['condition', 'infusiontype'])['id'].nunique().reset_index(name='n_subjects'))

# %% [markdown]
# ## Plot 2h — threshold sweep below 80th percentile\n\nExtend the sweep to lower percentiles to see if there is a better separation point below 80.

# %%
# Extend sweep below 80 using already-cached snips (no recomputation needed)
low_percentiles = [60, 65, 70, 75, 80]

low_rows = []
for percentile in low_percentiles:
    for group_name, session_df in group_sessions.items():
        for row in session_df.itertuples(index=False):
            cached = session_cache[row.stub]
            pct_vals, _ = get_time_above_simba_ci(
                cached['real_snips'], cached['null_snips'],
                percentile=percentile, start_bin=auc_start, end_bin=auc_end,
            )
            for trial_idx, pct_val in enumerate(pct_vals, start=1):
                low_rows.append({
                    'percentile': float(percentile),
                    'group': group_name,
                    'subject_id': row.subject_id,
                    'stub': row.stub,
                    'trial_number': int(trial_idx),
                    'pct_time_above_ci': float(pct_val),
                })

low_trial_df = pd.DataFrame(low_rows)

# Combine with existing sweep for unified plotting
combined_df = pd.concat([low_trial_df, sweep_trial_df], ignore_index=True)
all_percentiles = sorted(combined_df['percentile'].unique())

low_summary_rows = []
for percentile in all_percentiles:
    subset = combined_df.query('percentile == @percentile')
    dep = subset.query("group == 'deplete_10'")['pct_time_above_ci'].to_numpy()
    rep10 = subset.query("group == 'replete_10'")['pct_time_above_ci'].to_numpy()
    rep45 = subset.query("group == 'replete_45'")['pct_time_above_ci'].to_numpy()

    dep_q25 = float(np.nanpercentile(dep, 25))
    dep_mean = float(np.nanmean(dep))
    rep10_mean = float(np.nanmean(rep10))
    rep45_mean = float(np.nanmean(rep45))
    rep10_q75 = float(np.nanpercentile(rep10, 75))
    rep45_q75 = float(np.nanpercentile(rep45, 75))
    separation_score = dep_q25 - max(rep10_q75, rep45_q75)
    contrast_score = dep_mean - 0.5 * (rep10_mean + rep45_mean)

    low_summary_rows.append({
        'percentile': float(percentile),
        'deplete10_mean': dep_mean,
        'deplete10_q25': dep_q25,
        'replete10_mean': rep10_mean,
        'replete45_mean': rep45_mean,
        'separation_score': separation_score,
        'contrast_score': contrast_score,
    })

combined_summary = pd.DataFrame(low_summary_rows).sort_values('percentile').reset_index(drop=True)
print(combined_summary[['percentile', 'deplete10_mean', 'deplete10_q25',
                          'replete10_mean', 'replete45_mean',
                          'separation_score', 'contrast_score']].to_string(index=False))

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

ax = axes[0]
ax.plot(combined_summary['percentile'], combined_summary['deplete10_mean'] * 100,
        'o-', color='#d62728', label='deplete-10 (mean)')
ax.fill_between(
    combined_summary['percentile'],
    (combined_summary['deplete10_mean'] - combined_summary['deplete10_mean'].std()) * 100,
    (combined_summary['deplete10_mean'] + combined_summary['deplete10_mean'].std()) * 100,
    alpha=0.15, color='#d62728',
)
ax.plot(combined_summary['percentile'], combined_summary['replete10_mean'] * 100,
        's--', color='#1f77b4', label='replete-10 (mean)')
ax.plot(combined_summary['percentile'], combined_summary['replete45_mean'] * 100,
        '^--', color='#2ca02c', label='replete-45 (mean)')
ax.axvline(80, color='gray', lw=0.8, linestyle=':', label='80th')
ax.set_xlabel('CI percentile threshold')
ax.set_ylabel('% time above CI (%)')
ax.set_title('Group means vs threshold (full range)')
ax.legend(fontsize=8)

ax = axes[1]
ax.plot(combined_summary['percentile'], combined_summary['separation_score'] * 100,
        'ko-', label='separation (dep Q25 − rep Q75)')
ax.plot(combined_summary['percentile'], combined_summary['contrast_score'] * 100,
        'rs--', label='contrast (dep mean − avg rep mean)')
ax.axhline(0, color='gray', lw=0.8, linestyle=':')
ax.axvline(80, color='gray', lw=0.8, linestyle=':', label='80th')
best_row = combined_summary.loc[combined_summary['separation_score'].idxmax()]
ax.axvline(best_row['percentile'], color='orange', lw=1.5, linestyle='--',
           label=f"best sep: {best_row['percentile']:.0f}th")
ax.set_xlabel('CI percentile threshold')
ax.set_ylabel('Score (pp)')
ax.set_title('Separation & contrast scores')
ax.legend(fontsize=8)

fig.tight_layout()
plt.show()
print(f"\nBest separation at: {best_row['percentile']:.1f}th percentile "
      f"(score={best_row['separation_score']*100:+.1f} pp)")


# %% [markdown]
# ## Plot 2i — net score: pct_above − pct_below\n\nFor each upper percentile threshold *p*, also compute the fraction of bins **below** the symmetric lower bound (100 − p). The net score `pct_above − pct_below` is centred around zero: positive = more time above than below the null, negative = suppressed.\n\nExample: threshold = 80th → upper bound at p80, lower bound at p20.

# %%
def get_time_below_simba_ci(real_snips, shifted_snips, percentile=80, start_bin=50, end_bin=150):
    """Per-trial fraction of bins BELOW the lower null-CI bound.

    The lower bound is the (100 - percentile)th percentile of the null
    distribution, symmetric with the upper bound used by get_time_above_simba_ci.
    E.g. percentile=80 → lower bound at 20th percentile of null.
    """
    lower_pct = 100.0 - percentile
    threshold_per_bin = np.nanpercentile(
        np.asarray(shifted_snips, dtype=float)[:, start_bin:end_bin],
        lower_pct,
        axis=0,
    )
    proportions = []
    for snip in np.asarray(real_snips, dtype=float):
        window = np.asarray(snip[start_bin:end_bin], dtype=float)
        valid_mask = ~np.isnan(window)
        if not np.any(valid_mask):
            proportions.append(np.nan)
        else:
            proportions.append(float(np.mean(window[valid_mask] < threshold_per_bin[valid_mask])))
    return np.asarray(proportions, dtype=float)


# ── Re-sweep with both above and below metrics ────────────────────────────────
sweep_percentiles = sorted(combined_df['percentile'].unique())   # [60,65,...,97.5]

net_rows = []
for percentile in sweep_percentiles:
    for group_name, session_df in group_sessions.items():
        for row in session_df.itertuples(index=False):
            cached = session_cache[row.stub]
            pct_above, _ = get_time_above_simba_ci(
                cached['real_snips'], cached['null_snips'],
                percentile=percentile, start_bin=auc_start, end_bin=auc_end,
            )
            pct_below = get_time_below_simba_ci(
                cached['real_snips'], cached['null_snips'],
                percentile=percentile, start_bin=auc_start, end_bin=auc_end,
            )
            for trial_idx, (pa, pb) in enumerate(zip(pct_above, pct_below), start=1):
                net_rows.append({
                    'percentile': float(percentile),
                    'group': group_name,
                    'subject_id': row.subject_id,
                    'stub': row.stub,
                    'trial_number': int(trial_idx),
                    'pct_above': float(pa),
                    'pct_below': float(pb),
                    'net_score': float(pa) - float(pb),
                })

net_df = pd.DataFrame(net_rows)

# ── Summary per percentile ────────────────────────────────────────────────────
net_summary_rows = []
for percentile in sweep_percentiles:
    sub = net_df.query('percentile == @percentile')
    for grp in ['deplete_10', 'replete_10', 'replete_45']:
        vals = sub.query('group == @grp')['net_score'].to_numpy()
        net_summary_rows.append({
            'percentile': percentile,
            'group': grp,
            'mean_net': float(np.nanmean(vals)),
            'sem_net': float(np.nanstd(vals) / np.sqrt(np.sum(~np.isnan(vals)))),
            'q25_net': float(np.nanpercentile(vals, 25)),
            'q75_net': float(np.nanpercentile(vals, 75)),
        })

net_summary = pd.DataFrame(net_summary_rows)

# Print table
pivot = net_summary.pivot(index='percentile', columns='group', values='mean_net') * 100
pivot['dep10_minus_rep_avg'] = pivot['deplete_10'] - 0.5 * (pivot['replete_10'] + pivot['replete_45'])
print(pivot[['deplete_10', 'replete_10', 'replete_45', 'dep10_minus_rep_avg']].to_string(float_format='%.1f'))

# ── Plot ──────────────────────────────────────────────────────────────────────
grp_colors_net = {'deplete_10': '#d62728', 'replete_10': '#1f77b4', 'replete_45': '#2ca02c'}
grp_markers = {'deplete_10': 'o', 'replete_10': 's', 'replete_45': '^'}
grp_ls = {'deplete_10': '-', 'replete_10': '--', 'replete_45': '--'}

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Left: mean net score per group across percentiles
ax = axes[0]
for grp in ['deplete_10', 'replete_10', 'replete_45']:
    sub = net_summary.query('group == @grp').sort_values('percentile')
    ax.plot(sub['percentile'], sub['mean_net'] * 100,
            marker=grp_markers[grp], linestyle=grp_ls[grp],
            color=grp_colors_net[grp], label=grp.replace('_', '-'))
    ax.fill_between(sub['percentile'],
                    (sub['mean_net'] - sub['sem_net']) * 100,
                    (sub['mean_net'] + sub['sem_net']) * 100,
                    alpha=0.15, color=grp_colors_net[grp])
ax.axhline(0, color='gray', lw=0.8, linestyle=':')
ax.set_xlabel('Upper CI percentile (lower bound = 100 − p)')
ax.set_ylabel('Net score: pct_above − pct_below (%)')
ax.set_title('Mean net score per group')
ax.legend(fontsize=9)

# Right: bar chart at best-separating percentile (use contrast on net score)
contrast = net_summary.groupby('percentile').apply(
    lambda d: d.set_index('group')['mean_net'].get('deplete_10', np.nan)
              - 0.5 * (d.set_index('group')['mean_net'].get('replete_10', np.nan)
                       + d.set_index('group')['mean_net'].get('replete_45', np.nan))
).rename('contrast')
best_pct = float(contrast.idxmax())

ax = axes[1]
sub_best = net_df.query('percentile == @best_pct')
for gi, grp in enumerate(['deplete_10', 'replete_10', 'replete_45']):
    vals = sub_best.query('group == @grp')['net_score'].to_numpy() * 100
    mean_v = np.nanmean(vals)
    sem_v = np.nanstd(vals) / np.sqrt(np.sum(~np.isnan(vals)))
    ax.bar(gi, mean_v, yerr=sem_v, color=grp_colors_net[grp],
           alpha=0.8, capsize=4, label=grp.replace('_', '-'))
    ax.scatter(np.full(len(vals), gi) + np.random.default_rng(0).uniform(-0.2, 0.2, len(vals)),
               vals, color=grp_colors_net[grp], s=12, alpha=0.4, zorder=3)
ax.axhline(0, color='gray', lw=0.8, linestyle=':')
ax.set_xticks([0, 1, 2])
ax.set_xticklabels(['deplete-10', 'replete-10', 'replete-45'], fontsize=9)
ax.set_ylabel('Net score (%)')
ax.set_title(f'Net score at best percentile ({best_pct:.0f}th)')

fig.tight_layout()
plt.show()
print(f"\nBest contrast at: {best_pct:.1f}th percentile")


# %% [markdown]
# ## Plot 2j — median-based net score\n\nUse the null **median** (50th percentile) as the threshold. Because the lower bound is also the 50th percentile, the net score simplifies to:\n\n$$\\text{net} = 2 \\times \\text{pct\\_above\\_median} - 1$$\n\nRanges from −1 (always below null) to +1 (always above null), centred at 0 when signal matches the null.

# %%
# Compute pct_above_median for each session/trial, then derive net = 2*pct_above - 1
median_rows = []
for group_name, session_df in group_sessions.items():
    for row in session_df.itertuples(index=False):
        cached = session_cache[row.stub]
        pct_above_med, _ = get_time_above_simba_ci(
            cached['real_snips'], cached['null_snips'],
            percentile=50, start_bin=auc_start, end_bin=auc_end,
        )
        for trial_idx, pa in enumerate(pct_above_med, start=1):
            median_rows.append({
                'group': group_name,
                'subject_id': row.subject_id,
                'stub': row.stub,
                'trial_number': int(trial_idx),
                'pct_above_median': float(pa),
                'net_score': 2.0 * float(pa) - 1.0,
            })

median_df = pd.DataFrame(median_rows)

# ── Summary by group ──────────────────────────────────────────────────────────
grp_summary = (
    median_df.groupby('group')['net_score']
    .agg(mean='mean', std='std', count='count')
    .assign(sem=lambda d: d['std'] / np.sqrt(d['count']))
    .loc[['deplete_10', 'replete_10', 'replete_45']]
)
grp_summary['mean_pct'] = grp_summary['mean'] * 100
grp_summary['sem_pct'] = grp_summary['sem'] * 100
print("Median-based net score (net = 2×pct_above_median − 1):")
print(grp_summary[['mean_pct', 'sem_pct']].rename(
    columns={'mean_pct': 'mean (%)', 'sem_pct': 'SEM (%)'}
).to_string(float_format='%.1f'))

# ── Per-subject means ─────────────────────────────────────────────────────────
subj_means = (
    median_df.groupby(['group', 'subject_id'])['net_score'].mean().reset_index()
)

# ── Plots ─────────────────────────────────────────────────────────────────────
grp_order = ['deplete_10', 'replete_10', 'replete_45']
grp_colors_j = {'deplete_10': '#d62728', 'replete_10': '#1f77b4', 'replete_45': '#2ca02c'}
rng = np.random.default_rng(1)

fig, axes = plt.subplots(1, 2, figsize=(13, 5))

# Left: bar chart with individual trial points
ax = axes[0]
for gi, grp in enumerate(grp_order):
    vals = median_df.query('group == @grp')['net_score'].to_numpy() * 100
    mean_v = np.nanmean(vals)
    sem_v = np.nanstd(vals) / np.sqrt(np.sum(~np.isnan(vals)))
    ax.bar(gi, mean_v, yerr=sem_v, color=grp_colors_j[grp], alpha=0.75, capsize=5, width=0.6)
    jitter = rng.uniform(-0.22, 0.22, len(vals))
    ax.scatter(gi + jitter, vals, color=grp_colors_j[grp], s=8, alpha=0.35, zorder=3)
ax.axhline(0, color='gray', lw=0.9, linestyle='--')
ax.set_xticks(range(len(grp_order)))
ax.set_xticklabels(['deplete-10', 'replete-10', 'replete-45'], fontsize=10)
ax.set_ylabel('Net score = 2×pct_above_median − 1 (%)', fontsize=9)
ax.set_title('Median-based net score (all trials)')

# Right: per-subject means as strip plot
ax = axes[1]
for gi, grp in enumerate(grp_order):
    vals = subj_means.query('group == @grp')['net_score'].to_numpy() * 100
    mean_v = np.nanmean(vals)
    sem_v = np.nanstd(vals) / np.sqrt(len(vals))
    ax.bar(gi, mean_v, yerr=sem_v, color=grp_colors_j[grp], alpha=0.75, capsize=5, width=0.6)
    jitter = rng.uniform(-0.18, 0.18, len(vals))
    ax.scatter(gi + jitter, vals, color=grp_colors_j[grp], s=40, zorder=3)
ax.axhline(0, color='gray', lw=0.9, linestyle='--')
ax.set_xticks(range(len(grp_order)))
ax.set_xticklabels(['deplete-10', 'replete-10', 'replete-45'], fontsize=10)
ax.set_ylabel('Mean net score per subject (%)', fontsize=9)
ax.set_title('Median-based net score (per-subject means)')

fig.tight_layout()
plt.show()


# %%
# Trial-by-trial median net score for every animal, arranged by group
grp_order_j = ['deplete_10', 'replete_10', 'replete_45']
grp_label_map = {'deplete_10': 'Deplete-10', 'replete_10': 'Replete-10', 'replete_45': 'Replete-45'}
grp_colors_k = {'deplete_10': '#d62728', 'replete_10': '#1f77b4', 'replete_45': '#2ca02c'}
grp_bg = {'deplete_10': '#fdecea', 'replete_10': '#eaf1fb', 'replete_45': '#eafaea'}

# Collect ordered list of (group, subject_id) pairs
ordered_subjects = []
for grp in grp_order_j:
    sids = sorted(median_df.query('group == @grp')['subject_id'].unique())
    for sid in sids:
        ordered_subjects.append((grp, sid))

n_subj = len(ordered_subjects)
n_cols = 5
n_rows = int(np.ceil(n_subj / n_cols))

fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 3.2, n_rows * 2.8), sharey=True)
axes_flat = axes.flatten()

for idx, (grp, sid) in enumerate(ordered_subjects):
    ax = axes_flat[idx]
    sub = median_df.query('group == @grp and subject_id == @sid').sort_values('trial_number')
    trials = sub['trial_number'].to_numpy()
    scores = sub['net_score'].to_numpy() * 100

    color = grp_colors_k[grp]
    ax.set_facecolor(grp_bg[grp])
    ax.axhline(0, color='gray', lw=0.8, linestyle='--', zorder=1)
    ax.bar(trials, scores, color=color, alpha=0.7, width=0.8, zorder=2)
    ax.plot(trials, scores, 'o', color=color, ms=3, zorder=3)

    # Running mean overlay
    if len(scores) >= 3:
        kernel = np.ones(3) / 3
        rm = np.convolve(scores, kernel, mode='same')
        ax.plot(trials, rm, color='k', lw=1.2, zorder=4, alpha=0.6)

    grp_label = grp_label_map[grp]
    ax.set_title(f'{sid}', fontsize=9, fontweight='bold', pad=2)
    ax.text(0.03, 0.97, grp_label, transform=ax.transAxes,
            fontsize=7, va='top', ha='left', color=color, fontweight='bold')
    mean_val = np.nanmean(scores)
    ax.text(0.97, 0.97, f'μ={mean_val:.0f}%', transform=ax.transAxes,
            fontsize=7, va='top', ha='right', color='k')
    ax.set_xlim(0.5, max(trials) + 0.5)
    ax.tick_params(labelsize=7)
    if idx % n_cols == 0:
        ax.set_ylabel('Net score (%)', fontsize=7)
    if idx >= (n_rows - 1) * n_cols:
        ax.set_xlabel('Trial', fontsize=7)

# Hide unused axes
for idx in range(n_subj, len(axes_flat)):
    axes_flat[idx].set_visible(False)

fig.suptitle('Trial-by-trial median net score (2×pct_above_median − 1) — all animals',
             fontsize=11, y=1.01)
fig.tight_layout()
plt.show()


# %%
# ── Load deplete-45 sessions and compute median net scores ────────────────────
dep45_sessions = (
    full_file_key
    .query("TreatNum == 45 and `Physiological state` == 'Sodium Depleted'")
    [['Subject', 'Folder']].drop_duplicates()
    .rename(columns={'Subject': 'subject_id', 'Folder': 'stub'})
)
print(f"deplete_45: {len(dep45_sessions)} sessions")
print(dep45_sessions.to_string(index=False))

new_rows = []
for row in dep45_sessions.itertuples(index=False):
    if row.stub not in session_cache:
        pvec, _ = read_simba_probability(row.stub, simba_folder, probability_column=prob_col)
        s_ttls = get_ttls(row.stub, data_folder)
        try:
            s_cam1 = get_cam1_onsets_for_stub(row.stub, tank_folder)
        except Exception:
            s_cam1 = None
        s_real = make_simba_snips(
            pvec, s_ttls, fps=fps, pre_bins=pre_bins, post_bins=post_bins, cam1_onsets=s_cam1
        )
        s_multi = get_shifted_snip_means_multi_shift(
            pvec, s_ttls, fps=fps, pre_bins=pre_bins, post_bins=post_bins,
            shift_frames_pool=shift_pool_frames, n_shuffles=THRESHOLD_SWEEP_SHUFFLES,
            cam1_onsets=s_cam1, random_state=42,
        )
        session_cache[row.stub] = {'real_snips': s_real, 'null_snips': s_multi}

    cached = session_cache[row.stub]
    pct_above_med, _ = get_time_above_simba_ci(
        cached['real_snips'], cached['null_snips'],
        percentile=50, start_bin=auc_start, end_bin=auc_end,
    )
    for trial_idx, pa in enumerate(pct_above_med, start=1):
        new_rows.append({
            'group': 'deplete_45',
            'subject_id': row.subject_id,
            'stub': row.stub,
            'trial_number': int(trial_idx),
            'pct_above_median': float(pa),
            'net_score': 2.0 * float(pa) - 1.0,
        })

# Append to median_df (avoid duplicates if re-run)
median_df = pd.concat(
    [median_df[median_df['group'] != 'deplete_45'], pd.DataFrame(new_rows)],
    ignore_index=True,
)

grp_summary_all = (
    median_df.groupby('group')['net_score']
    .agg(mean='mean', std='std', count='count')
    .assign(sem=lambda d: d['std'] / np.sqrt(d['count']))
)
grp_summary_all['mean_pct'] = grp_summary_all['mean'] * 100
grp_summary_all['sem_pct'] = grp_summary_all['sem'] * 100
print("\nMedian net score — all groups:")
print(grp_summary_all[['mean_pct', 'sem_pct']].rename(
    columns={'mean_pct': 'mean (%)', 'sem_pct': 'SEM (%)'}
).to_string(float_format='%.1f'))

# ── Trial-by-trial plot — all four groups ─────────────────────────────────────
grp_order_4 = ['deplete_10', 'deplete_45', 'replete_10', 'replete_45']
grp_label_map_4 = {
    'deplete_10': 'Deplete-10', 'deplete_45': 'Deplete-45',
    'replete_10': 'Replete-10', 'replete_45': 'Replete-45',
}
grp_colors_4 = {
    'deplete_10': '#d62728', 'deplete_45': '#ff7f0e',
    'replete_10': '#1f77b4', 'replete_45': '#2ca02c',
}
grp_bg_4 = {
    'deplete_10': '#fdecea', 'deplete_45': '#fff4e5',
    'replete_10': '#eaf1fb', 'replete_45': '#eafaea',
}

ordered_subjects_4 = []
for grp in grp_order_4:
    sids = sorted(median_df.query('group == @grp')['subject_id'].unique())
    for sid in sids:
        ordered_subjects_4.append((grp, sid))

n_subj_4 = len(ordered_subjects_4)
n_cols_4 = 5
n_rows_4 = int(np.ceil(n_subj_4 / n_cols_4))

fig, axes = plt.subplots(n_rows_4, n_cols_4, figsize=(n_cols_4 * 3.2, n_rows_4 * 2.8), sharey=True)
axes_flat = axes.flatten()
rng2 = np.random.default_rng(2)

for idx, (grp, sid) in enumerate(ordered_subjects_4):
    ax = axes_flat[idx]
    sub = median_df.query('group == @grp and subject_id == @sid').sort_values('trial_number')
    trials = sub['trial_number'].to_numpy()
    scores = sub['net_score'].to_numpy() * 100
    color = grp_colors_4[grp]

    ax.set_facecolor(grp_bg_4[grp])
    ax.axhline(0, color='gray', lw=0.8, linestyle='--', zorder=1)
    ax.bar(trials, scores, color=color, alpha=0.7, width=0.8, zorder=2)
    ax.plot(trials, scores, 'o', color=color, ms=3, zorder=3)

    if len(scores) >= 3:
        kernel = np.ones(3) / 3
        ax.plot(trials, np.convolve(scores, kernel, mode='same'),
                color='k', lw=1.2, zorder=4, alpha=0.6)

    ax.set_title(sid, fontsize=9, fontweight='bold', pad=2)
    ax.text(0.03, 0.97, grp_label_map_4[grp], transform=ax.transAxes,
            fontsize=7, va='top', ha='left', color=color, fontweight='bold')
    ax.text(0.97, 0.97, f'μ={np.nanmean(scores):.0f}%', transform=ax.transAxes,
            fontsize=7, va='top', ha='right', color='k')
    ax.set_xlim(0.5, max(trials) + 0.5)
    ax.tick_params(labelsize=7)
    if idx % n_cols_4 == 0:
        ax.set_ylabel('Net score (%)', fontsize=7)
    if idx >= (n_rows_4 - 1) * n_cols_4:
        ax.set_xlabel('Trial', fontsize=7)

for idx in range(n_subj_4, len(axes_flat)):
    axes_flat[idx].set_visible(False)

fig.suptitle('Trial-by-trial median net score — all animals, all groups',
             fontsize=11, y=1.01)
fig.tight_layout()
plt.show()


# %% [markdown]
# ## Plot 2k — mean z-score in infusion window (all groups)\n\nZ-score each session's SIMBA probability vector using the session mean and SD of the entire recording (matching `assemble_all_data.py`). For each trial, take the mean z-score across the infusion window (`auc_start:auc_end`). Zero = typical session activity; positive = elevated above session baseline.

# %%
# Build full list of sessions across all four groups
all_group_sessions = {
    **group_sessions,                       # deplete_10, replete_10, replete_45
    'deplete_45': dep45_sessions,           # added in previous cell
}

zscore_rows = []
for grp, session_df in all_group_sessions.items():
    for row in session_df.itertuples(index=False):
        # Load full probability vector to get session-level mean/std
        pvec, _ = read_simba_probability(row.stub, simba_folder, probability_column=prob_col)
        prob_mean = float(np.mean(pvec))
        prob_std = float(np.std(pvec))

        # Z-score each real snip using session-level stats
        real = np.asarray(session_cache[row.stub]['real_snips'], dtype=float)
        z_snips = (real - prob_mean) / prob_std   # shape: (n_trials, n_bins)

        # Mean z-score across the infusion window per trial
        for trial_idx, z_snip in enumerate(z_snips, start=1):
            window = z_snip[auc_start:auc_end]
            mean_z = float(np.nanmean(window))
            zscore_rows.append({
                'group': grp,
                'subject_id': row.subject_id,
                'stub': row.stub,
                'trial_number': int(trial_idx),
                'mean_z': mean_z,
            })

zscore_df = pd.DataFrame(zscore_rows)

# ── Summary by group ──────────────────────────────────────────────────────────
zgrp_order = ['deplete_10', 'deplete_45', 'replete_10', 'replete_45']
zgrp_summary = (
    zscore_df.groupby('group')['mean_z']
    .agg(mean='mean', std='std', count='count')
    .assign(sem=lambda d: d['std'] / np.sqrt(d['count']))
    .loc[zgrp_order]
)
print("Mean z-score in infusion window (all groups):")
print(zgrp_summary[['mean', 'sem']].rename(columns={'mean': 'mean z', 'sem': 'SEM'})
      .to_string(float_format='%.3f'))

# ── Trial-by-trial plot ───────────────────────────────────────────────────────
zgrp_colors = {
    'deplete_10': '#d62728', 'deplete_45': '#ff7f0e',
    'replete_10': '#1f77b4', 'replete_45': '#2ca02c',
}
zgrp_bg = {
    'deplete_10': '#fdecea', 'deplete_45': '#fff4e5',
    'replete_10': '#eaf1fb', 'replete_45': '#eafaea',
}
zgrp_label = {
    'deplete_10': 'Deplete-10', 'deplete_45': 'Deplete-45',
    'replete_10': 'Replete-10', 'replete_45': 'Replete-45',
}

ordered_z = []
for grp in zgrp_order:
    for sid in sorted(zscore_df.query('group == @grp')['subject_id'].unique()):
        ordered_z.append((grp, sid))

n_z = len(ordered_z)
n_cols_z = 5
n_rows_z = int(np.ceil(n_z / n_cols_z))

fig, axes = plt.subplots(n_rows_z, n_cols_z,
                         figsize=(n_cols_z * 3.2, n_rows_z * 2.8), sharey=True)
axes_flat = axes.flatten()
rng3 = np.random.default_rng(3)

for idx, (grp, sid) in enumerate(ordered_z):
    ax = axes_flat[idx]
    sub = zscore_df.query('group == @grp and subject_id == @sid').sort_values('trial_number')
    trials = sub['trial_number'].to_numpy()
    zvals = sub['mean_z'].to_numpy()
    color = zgrp_colors[grp]

    ax.set_facecolor(zgrp_bg[grp])
    ax.axhline(0, color='gray', lw=0.8, linestyle='--', zorder=1)
    ax.bar(trials, zvals, color=color, alpha=0.7, width=0.8, zorder=2)
    ax.plot(trials, zvals, 'o', color=color, ms=3, zorder=3)

    if len(zvals) >= 3:
        ax.plot(trials, np.convolve(zvals, np.ones(3) / 3, mode='same'),
                color='k', lw=1.2, zorder=4, alpha=0.6)

    ax.set_title(sid, fontsize=9, fontweight='bold', pad=2)
    ax.text(0.03, 0.97, zgrp_label[grp], transform=ax.transAxes,
            fontsize=7, va='top', ha='left', color=color, fontweight='bold')
    ax.text(0.97, 0.97, f'μ={np.nanmean(zvals):.2f}z', transform=ax.transAxes,
            fontsize=7, va='top', ha='right', color='k')
    ax.set_xlim(0.5, max(trials) + 0.5)
    ax.tick_params(labelsize=7)
    if idx % n_cols_z == 0:
        ax.set_ylabel('Mean z-score', fontsize=7)
    if idx >= (n_rows_z - 1) * n_cols_z:
        ax.set_xlabel('Trial', fontsize=7)

for idx in range(n_z, len(axes_flat)):
    axes_flat[idx].set_visible(False)

fig.suptitle('Mean z-score in infusion window — all animals, all groups',
             fontsize=11, y=1.01)
fig.tight_layout()
plt.show()

# ── Group-level bar chart ─────────────────────────────────────────────────────
fig2, ax2 = plt.subplots(figsize=(7, 4))
rng4 = np.random.default_rng(4)
subj_z = zscore_df.groupby(['group', 'subject_id'])['mean_z'].mean().reset_index()

for gi, grp in enumerate(zgrp_order):
    vals = subj_z.query('group == @grp')['mean_z'].to_numpy()
    mean_v = np.nanmean(vals)
    sem_v = np.nanstd(vals) / np.sqrt(len(vals))
    ax2.bar(gi, mean_v, yerr=sem_v, color=zgrp_colors[grp],
            alpha=0.75, capsize=5, width=0.6, label=zgrp_label[grp])
    ax2.scatter(gi + rng4.uniform(-0.18, 0.18, len(vals)), vals,
                color=zgrp_colors[grp], s=40, zorder=3)

ax2.axhline(0, color='gray', lw=0.9, linestyle='--')
ax2.set_xticks(range(len(zgrp_order)))
ax2.set_xticklabels([zgrp_label[g] for g in zgrp_order], fontsize=10)
ax2.set_ylabel('Mean z-score (per-subject mean across trials)')
ax2.set_title('Mean z-score in infusion window — per-subject means')
fig2.tight_layout()
plt.show()

