# Movement Data Assembly Pipeline

## Overview
This document details how `movement_data.pickle` is assembled from DLC (DeepLabCut) tracking data in the `assemble_dlc_data.ipynb` notebook using functions from `extract_behav_parameters.py`.

## Pipeline Steps

### 1. DLC Data Import
**Function:** `read_DLC_csv(stub, dlcfolder)`

- Reads DeepLabCut CSV files for each recording session (stub)
- Parses column headers to create format: `bodypart_x`, `bodypart_y`, `bodypart_likelihood`
- Returns DataFrame with all bodypart coordinates and likelihood scores

**Bodyparts tracked:**
```python
bodyparts = ["nose", "r_ear", "l_ear", "head_base", "back1", "back2",
             "back3", "back4", "tail_base", "tail1", "tail2", "tail3", "tail_tip"]
```

### 2. Low Likelihood Interpolation
**Function:** `interpolate_low_likehood(df, threshold=0.6)`

**Process:**
1. For each bodypart, identify frames where `likelihood < 0.6`
2. Set x and y coordinates to `NaN` for low-confidence frames
3. Apply linear interpolation to fill gaps (using pandas `.interpolate(method='linear', limit_direction='both')`)

**Purpose:** Remove unreliable tracking data and replace with interpolated estimates

### 3. Movement Calculation
**Function:** `calc_bodypart_movement(df, smooth_method="gaussian", smooth_window=10)`

This is the core movement calculation with several steps:

#### 3.1. Frame-to-Frame Distance
For each bodypart:
```python
dx = df[f"{bodypart}_x"].diff()
dy = df[f"{bodypart}_y"].diff()
distance = sqrt(dx² + dy²)
```
- Calculates Euclidean distance traveled between consecutive frames
- Units: pixels per frame

#### 3.2. Aggregation Across Bodyparts
```python
total_movement = mean(all_bodypart_distances)
```
- Averages movement across all 13 bodyparts
- Creates single movement vector representing overall body displacement

#### 3.3. Smoothing (Optional)
**Method:** Gaussian smoothing (default in notebook)
**Parameters:**
- `smooth_method = "gaussian"`
- `smooth_window = 10` frames
- `sigma = smooth_window / 4.0 = 2.5`

**Implementation:**
```python
from scipy.ndimage import gaussian_filter1d
smoothed_values = gaussian_filter1d(values, sigma=2.5)
```

**Purpose:** Remove high-frequency jitter from tracking noise

**Alternative smoothing methods available:**
- `"moving_avg"`: Simple moving average
- `"savgol"`: Savitzky-Golay filter (preserves peaks)

#### 3.4. Normalization to [0, 1]
**Method:** Min-max scaling
```python
normalized = (total_movement - min_val) / (max_val - min_val)
```

**Result:**
- `0` = no movement (minimum detected in entire session)
- `1` = maximum movement detected in session
- All values scaled proportionally between these extremes

**Important:** This normalization is applied **per session** before trial extraction.

### 4. Trial Extraction (Snipping)
**Function:** `get_behav_snips(behav_vector, solenoid_ts)`

**Parameters:**
- `behav_vector`: Normalized movement from step 3 (continuous, ~10 Hz sampling)
- `solenoid_ts`: Timestamps of solenoid openings (trial onsets)

**Snip dimensions:**
```python
start = int(solenoid_ts[i] * 10) - 50  # 5 seconds before
end = int(solenoid_ts[i] * 10) + 150    # 15 seconds after
snip = behav_vector[start:end]          # Total: 20 seconds, 200 frames
```

**Trial structure:**
- Frames 0-49: Baseline (5 seconds pre-trial)
- Frame 50: Solenoid opens (trial start)
- Frames 50-199: Post-trial period (15 seconds)

### 5. Z-Scoring (NOT Applied in movement_data.pickle)

**Important:** The notebook shows z-scoring is **commented out** in the final assembly:

```python
snips_tmp = get_behav_snips(
    solenoid_ts=get_ttls(stub, DATAFOLDER),
    behav_vector=get_movement(stub, DLCFOLDER),
    # zscore_to_entire_snips=True,     # COMMENTED OUT
    # zscore_to_baseline=True          # COMMENTED OUT
)
```

**Available z-scoring options (not used):**
- **`zscore_to_baseline=True`**: `(snip - baseline_mean) / baseline_std`
  - Uses frames 0-49 as baseline for each trial
  - Implemented via `trompy.zscore(snips, baseline_points=50)`
  
- **`zscore_to_entire_snips=True`**: `(snip - session_mean) / session_std`
  - Uses entire session statistics (excluding first/last 100 seconds)
  - Applied uniformly to all trials in session

**Conclusion:** Movement data in `movement_data.pickle` is **NOT z-scored**, only normalized [0,1].

### 6. Metadata Assembly
**Function:** `assemble_snips_and_x(meta_df, DLCFOLDER, get_movement)`

For each recording session:
1. Process all trials (snips)
2. Create metadata DataFrame with:
   - `trial`: Trial number within session
   - `id`: Rat ID
   - `condition`: Experimental condition (deplete/replete/thirsty/replete_exp)
   - `infusiontype`: NaCl concentration (10NaCl or 45NaCl)

### 7. Filtering and Final Export

**Condition filtering (Cell 10):**
```python
mask = (x_array.condition != "thirsty") & (x_array.condition != "replete_exp")
snips_red = snips[mask]
x_array_red = x_array[mask].reset_index(drop=True)
```
- Excludes "thirsty" and "replete_exp" conditions
- Keeps only "deplete" and "replete" trials

**Final pickle contents:**
```python
data_to_save = {
    "snips_movement": snips_red,    # Shape: (n_trials, 200)
    "x_movement": x_array_red,      # DataFrame with trial metadata
}
```

---

## Summary of Processing Parameters

| Processing Step | Parameter | Value | Notes |
|----------------|-----------|--------|-------|
| **Likelihood threshold** | `threshold` | 0.6 | Frames below are interpolated |
| **Smoothing method** | `smooth_method` | "gaussian" | Applied to movement vector |
| **Smoothing window** | `smooth_window` | 10 frames | = 1 second @ 10 Hz |
| **Smoothing sigma** | `sigma` | 2.5 | = window/4 |
| **Normalization** | Method | Min-max [0, 1] | Per session, before snipping |
| **Z-scoring** | Applied | **NO** | Options commented out |
| **Snip duration** | Total | 20 seconds | 200 frames @ 10 Hz |
| **Baseline period** | Frames 0-49 | 5 seconds | Pre-trial |
| **Trial period** | Frames 50-199 | 15 seconds | Post-solenoid |

---

## Key Insights

### What the movement values represent:
- **Raw calculation**: Average pixel displacement per frame across all bodyparts
- **After smoothing**: Gaussian-smoothed to remove jitter
- **After normalization**: Scaled to [0, 1] relative to session min/max
- **Final values**: NOT z-scored, NOT baseline-subtracted
  - `0.0` = least movement detected in that session
  - `1.0` = most movement detected in that session
  - Values are comparable within session but not across sessions

### Why normalized [0, 1]?
- Accounts for session-to-session variability (camera position, zoom, etc.)
- Makes movement comparable across different recording setups
- Provides interpretable scale where threshold (e.g., 0.02) represents 2% of max movement
- Allows setting universal thresholds for "moving" vs "stationary"

### What's NOT in the data:
- ❌ Z-scores (not applied)
- ❌ Baseline correction (not applied)
- ❌ Raw pixel units (normalized away)
- ❌ Angular velocity (that's in `angvel_data_z.pickle`)
- ❌ Individual bodypart movements (averaged together)

### Threshold interpretation:
If using `threshold = 0.02` for "time moving" calculation:
- Movement > 0.02 = moving
- 0.02 represents 2% of maximum movement in that session
- Accounts for fact that max movement varies between sessions

---

## Code References

### Main Functions Used:
1. `read_DLC_csv()` - Line 28, extract_behav_parameters.py
2. `interpolate_low_likehood()` - Line 63, extract_behav_parameters.py  
3. `calc_bodypart_movement()` - Line 136, extract_behav_parameters.py
4. `get_behav_snips()` - Line 290, extract_behav_parameters.py
5. `assemble_snips_and_x()` - Cell 7, assemble_dlc_data.ipynb

### Notebook Structure:
- **Cell 5**: `get_movement()` wrapper function with smoothing parameters
- **Cell 7**: Main assembly loop - processes all sessions
- **Cell 10**: Filtering and export to `movement_data.pickle`

---

## Comparison with New Assembly Pipeline

The newer `assemble_all_data.py` script creates `assembled_data.pickle` which:
- ✅ Still uses same movement calculation (`calc_bodypart_movement`)
- ✅ Still applies gaussian smoothing (window=10)
- ✅ Still normalizes to [0, 1]
- ✅ Now stores processing parameters in metadata
- ✅ Optionally calculates raw pixel movement alongside normalized
- ✅ Integrates photometry, clustering, and transition data in one file

The core movement calculation remains identical.
