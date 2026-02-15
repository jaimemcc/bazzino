
from pathlib import Path
import numpy as np
import pandas as pd
import tdt
from scipy.ndimage import gaussian_filter1d
from scipy.signal import savgol_filter

import trompy as tp

def get_ttls(stub, datafolder, tankfolder=None):

    datafolder = Path(datafolder)
    if tankfolder is not None:
        tankfolder = Path(tankfolder)

    # Get the TTLs
    if (datafolder / "ttls.csv").exists():
        ttls_df = pd.read_csv(datafolder / "ttls.csv")
        if stub in ttls_df.columns:
            sol = ttls_df.loc[:, stub].values
            return sol
    else:
        print("No pre-saved ttls.csv file found. Reading from TDT tank.")
        data = tdt.read_block(tankfolder / stub, evtype=["epocs"])
        sol = data.epocs.sol_.onset
    
    return sol

def read_DLC_csv(stub, dlcfolder):

    date = stub.split("-")[1]
    pattern_str = f"PB_NAapp-{date}_{stub}*.csv"
    
    matching_files = list(dlcfolder.glob(pattern_str))
    
    print(dlcfolder)

    filename = None
    if not matching_files:
        print(f"Error: No DLC file found for stub {stub} with pattern {pattern_str}")
        return None
    
    elif len(matching_files) > 1:
        print(f"Warning: Multiple DLC files found for stub {stub} with pattern {pattern_str}: {matching_files}")
        # Decide how to handle multiple matches: take the first, last, or error.
        # For now, let's take the first one.
        filename = matching_files[0]
        print(f"Using file: {filename}")
    else:
        filename = matching_files[0]
        print(f"Found file: {filename}")
    
    header_df = pd.read_csv(filename, skiprows=1, nrows=2, header=None)
    
    row2_values = header_df.iloc[0].astype(str) # This is original row 2
    row3_values = header_df.iloc[1].astype(str) # This is original row 3
    
    new_column_names = [f"{val2.lower().replace(' ', '')}_{val3}" for val2, val3 in zip(row2_values, row3_values)]
    
    df = pd.read_csv(filename, skiprows=3, header=None, names=new_column_names)
    
    return df

def interpolate_low_likehood(df, threshold=0.5):
    # Work on a copy to avoid mutating the caller's DataFrame
    df = df.copy()
    # Convert all columns to numeric, coercing errors. This is important.
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    # Identify unique bodyparts mentioned in the columns
    bodyparts = set()
    for col_name in df.columns:
        parts = col_name.split('_')
        if len(parts) > 1: # e.g., leftear_x, leftear_likelihood
            bodyparts.add(parts[0]) 
    
    # For each bodypart, interpolate x and y based on likelihood
    for bp in bodyparts:
        x_col = f"{bp}_x"
        y_col = f"{bp}_y"
        likelihood_col = f"{bp}_likelihood"

        if x_col in df.columns and y_col in df.columns and likelihood_col in df.columns:
            # Condition where likelihood is below threshold
            condition = df[likelihood_col] < threshold
            # print(np.sum(condition), "values below threshold for", bp)
            
            # Set x and y to NaN based on the condition
            df.loc[condition, x_col] = np.nan
            df.loc[condition, y_col] = np.nan
            
            # print(np.sum(condition), "values set to NaN for", x_col, y_col)
            
            # Interpolate the x and y columns (linear interpolation by default)
            df[x_col] = df[x_col].interpolate(method='linear', limit_direction='both')
            df[y_col] = df[y_col].interpolate(method='linear', limit_direction='both')
        # else:
            # print(f"Warning: Missing x, y, or likelihood columns for bodypart '{bp}' in {filename.name}")

    return df

def calc_angvel_earsonly(df, rightear="rightear", leftear="leftear", 
                        smooth=True, smooth_sigma=2.0, absolute=True):
    """
    Calculate angular velocity using relative ear positions (legacy version).
    
    Uses the vector from left ear to right ear as the head direction.
    Kept for backward compatibility. Can optionally smooth position data before calculation.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with columns for rightear and leftear positions (x and y)
    rightear : str, default "rightear"
        Name of the right ear bodypart
    leftear : str, default "leftear"
        Name of the left ear bodypart
    smooth : bool, default True
        If True, smooth position data before calculating angles using Gaussian filter.
    smooth_sigma : float, default 2.0
        Standard deviation for Gaussian smoothing. Only used if smooth=True.
    absolute : bool, default True
        If True, return absolute value of angular velocity (magnitude only).
    
    Returns
    -------
    pd.DataFrame
        Original dataframe with added columns including d_angle_deg (absolute if absolute=True)
    """
    # Make a copy to avoid modifying the input
    df_work = df.copy()
    
    # Apply smoothing to position data if requested
    if smooth:
        bodyparts_to_smooth = [rightear, leftear]
        
        for bp in bodyparts_to_smooth:
            x_col = f"{bp}_x"
            y_col = f"{bp}_y"
            
            if x_col in df_work.columns and y_col in df_work.columns:
                x_vals = pd.to_numeric(df_work[x_col], errors='coerce')
                y_vals = pd.to_numeric(df_work[y_col], errors='coerce')
                
                if x_vals.notna().any():
                    x_mask = x_vals.notna()
                    y_mask = y_vals.notna()
                    
                    x_smooth = x_vals.copy()
                    y_smooth = y_vals.copy()
                    
                    x_smooth[x_mask] = gaussian_filter1d(x_vals[x_mask], sigma=smooth_sigma)
                    y_smooth[y_mask] = gaussian_filter1d(y_vals[y_mask], sigma=smooth_sigma)
                    
                    df_work[x_col] = x_smooth
                    df_work[y_col] = y_smooth
    
    # Calculate angular velocity
    result = (
        df_work
        .assign(
            _rel_rightear_x_orig = lambda x_df: x_df[f"{rightear}_x"] - x_df[f"{leftear}_x"],
            _rel_rightear_y_orig = lambda x_df: x_df[f"{rightear}_y"] - x_df[f"{leftear}_y"]
        )
        .assign(
            ear_distance = lambda x_df: np.sqrt(x_df._rel_rightear_x_orig**2 + x_df._rel_rightear_y_orig**2)
        )
        .assign(
            rel_rightear_x = lambda x_df: np.where(x_df.ear_distance >= 90, np.nan, x_df._rel_rightear_x_orig),
            rel_rightear_y = lambda x_df: np.where(x_df.ear_distance >= 90, np.nan, x_df._rel_rightear_y_orig)
        )
        .assign(
            angle_rad = lambda x_df: np.arctan2(x_df.rel_rightear_y, x_df.rel_rightear_x)
        )
        .assign(
            _d_angle_raw = lambda x_df: x_df.angle_rad.diff()
        )
        .assign(
            d_angle = lambda x_df: x_df._d_angle_raw.fillna(0),
            d_angle_wrapped = lambda x_df: (x_df._d_angle_raw + np.pi) % (2 * np.pi) - np.pi
        )
        .assign(
            d_angle_deg = lambda x_df: np.rad2deg(x_df.d_angle_wrapped)
        )
    )
    
    # Apply absolute value if requested
    if absolute:
        result = result.assign(
            d_angle_deg = lambda x_df: np.abs(x_df['d_angle_deg'])
        )
    
    return result.drop(columns=['_rel_rightear_x_orig', '_rel_rightear_y_orig', '_d_angle_raw'], errors='ignore')

def calc_angular_velocity(df, rightear="rightear", leftear="leftear", head_base="head_base",
                         smooth=True, smooth_sigma=2.0, absolute=True):
    """
    Calculate angular velocity using head_base as reference point.
    
    Uses the vector from head_base to the midpoint between the two ears as the head direction.
    This provides a more stable estimate of head orientation based on the overall head position.
    Can optionally smooth position data before calculation to reduce jitter.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with columns for rightear, leftear, and head_base positions (x and y)
    rightear : str, default "rightear"
        Name of the right ear bodypart
    leftear : str, default "leftear"
        Name of the left ear bodypart
    head_base : str, default "head_base"
        Name of the head base (center of head) bodypart
    smooth : bool, default True
        If True, smooth position data before calculating angles using Gaussian filter.
        Smoothing reduces DLC jitter and leads to cleaner angular velocity signals.
    smooth_sigma : float, default 2.0
        Standard deviation for Gaussian smoothing. Only used if smooth=True.
        Higher values = more smoothing. Typical range: 0.5-3.0
    absolute : bool, default True
        If True, return absolute value of angular velocity (magnitude only).
        If False, return signed angular velocity (preserves direction).
    
    Returns
    -------
    pd.DataFrame
        Original dataframe with added columns:
        - ear_distance: distance between the two ears
        - ear_midpoint_x, ear_midpoint_y: midpoint between ears
        - rel_head_x, rel_head_y: head direction vector (from head_base to ear midpoint)
        - angle_rad: head angle in radians
        - d_angle_deg: frame-to-frame angular velocity in degrees (absolute if absolute=True)
    """
    # Make a copy to avoid modifying the input
    df_work = df.copy()
    
    # Apply smoothing to position data if requested
    if smooth:
        bodyparts_to_smooth = [rightear, leftear, head_base]
        
        for bp in bodyparts_to_smooth:
            x_col = f"{bp}_x"
            y_col = f"{bp}_y"
            
            if x_col in df_work.columns and y_col in df_work.columns:
                x_vals = pd.to_numeric(df_work[x_col], errors='coerce')
                y_vals = pd.to_numeric(df_work[y_col], errors='coerce')
                
                # Only smooth non-NaN values
                if x_vals.notna().any():
                    x_mask = x_vals.notna()
                    y_mask = y_vals.notna()
                    
                    # Apply Gaussian filter only to valid values
                    x_smooth = x_vals.copy()
                    y_smooth = y_vals.copy()
                    
                    x_smooth[x_mask] = gaussian_filter1d(x_vals[x_mask], sigma=smooth_sigma)
                    y_smooth[y_mask] = gaussian_filter1d(y_vals[y_mask], sigma=smooth_sigma)
                    
                    df_work[x_col] = x_smooth
                    df_work[y_col] = y_smooth
    
    # Calculate angular velocity
    result = (
        df_work
        .assign(
            # Calculate midpoint between the two ears
            ear_midpoint_x = lambda x_df: (x_df[f"{rightear}_x"] + x_df[f"{leftear}_x"]) / 2,
            ear_midpoint_y = lambda x_df: (x_df[f"{rightear}_y"] + x_df[f"{leftear}_y"]) / 2
        )
        .assign(
            # Calculate distance between ears (for data quality check)
            ear_distance = lambda x_df: np.sqrt(
                (x_df[f"{rightear}_x"] - x_df[f"{leftear}_x"])**2 + 
                (x_df[f"{rightear}_y"] - x_df[f"{leftear}_y"])**2
            )
        )
        .assign(
            # Vector from head_base to ear midpoint (head direction)
            _rel_head_x_orig = lambda x_df: x_df.ear_midpoint_x - x_df[f"{head_base}_x"],
            _rel_head_y_orig = lambda x_df: x_df.ear_midpoint_y - x_df[f"{head_base}_y"]
        )
        .assign(
            # Filter out unlikely values (e.g., ear distance too large)
            rel_head_x = lambda x_df: np.where(x_df.ear_distance >= 90, np.nan, x_df._rel_head_x_orig),
            rel_head_y = lambda x_df: np.where(x_df.ear_distance >= 90, np.nan, x_df._rel_head_y_orig)
        )
        .assign(
            # Calculate angle from head_base to ear midpoint
            angle_rad = lambda x_df: np.arctan2(x_df.rel_head_y, x_df.rel_head_x)
        )
        .assign(
            # Calculate frame-to-frame angle changes
            _d_angle_raw = lambda x_df: x_df.angle_rad.diff()
        )
        .assign(
            d_angle = lambda x_df: x_df._d_angle_raw.fillna(0),
            d_angle_wrapped = lambda x_df: (x_df._d_angle_raw + np.pi) % (2 * np.pi) - np.pi
        )
        .assign(
            # Convert to degrees
            d_angle_deg = lambda x_df: np.rad2deg(x_df.d_angle_wrapped)
        )
    )
    
    # Apply absolute value if requested
    if absolute:
        result = result.assign(
            d_angle_deg = lambda x_df: np.abs(x_df['d_angle_deg'])
        )
    
    return result.drop(columns=['_rel_head_x_orig', '_rel_head_y_orig', '_d_angle_raw'], errors='ignore')
    
def calc_bodypart_movement(df, weight_by_zscore=False, smooth_method=None, smooth_window=5, 
                          include_bodyparts=None, exclude_bodyparts=None, normalize=True):
    """
    Calculate overall bodypart movement across time.
    
    Computes frame-to-frame displacement for each bodypart (in pixels), optionally weights
    by z-scored movement, applies smoothing to remove jitter, and optionally normalizes 
    to [0, 1] range.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with columns like 'bodypart_x', 'bodypart_y', 'bodypart_likelihood'
    weight_by_zscore : bool, default False
        If True, each bodypart's movement is z-scored before aggregation.
        This emphasizes bodyparts that move more than usual relative to their
        own movement baseline.
    smooth_method : str, optional
        Smoothing method to apply. Options:
        - None: no smoothing
        - 'gaussian': Gaussian blur (recommended for general use)
        - 'moving_avg': Simple moving average
        - 'savgol': Savitzky-Golay filter (preserves peaks, requires odd window)
    smooth_window : int, default 5
        Window size for smoothing. Should be odd for 'savgol'.
        If even with 'savgol', will be incremented by 1.
    include_bodyparts : list of str, optional
        If provided, only these bodyparts will be included in the calculation.
        Example: ['rightear', 'leftear']. If None, all detected bodyparts are used.
    exclude_bodyparts : list of str, optional
        If provided, these bodyparts will be excluded from the calculation.
        Example: ['rightear']. Cannot be used together with include_bodyparts.
    normalize : bool, default True
        If True, normalize movement to [0, 1] range using min-max scaling.
        If False, return raw average pixels per bodypart per frame.
    
    Returns
    -------
    pd.Series
        Movement for each frame. If normalize=True (default), values are in [0, 1].
        If normalize=False, values are in pixels (average across bodyparts).
    
    Raises
    ------
    ValueError
        If both include_bodyparts and exclude_bodyparts are provided.
    
    Notes
    -----
    - Bodyparts are identified by parsing column names for patterns like
      'bodypart_x' and 'bodypart_y'.
    - Likelihood columns are ignored.
    - Frame-to-frame distance is calculated as sqrt(dx^2 + dy^2).
    - Values are averaged across bodyparts.
    - With smoothing, edge frames may be affected depending on method.
    """
    if include_bodyparts is not None and exclude_bodyparts is not None:
        raise ValueError("Cannot specify both include_bodyparts and exclude_bodyparts")
    
    df = df.copy()
    
    # Convert all columns to numeric
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    
    # Extract unique bodyparts by finding x/y column pairs
    bodyparts = set()
    for col_name in df.columns:
        if col_name.endswith('_likelihood'):
            continue
        parts = col_name.split('_')
        if len(parts) >= 2 and parts[-1] in ('x', 'y'):
            # Reconstruct bodypart name (handles names with underscores)
            bodypart = '_'.join(parts[:-1])
            bodyparts.add(bodypart)
    
    # Filter bodyparts based on include/exclude lists
    if include_bodyparts is not None:
        bodyparts = bodyparts.intersection(set(include_bodyparts))
    elif exclude_bodyparts is not None:
        bodyparts = bodyparts.difference(set(exclude_bodyparts))
    
    # Calculate movement for each bodypart
    bodypart_movements = []
    for bp in sorted(bodyparts):
        x_col = f"{bp}_x"
        y_col = f"{bp}_y"
        
        if x_col in df.columns and y_col in df.columns:
            dx = df[x_col].diff()
            dy = df[y_col].diff()
            distance = np.sqrt(dx**2 + dy**2)
            
            if weight_by_zscore:
                # Z-score this bodypart's movement
                mean = distance.mean()
                std = distance.std()
                if std > 0:
                    distance = (distance - mean) / std
                else:
                    distance = pd.Series(0, index=df.index)
            
            bodypart_movements.append(distance)
    
    # Aggregate: average movement across bodyparts
    if bodypart_movements:
        total_movement = pd.concat(bodypart_movements, axis=1).mean(axis=1)
        
        # Apply smoothing if requested
        if smooth_method is not None:
            values = total_movement.values
            
            if smooth_method == 'gaussian':
                # Gaussian smoothing
                sigma = smooth_window / 4.0  # Convert window to sigma
                smoothed_values = gaussian_filter1d(values, sigma=sigma)
                total_movement = pd.Series(smoothed_values, index=total_movement.index)
            
            elif smooth_method == 'moving_avg':
                # Simple moving average
                total_movement = total_movement.rolling(
                    window=smooth_window, 
                    center=True,
                    min_periods=1
                ).mean()
            
            elif smooth_method == 'savgol':
                # Savitzky-Golay filter (preserves peaks better)
                # Make sure window is odd
                window = smooth_window if smooth_window % 2 == 1 else smooth_window + 1
                # Polyorder should be less than window size
                polyorder = min(3, window - 2)
                smoothed_values = savgol_filter(values, window, polyorder, mode='nearest')
                total_movement = pd.Series(smoothed_values, index=total_movement.index)
        
        # Normalize to [0, 1] using min-max scaling if requested
        if normalize:
            min_val = total_movement.min()
            max_val = total_movement.max()
            
            if max_val > min_val:
                normalized = (total_movement - min_val) / (max_val - min_val)
            else:
                # All values are the same, return zeros
                normalized = pd.Series(0.0, index=total_movement.index)
            return normalized
        else:
            # Return raw pixel movement
            return total_movement
    else:
        # No valid bodyparts found
        if normalize:
            return pd.Series(0.0, index=df.index)
        else:
            return pd.Series(0.0, index=df.index)

def get_behav_snips(behav_vector,
                    solenoid_ts,
                    zscore_to_baseline=False,
                    zscore_to_entire_snips=False
                    ):

    snips = []
    for i in range(len(solenoid_ts)-1):
        start = int(solenoid_ts[i] * 10) - 50
        end = int(solenoid_ts[i] * 10) + 150
        snips.append(behav_vector[start:end])
        
    snips = np.array(snips)
    
    if zscore_to_baseline:
        snips = tp.zscore(snips, baseline_points=50)
    elif zscore_to_entire_snips:
        snips = (snips - np.mean(behav_vector[1000:-1000])) / np.std(behav_vector[1000:-1000])
    else:
        pass
    
    return snips



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