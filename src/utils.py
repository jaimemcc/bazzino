import numpy as np
import trompy as tp

def make_realigned_trials(x_array, fits, verbose=True):

# Realigning behavioral data to dopamine transitions
    realigned_trials = []

    for rat in x_array.id.unique():
        
        x_array_r = x_array.query("id == @rat & condition == 'deplete' & infusiontype == '45NaCl'")

        if rat not in fits.id.unique():
            if verbose:
                print(f"Rat {rat} not found in fitted parameters, skipping.")
            realigned_trials.append([np.nan] * len(x_array_r.trial))
        else:
            transition_point = int(fits.query("id == @rat").x0_orig.values[0])
            realigned_trials.append(x_array_r.trial - transition_point)

    return (
        x_array
        .assign(trial_aligned=tp.flatten_list(realigned_trials))
        .dropna()
        .reset_index(drop=True)
        # .astype(int)  # Ensure the trial_aligned column is of integer type
        )

def get_time_moving(snips, threshold=1):

    moving = []
    for i in range(snips.shape[0]):
        snip = snips[i, 50:150]
        tmp = len([x for x in snip if x > threshold]) / len(snip)
        moving.append(tmp)

    return np.array(moving)


def recalculate_time_moving(df, snips_movement, threshold=0.02, start_bin=50, end_bin=150):
    """
    Recalculate time_moving with a different threshold and update the dataframe.
    
    Parameters
    ----------
    df : pd.DataFrame
        Dataframe with trial information (must have same number of rows as snips_movement)
    snips_movement : np.ndarray
        Movement snips array with shape (n_trials, n_bins)
    threshold : float, default=0.02
        Movement threshold value (normalized movement units)
    start_bin : int, default=50
        Starting bin index for calculation window
    end_bin : int, default=150
        Ending bin index for calculation window
        
    Returns
    -------
    pd.DataFrame
        Copy of df with updated time_moving column
    """
    moving = []
    for i in range(snips_movement.shape[0]):
        snip = snips_movement[i, start_bin:end_bin]
        tmp = len([x for x in snip if x > threshold]) / len(snip)
        moving.append(tmp)
    
    return df.assign(time_moving=np.array(moving))

def scale_vlim_to_data(snips, percentile=99):
    # Get the vlim values for the heatmap based on the data distribution
    vlim = np.percentile(np.abs(snips), percentile)
    return -vlim, vlim

def otsu_threshold(values, nbins=256):
    x = np.asarray(values).ravel()
    x = x[np.isfinite(x)]
    hist, bin_edges = np.histogram(x, bins=nbins)
    bin_centers = 0.5*(bin_edges[:-1] + bin_edges[1:])
    weight1 = np.cumsum(hist).astype(float)
    weight2 = weight1[-1] - weight1
    mean1 = np.cumsum(hist * bin_centers) / np.maximum(weight1, 1)
    mean2 = (np.cumsum((hist * bin_centers)[::-1]) / np.maximum(weight2[::-1], 1))[::-1]
    # Between-class variance:
    var_between = weight1[:-1] * weight2[:-1] * (mean1[:-1] - mean2[:-1])**2
    idx = np.nanargmax(var_between)
    return float(bin_centers[:-1][idx])