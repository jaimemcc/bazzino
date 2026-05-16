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
# %load_ext autoreload
# %autoreload 2

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

import seaborn as sns
# import tdt
import trompy as tp

import dill

from extract_behav_parameters import get_ttls, read_DLC_csv, interpolate_low_likehood, calc_angular_velocity, calc_bodypart_movement, get_behav_snips, smooth_array

colors = ["#67AFD2", "#016895", "#F4795B", "#C74632"]
custom_cmap = LinearSegmentedColormap.from_list("custom_diverging", [colors[1], "white", colors[3]])


# %%
DATAFOLDER = Path("..//data")
TANKFOLDER = Path("D://TestData//bazzino//from_paula")
DLCFOLDER = TANKFOLDER / "Sodium_Appetite_DLC" #office computer, Paula's tracking
DLCFOLDER = Path("D:/TestData/bazzino/output_csv_shuffle4") #office computer, Dionne's tracking
# DLCFOLDER = Path("C:/Users/jmc010/Data/bazzino/Output DLC shuffle 4 csv files") # laptop

# %%
stub = "PB71-221123-113609"
# get_ttls(stub, DATAFOLDER, TANKFOLDER)
dlcdata = read_DLC_csv(stub, DLCFOLDER)
df = interpolate_low_likehood(dlcdata, threshold=0.6)
# df = calc_angular_velocity(df, rightear="r_ear", leftear="l_ear")
calc_bodypart_movement(df,
                       normalize=False,
                       calibration_factor=3
                       )

# %%
dlcdata.columns
bodyparts = ["nose", "r_ear", "l_ear", "head_base", "back1", "back2",
             "back3", "back4", "tail_base", "tail1", "tail2", "tail3", "tail_tip"]


# %%
def get_angular_velocity(stub, dlcfolder, zscore=False):
    
    df = read_DLC_csv(stub, dlcfolder) 
    df = interpolate_low_likehood(df, threshold=0.5)
    df = calc_angular_velocity(df, rightear="r_ear", leftear="l_ear")
    
    return df.d_angle_deg

def get_movement(stub, dlcfolder, bodyparts=None):
    
    df = read_DLC_csv(stub, dlcfolder)
    df = interpolate_low_likehood(df, threshold=0.6)

    if bodyparts is None:
        movement = calc_bodypart_movement(df, 
                                        # weight_by_zscore=True,
                                        smooth_method="gaussian", smooth_window=10,
                                        )
    else:
        movement = calc_bodypart_movement(df, 
                                        # weight_by_zscore=True,
                                        include_bodyparts=bodyparts,
                                        smooth_method="gaussian", smooth_window=10,
                                        )
    
    return movement


# %%
## For testing purposes
stub = "PB71-221123-113609" 
# stub = "PB75-221123-140659" #0.45, deplete
stub = "PB61-221024-110730" #0.10, deplete

stub = "PB26-220608-113244" #0.45, deplete
stub = "PB75-221123-140659" #0.45, deplete
stub = "PB73-221123-131413" #0.45, deplete
stub = "PB48-220926-121542" #0.45, deplete
stub = "PB26-220608-113244" #0.45, deplete
stub = "PB71-221123-113609" #0.45, deplete
stub = "PB26-220608-113244" #0.45, deplete

snips = get_behav_snips(solenoid_ts=get_ttls(stub, DATAFOLDER, TANKFOLDER),
                        behav_vector=get_movement(stub, DLCFOLDER),
                        # zscore_to_entire_snips=True,
                        # zscore_to_baseline=True
)

def scale_vlim_to_data(snips, percentile=99):
    # Get the vlim values for the heatmap based on the data distribution
    vlim = np.percentile(np.abs(snips), percentile)
    return -vlim, vlim

pc = 90

def get_time_moving(snips, threshold=1):

    moving = []
    for i in range(snips.shape[0]):
        snip = snips[i, 50:150]
        tmp = len([x for x in snip if x > threshold]) / len(snip)
        moving.append(tmp)

    return np.array(moving)

f, ax = plt.subplots(figsize=(6, 2), ncols=3)
sns.heatmap(snips,
            vmin=scale_vlim_to_data(snips, percentile=pc)[0],
            vmax=scale_vlim_to_data(snips, percentile=pc)[1],
            ax=ax[0],
            cmap=custom_cmap
            )

ax[1].plot(np.mean(snips[:,50:150], axis=1))
ax[2].plot(get_time_moving(snips, threshold=0.018))


# %%
def tweak_x_array(x_array):

    df = pd.concat(x_array, axis=0)

    return (df
            .replace({"condition": {"Sodium Depleted": "deplete",
                                    "Sodium Replete": "replete",
                                    "Sodium Replete Experienced": "replete_exp",
                                    "Thirsty": "thirsty",                                    
                                    },
                      "infusiontype": {45: "45NaCl",
                                       1: "10NaCl"
                                       }
                      })
            .reset_index(drop=True)
            )

def assemble_snips_and_x(meta_df, dlc_folder, snips_function):
    # metadata = meta_df.iloc[28:31, :]
    metadata = meta_df
    
    snips_array, x_array = [], []
    
    for row in metadata.iterrows():
        stub = row[1]["Folder"]
        
        print(stub)
        try:
            snips_tmp = get_behav_snips(solenoid_ts=get_ttls(stub, DATAFOLDER),
                        behav_vector=snips_function(stub, DLCFOLDER),
                        # zscore_to_entire_snips=True,
                        zscore_to_baseline=True
                        )
            nsnips = len(snips_tmp)
            print(nsnips)
            snips_array.append(snips_tmp)
            
            x_array.append(pd.DataFrame(data={"trial": np.arange(nsnips),
                                                "id": row[1]["Subject"],
                                                "condition": row[1]["Physiological state"],
                                                "infusiontype": row[1]["TreatNum"]
                                                }
                                        )
                            )
        except:
            print("Error with tank for", row[1]["Subject"], row[1]['Physiological state'])
            
    snips_array = np.concatenate(snips_array, axis=0)
    x = tweak_x_array(x_array)

    return snips_array, x

meta_df = pd.concat([pd.read_csv(DATAFOLDER / "10NaCl_FileKey.csv"),
                      pd.read_csv(DATAFOLDER / "45NaCl_FileKey.csv")]
                    )

subjects_df = (pd
               .concat([pd.read_csv(DATAFOLDER / "10NaCl_SubjectKey.csv").iloc[:, :2],
                         pd.read_csv(DATAFOLDER / "45NaCl_SubjectKey.csv").iloc[:, :2]]
                        )
               .reset_index()
               .rename(columns={"Subject": "id",
                               "Sex": "sex"})
               .drop(columns=["index"])
)

# snips, x_array = assemble_snips_and_x(meta_df, DLCFOLDER, get_movement)
snips, x_array = assemble_snips_and_x(meta_df, DLCFOLDER, get_angular_velocity)
x_array = pd.merge(x_array, subjects_df[['id', 'sex']], on='id', how='left')


# %%
# for testing

data_to_save = {"snips_ang_vel_z": snips,
                "x_array": x_array,
                }

with open(DATAFOLDER / "angvel_data_z.pickle", "wb") as f:
    dill.dump(data_to_save, f)



# %%
# to make datasets for each body part

data_to_save = {}
data_to_save["bodyparts"] = bodyparts
for bodypart in bodyparts:
    snips, x_array = assemble_snips_and_x(meta_df, DLCFOLDER, lambda stub, dlcfolder: get_movement(stub, dlcfolder, bodyparts=[bodypart]))
    data_to_save[bodypart] = snips

data_to_save["x_array"] = x_array

with open(DATAFOLDER / "bodypart_data.pickle", "wb") as f:
    dill.dump(data_to_save, f)


# %%
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
