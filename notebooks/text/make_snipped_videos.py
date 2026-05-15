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
import sys

# Register dill/pathlib compatibility shim BEFORE importing dill
sys.path.insert(0, str(Path("../src").resolve()))
from pickle_compat import enable_dill_pathlib_compat
enable_dill_pathlib_compat()

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import tdt
import trompy as tp

import dill

# %%
DATAFOLDER = Path("..//data")
TANKFOLDER = Path("D://TestData//bazzino//from_paula")
DLCFOLDER = TANKFOLDER / "Sodium_Appetite_DLC"
RESULTSFOLDER = Path("..//results")


# %%
def get_ttls(stub):

    # Get the TTLs
    data = tdt.read_block(TANKFOLDER / stub, evtype=["epocs"])
    sol = data.epocs.sol_.onset
    
    return sol


# %%
## For testing purposes
# stub = "PB71-221123-113609"
# get_ttls(stub)

# df = pd.DataFrame(data=get_ttls(stub), columns=[stub])

# %%
# also make accompanying df that contains details of the rat and the condition and the time in session

def assemble_ttls(csv_path, tank_folder):
    metadata = pd.read_csv(csv_path)
    
    dfs_all = []
    
    for row in metadata.iterrows():
        stub = row[1]["Folder"]
        
        print(stub, len(dfs_all))
        try:
            df_tmp = pd.DataFrame(data=get_ttls(stub), columns=[stub])        
            dfs_all.append(df_tmp)
     
        except FileNotFoundError:
            print("Error with tank for", row[1]["Subject"], row[1]['Physiological state'])
    
    df_concat = pd.concat(dfs_all, axis=1)
    
    return df_concat

dfs_10 = assemble_ttls(DATAFOLDER / "10NaCl_FileKey.csv", TANKFOLDER)
dfs_45 = assemble_ttls(DATAFOLDER / "45NaCl_FileKey.csv", TANKFOLDER)

df_ttls = pd.concat([dfs_10, dfs_45], axis=1)


# %%
df_ttls.to_csv(DATAFOLDER / "ttls.csv")

# %%
import cv2

def save_video_snippet(input_path, output_path, start_frame, end_frame, text=None, text_start_frame=100):
    """
    Extract a portion of an AVI file and save it.
    
    Parameters:
    - input_path: path to input AVI file
    - output_path: path to output AVI file
    - start_frame: starting frame number (0-indexed)
    - end_frame: ending frame number (exclusive)
    - text: text to overlay on the video (optional)
    - text_start_frame: frame number to start displaying text (relative to start_frame)
    """
    # Convert Path objects to strings
    input_path = str(input_path)
    output_folder = output_path / stub
    if not output_folder.exists():
        output_folder.mkdir(parents=True, exist_ok=True)
        
    output_path = str(output_path / stub / "{}_{}.avi".format(stub, start+10))
    print(output_path)

    # Open the input video
    cap = cv2.VideoCapture(input_path)
    
    if not cap.isOpened():
        print(f"Error: Could not open input video: {input_path}")
        return
    
    # Get video properties
    fps = cap.get(cv2.CAP_PROP_FPS)
    width = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
    height = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))
    
    print(f"Input video: {width}x{height} @ {fps} fps")
    
    # Define codec and create VideoWriter
    # Try XVID codec (more compatible than MJPEG for AVI)
    fourcc = cv2.VideoWriter_fourcc(*'XVID')
    out = cv2.VideoWriter(output_path, fourcc, fps, (width, height))
    
    if not out.isOpened():
        print(f"Error: Could not create VideoWriter. Trying alternative codec...")
        # Fallback to MJPEG
        fourcc = cv2.VideoWriter_fourcc(*'MJPG')
        out = cv2.VideoWriter(output_path, fourcc, fps, (width, height))
        if not out.isOpened():
            print("Error: Could not create VideoWriter with either codec")
            cap.release()
            return
    
    # Set starting frame
    cap.set(cv2.CAP_PROP_POS_FRAMES, start_frame)
    
    # Extract and write frames
    frame_count = start_frame
    snippet_frame = 0
    while frame_count < end_frame:
        ret, frame = cap.read()
        if not ret:
            break
        
        # Add text if specified and frame count exceeds text_start_frame
        if text is not None and snippet_frame >= text_start_frame:
            cv2.putText(frame, text, (50, 50), cv2.FONT_HERSHEY_SIMPLEX, 
                        1, (0, 255, 0), 2)
        
        out.write(frame)
        frame_count += 1
        snippet_frame += 1
    
    # Release everything
    cap.release()
    out.release()
    
    print(f"Video saved: {output_path} ({end_frame - start_frame} frames)")

# Example usage:

def get_videopath(stub):
    
    date = stub.split("-")[1]
    
    return "PB_NAapp-" + date + "_" + stub + "_Cam1.avi"

stub = "PB71-221123-113609"
ttls = get_ttls(stub)
start = int(ttls[0] - 10) * 10
end = int(ttls[0] + 10) * 10

save_video_snippet(TANKFOLDER / stub / get_videopath(stub), RESULTSFOLDER, start_frame=start, end_frame=end, text="Infusion")


# %%
for ttl in ttls:
    start = int(ttl - 10) * 10
    end = int(ttl + 10) * 10
    
    save_video_snippet(TANKFOLDER / stub / get_videopath(stub), RESULTSFOLDER, start_frame=start, end_frame=end, text="Infusion")
    
    

# %%
