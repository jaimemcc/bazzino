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
import tdt

# %%
TANK = Path("C:/Users/jmc010/Data/bazzino_data/tanks/PB23-220608-131619")

# %%
data = tdt.read_block(TANK)

# %%
data.epocs.keys()

# %%
import trompy as tp

blue = data.streams['x65A'].data
uv = data.streams['x05A'].data
fs = data.streams['x65A'].fs

sol = data.epocs["sol_"].onset
cam = data.epocs["Cam1"].onset

filtered = tp.processdata(blue, uv)

# %%
snips, _ = tp.snipper(filtered,
                        sol,
                        fs=fs,
                        baseline_length=5,
                        trial_length=20,
                        bins=20)

# %%
import matplotlib.pyplot as plt

f, ax = plt.subplots()
ax.plot(snips.mean(axis=0))

ax.axvspan(5, 8, color="red", alpha=0.1)



# %%
import numpy as np
import matplotlib.pyplot as plt

bins = np.arange(0.095, 0.105, 0.001)
plt.hist(np.diff(cam), bins=bins)

# %%
