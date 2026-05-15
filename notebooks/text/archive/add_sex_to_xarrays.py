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
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
import seaborn as sns
import trompy as tp

from scipy import stats

import dill


savefigs = True

DATAFOLDER = Path("..//data")

# %%
with open(DATAFOLDER / "x_array_2clusters_alltrials.pickle", "rb") as f:
    x_array = dill.load(f)

# %%
subjects_10NaCl = pd.read_csv(DATAFOLDER / "10NaCl_SubjectKey.csv")
subjects_45NaCl = pd.read_csv(DATAFOLDER / "45NaCl_SubjectKey.csv")

# %%
subject_df = (pd.concat([subjects_10NaCl.iloc[:, :2], subjects_45NaCl.iloc[:, :2]], axis=0)
              .reset_index()
              .rename(columns={"Subject": "id",
                               "Sex": "sex"})
              .drop(columns=["index"])
)

# %%
x_array = pd.merge(x_array, subject_df[['id', 'sex']], on='id', how='left')

# %%
with open(DATAFOLDER / "x_array_2clusters_alltrials.pickle", "rb") as f:
    x_array = dill.load(f)

# %%
