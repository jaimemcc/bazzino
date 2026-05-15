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
import tdt
import trompy as tp

import dill

# %%
DATAFOLDER = Path("..//data")
TANKFOLDER = Path("D://TestData//bazzino//to McCutcheon from Paula Bazzino")
DLCFOLDER = TANKFOLDER / "Sodium_Appetite_DLC"

# %%
with open(DATAFOLDER / "angvel_data.pickle", "rb") as f:
    data = dill.load(f)

snips_vel = data["snips_vel"]
x_vel = data["x_vel"]

with open(DATAFOLDER / "snips_data.pickle", "rb") as f:
    data = dill.load(f)

snips_photo = data["snips_all"]
x_array = data["x_array"].reset_index(drop=True)



# %%
df1 = x_vel.iloc[:, :4].reset_index(drop=True)
df2 = x_array.iloc[:, :4].reset_index(drop=True)

# %%

# 1. Check if they are exactly equal (quick check)
if df1.equals(df2):
    print("The DataFrames are identical.")
else:
    print("The DataFrames are NOT identical. Finding differences...")

    # 2. Ensure they have the same shape and columns before detailed comparison
    if df1.shape != df2.shape:
        print(f"Shapes are different: df1 is {df1.shape}, df2 is {df2.shape}")
        # Further analysis might be needed if shapes differ significantly.
    elif not df1.columns.equals(df2.columns):
        print("Column names or order are different.")
        print(f"df1 columns: {df1.columns.tolist()}")
        print(f"df2 columns: {df2.columns.tolist()}")
    else:
        # 3. Perform element-wise comparison
        # This creates a boolean DataFrame of the same shape, True where elements differ
        diff_mask = (df1 != df2) & ~(df1.isnull() & df2.isnull()) # Handle NaNs correctly (NaN != NaN is True)

        # 4. Identify rows where at least one value is different
        # .any(axis=1) checks if any True exists across the columns for each row
        differing_rows_mask = diff_mask.any(axis=1)

        if differing_rows_mask.any():
            print("\nRows that differ in df1:")
            print(df1[differing_rows_mask])
            
            print("\nCorresponding rows in df2:")
            print(df2[differing_rows_mask])
            
            print("\nBoolean mask of differences (True where elements differ):")
            print(diff_mask[differing_rows_mask])
            
            # For a more combined view, you can concatenate them
            # df_all = pd.concat([df1[differing_rows_mask], df2[differing_rows_mask]], keys=['df1', 'df2'], axis=0)
            # df_comparison = df_all.swaplevel(axis=0).sort_index()
            # print("\nSide-by-side comparison of differing rows:")
            # print(df_comparison)

        else:
            print("No content differences found after aligning (could be dtype or specific NaN differences if .equals was False).")

# To investigate dtype differences if .equals() is False but no content diffs are found:

# %%
# ...existing code...
# Cell c9da2892:
# df1 = x_vel.reset_index(drop=True)
# df2 = x_array.reset_index(drop=True)
# ...existing code...

# New cell (e.g., f1d991c6)

# First, ensure df1 and df2 are defined from your previous cells
# For clarity, let's re-state their origin if this cell is run independently
# df1 = x_vel.iloc[:, :4].reset_index(drop=True) 
# df2 = x_array.iloc[:, :4].reset_index(drop=True)

if len(df1) == len(df2):
    print("DataFrames have the same number of rows.")
    print("Use the comparison method from the previous cell to find content differences within rows.")
else:
    # Determine which DataFrame is longer
    if len(df1) > len(df2):
        df_longer = df1
        df_shorter = df2
        longer_name = "df1 (from x_vel)"
        shorter_name = "df2 (from x_array)"
    else:
        df_longer = df2
        df_shorter = df1
        longer_name = "df2 (from x_array)"
        shorter_name = "df1 (from x_vel)"

    print(f"{longer_name} has {len(df_longer) - len(df_shorter)} more rows than {shorter_name}.")

    # Perform a merge to find rows in df_longer that are not in df_shorter.
    # We merge on all columns to ensure rows are identical.
    # It's assumed that df_longer and df_shorter have the same column names for comparison.
    # If column names differ, you'd need to align them or specify left_on and right_on.
    
    # Get the list of columns to merge on (all columns of the DataFrames)
    # This assumes df_longer.columns and df_shorter.columns are identical
    # which they should be after iloc[:, :4] if original DFs had comparable structures.
    merge_columns = df_longer.columns.tolist()

    merged_df = pd.merge(df_longer, df_shorter, on=merge_columns, how='left', indicator=True)

    # Filter for rows that are only in the longer DataFrame
    # These rows will have '_merge' == 'left_only'
    extra_rows_in_longer = merged_df[merged_df['_merge'] == 'left_only']

    # Select only the original columns (drop the '_merge' column)
    extra_rows_content = extra_rows_in_longer[merge_columns]

    if not extra_rows_content.empty:
        print(f"\nContent of the rows present in {longer_name} but not in {shorter_name}:")
        print(extra_rows_content)
    else:
        # This case would be unusual if the lengths are different and common rows are expected to match.
        # It might imply that all rows in the shorter DataFrame are also in the longer one,
        # and the additional rows in the longer DataFrame are exact duplicates of rows that *are* also in the shorter one.
        # Or, it could mean all rows from the longer DataFrame found a match in the shorter one, which contradicts the length difference.
        print(f"No unique extra rows found in {longer_name} using this method. This could indicate that while lengths differ, all rows in the longer DataFrame found a match in the shorter one (possibly due to duplicate rows within the shorter DataFrame that match the 'extra' rows of the longer one), or that the 'extra' rows are not actually unique to the longer DataFrame when compared to the shorter one.")

# %%
# This code should be in a new cell, AFTER the cell that defines 'extra_rows_content'.
# It assumes df1, df2, df_longer, extra_rows_content, x_array (as loaded initially)
# are defined in the notebook's global scope from previous cells.

# --- IMPORTANT: Define your corresponding NumPy array here ---
# Replace 'None' with your actual NumPy array that corresponds to x_array.
# It MUST have the same number of rows as x_array before this cleaning step.
# For example:

snips_data_for_x_array = snips_photo
# -------------------------------------------------------------

if 'extra_rows_content' not in locals() or 'df_longer' not in locals() or 'df2' not in locals():
    print("Error: Essential variables ('extra_rows_content', 'df_longer', 'df2') not found.")
    print("Please ensure the previous cell (that identifies differing rows) has been run successfully.")
else:
    indices_to_delete = extra_rows_content.index # These are indices from df_longer

    # Check if x_array (represented by df2) was the longer DataFrame
    if df_longer is df2: 
        print(f"x_array (via df2) was identified as the longer DataFrame. Proceeding to remove {len(indices_to_delete)} rows.")

        # 1. Create a cleaned version of the DataFrame
        # df_longer is the DataFrame (x_array.reset_index(drop=True)) from which extra rows were identified.
        x_array_cleaned = x_array.drop(indices_to_delete).reset_index(drop=True) 
        
        print(f"Original x_array (represented by df_longer) row count: {len(df_longer)}")
        print(f"Cleaned x_array_cleaned row count: {len(x_array_cleaned)}")
        
        # To update your original 'x_array' variable if desired:
        # x_array = x_array_cleaned
        # print("Notebook variable 'x_array' can be updated to 'x_array_cleaned'.")


        # 2. Remove corresponding rows from the NumPy array
        if snips_data_for_x_array is not None:
            if len(snips_data_for_x_array) == len(df_longer):
                snips_data_for_x_array_cleaned = np.delete(snips_data_for_x_array, indices_to_delete.to_numpy(), axis=0)
                print(f"Original NumPy array 'snips_data_for_x_array' row count: {len(snips_data_for_x_array)}")
                print(f"Cleaned NumPy array 'snips_data_for_x_array_cleaned' row count: {len(snips_data_for_x_array_cleaned)}")

                # To update your original NumPy array variable if desired:
                # snips_data_for_x_array = snips_data_for_x_array_cleaned
                # print("Notebook variable 'snips_data_for_x_array' can be updated to 'snips_data_for_x_array_cleaned'.")
            else:
                print(f"Error: Length of 'snips_data_for_x_array' ({len(snips_data_for_x_array)}) "
                      f"does not match the length of x_array before deletion ({len(df_longer)}). "
                      "Rows not deleted from NumPy array.")
        else:
            print("Warning: 'snips_data_for_x_array' is not defined or is None.")
            print("Cannot delete rows from the corresponding NumPy array. Please define it.")

    elif df_longer is df1: # This means x_vel (represented by df1) was longer
        print(f"x_vel (via df1) was identified as the longer DataFrame. {len(indices_to_delete)} extra rows were found in it.")
        print("The request was to modify 'x_array'. No changes made to 'x_vel' or its corresponding NumPy array in this step.")
    else:
        print("Error: Could not definitively determine if df_longer corresponds to x_array (df2) or x_vel (df1). No rows deleted.")

# %%
with open(DATAFOLDER / "bazzino_data.pickle", "wb") as f:
    dill.dump({
        "x_array": x_array_cleaned,
        "snips_photo": snips_data_for_x_array_cleaned,
        "snips_vel": snips_vel
    }, f)

# %%
with open(DATAFOLDER / "snips_photo_forDLC.pickle", "wb") as f:
    dill.dump(snips_data_for_x_array_cleaned, f)
    print("Cleaned NumPy array 'snips_data_for_x_array_cleaned' saved to disk.")
    
with open(DATAFOLDER / "x_array_forDLC.pickle", "wb") as f:
    dill.dump(x_array_cleaned, f)
    print("Cleaned x_array 'snips_data_for_x_array_cleaned' saved to disk.")

# %%
x_array_cleaned.head()

# %%
x_array_cleaned.shape

# %%
snips_data_for_x_array_cleaned.shape

# %%
