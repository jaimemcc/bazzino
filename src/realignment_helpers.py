import pandas as pd

def get_realigned_data(df, realignment_column, rats_to_exclude=[], verbose=True):
    
    subset_aligned = (
        df
        .query("condition == 'deplete' & infusiontype == '45NaCl'")
        .dropna(subset=[realignment_column])
        .reset_index(drop=True)
        .sort_values(['id', 'trial'])
        .query("id not in @rats_to_exclude")
    )
    
    # Get animal info
    animals = sorted(subset_aligned['id'].unique())
    n_required = len(animals)
    if verbose:
        print(f"Animals with both transitions and deplete+45NaCl trials: {animals}")
        print(f"Number of animals: {n_required}")
        
        # Show summary statistics
        print(f"\nTrial counts per animal:")
        for animal in animals:
            n_trials = len(subset_aligned.query("id == @animal"))
            print(f"  {animal}: {n_trials} trials")

    return subset_aligned

def only_keep_complete_trials(df, realignment_column, verbose=False):
    # Get the number of unique animals in the dataset
    n_animals = df.id.nunique()
    
    # Group by realignment column and filter groups that have all animals represented
    df_complete = (
        df
        .groupby(realignment_column, group_keys=False)
        .filter(lambda g: len(g) == n_animals)
    )
                
    if verbose:
        print(f"\nComplete realigned trial rows: {len(df_complete)}")
        
        print(df_complete
                .query(f"{realignment_column} == 0")
                .loc[:, ["id", "trial", realignment_column]]) 
    
    return df_complete
    
    