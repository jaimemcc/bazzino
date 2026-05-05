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
            
    # Keep only realigned trials where all selected rats are represented
    subset_aligned_complete = (
        subset_aligned
        .groupby(realignment_column, group_keys=False)
        .filter(lambda g: len(g) == n_required)
        .copy()
    )
    
    if verbose:
        print(f"\nComplete realigned trial rows: {len(subset_aligned_complete)}")
        
        print(subset_aligned_complete
                .query(f"{realignment_column} == 0")
                .loc[:, ["id", "trial", realignment_column]]) 
    
    return subset_aligned_complete
    
    