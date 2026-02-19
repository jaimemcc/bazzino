# Findings: Why Realigned Data Looks Different

## Summary
Investigating why the realigned data in `assembled_data.pickle` shows less steep transitions compared to the original `movement_data.pickle`.

## Confirmed: Both Use Same Preprocessing

### Z-Scoring Status (VERIFIED)
Both old and new data:
- ✅ NO z-scoring applied (both options commented out in old, False in new)
- ✅ Only normalized to [0,1]
- ✅ No negative values present
- ✅ Documentation in `movement_data_assembly_pipeline.md` is CORRECT

### Matching Processing Parameters

| Parameter | Old | New | Status |
|-----------|-----|-----|--------|
| DLC likelihood threshold | 0.6 | 0.6 | ✅ Match |
| Smoothing method | gaussian | gaussian | ✅ Match |
| Smoothing window | 10 | 10 | ✅ Match |
| Normalization | [0,1] | [0,1] | ✅ Match |
| Z-score to baseline | False | False | ✅ Match |
| Snip pre-time | 5s (50 bins) | 5s (50 bins) | ✅ Match |
| Snip post-time | 15s (150 bins) | 15s (150 bins) | ✅ Match |

## Other Potential Causes

Since preprocessing is identical, the reduced steepness must come from somewhere else:

## Current Investigation: Snips Differ

**Confirmed:** The movement snips themselves differ between old and new data for deplete + 45NaCl trials.

This means either:
1. **Reordering**: Same trials in both datasets but in different order
2. **Wrong trials**: Different trials being included/excluded
3. **Different values**: Same trial IDs but different snip values (different processing)

### Single-Session Diagnostic Approach

To identify the exact processing parameters used in the old data, we:

1. **Select a test rat** with a successful fit present in both datasets
2. **Reprocess raw DLC data** with all possible parameter combinations:
   - **4 bodypart options**: all (13 bodyparts), head_only, back_only, tail_only
   - **2 weight_by_zscore options**: False (default) or True (weight bodyparts by z-scored movement)
   - **4 processing options**: normalized only, zscore_baseline, zscore_entire, no_normalization
   - **Total: 32 combinations**
3. **Compare each combination** to old data using MSE, MAE, correlation
4. **Identify best match** = exact parameters used in old pipeline

### Key Parameter: weight_by_zscore

The `calc_bodypart_movement` function has a `weight_by_zscore` parameter that:
- When `False` (default): Simply averages movement across all bodyparts
- When `True`: Z-scores each bodypart's movement before aggregating
  - Emphasizes bodyparts that move more relative to their baseline
  - Changes the relative contribution of different body regions

This parameter was not initially tested but could explain differences if the old pipeline used `weight_by_zscore=True`.

### Why This Approach Works

If the best match has MSE < 1e-6, we've found the exact processing used in the original `movement_data.pickle`.

## Additional Checks in Diagnostic Notebook

The notebook also includes checks for:
- **Trial matching**: Identifying which trials are in old vs new datasets
- **Order comparison**: Detecting if trials are just reordered vs actually different
- **Rat-by-rat analysis**: Finding which specific rats have mismatches
- **Transition point comparison**: Verifying fitted x0 values match
- **Realignment visualization**: Comparing cluster transition steepness

## Most Likely Causes (if not processing parameters)

If single-session match is perfect but full dataset still differs:

### 1. Equalization Step Reordering
The new pipeline's `equalize_datasets()` function merges photometry and behavior data:
- Inner merge may reorder trials (pandas merges don't preserve order)
- Could drop trials without matching photometry data
- May create different trial sequences even with same trials present

### 2. Different Trial Selection
- Old pipeline: Movement data only
- New pipeline: Must have matching photometry AND movement
- Trials present in old but missing in new would dilute realignment

### 3. Metadata/Assembly Differences
- Different session processing order (10NaCl then 45NaCl vs reverse)
- Different error handling dropping different sessions
- FileKey CSV differences between old and new runs
