# Feature and SHAP interpretation notes

This note captures the current working interpretation of the engineered features and the SHAP-based analysis. It is meant to be read alongside the feature extractor in [src/hybrid_feature_extractor.py](src/hybrid_feature_extractor.py).

## 1. Feature families and what they represent

### Body geometry
- Features such as `Mouse_nose_to_tail`, `Mouse_head_to_tail`, and `Mouse_Ear_distance`.
- These are direct Euclidean distances between pairs of landmarks in a single frame.
- Larger values mean the landmarks are farther apart, which usually corresponds to a more stretched or spatially extended body configuration.
- Smaller values mean a more compact posture.

### Body-part movement
- Features such as `Movement_mouse_nose`, `Movement_mouse_tail_base`, `Movement_mouse_left_ear`, `Movement_mouse_right_ear`, and `Movement_mouse_head_base`.
- These are frame-to-frame displacements of individual body parts.
- They are computed as the distance moved between the previous frame and the current frame.
- Larger values mean the body part moved more between consecutive frames.
- Smaller values mean little or no frame-to-frame displacement.

### Aggregate movement
- Features such as `Total_movement_all_bodyparts_M1` and `Total_movement_M1_*`.
- These summarize movement across several body parts.
- The base aggregate is the sum of the five individual movement features.
- Larger values mean more overall movement across the tracked body.
- These are useful as an overall movement proxy, but they are not equivalent to any single landmark’s movement.

### Hull geometry
- Features derived from pairwise distances between points in the animal body representation.
- Examples include `M1_largest_euclidean_distance_hull`, `M1_smallest_euclidean_distance_hull`, `M1_mean_euclidean_distance_hull`, and `M1_sum_euclidean_distance_hull`.
- These capture the overall spread or geometry of the body in the frame.
- Larger values mean a broader or more stretched body configuration.
- They are mostly about body shape and spread, not about motion per se.

### Rolling-window summaries
- Features such as `Tail_base_movement_M1_mean_2`, `Tail_base_movement_M1_sum_5`, `Head_base_movement_M1_median_7.5`, `Nose_movement_M1_mean_10`, and similar variants.
- These are smoothed summaries of the underlying movement or geometry signal over a short window.
- `mean` and `median` capture the typical level over the window.
- `sum` captures accumulated signal over the window and is more sensitive to sustained movement or sustained spread.
- A larger rolling-window value usually means that the underlying feature was elevated over the recent window, not just at a single frame.

### Derived scores
- Features ending in `_deviation` or `_percentile_rank`.
- These are not raw movement values; they describe how a value compares with its own average or with the overall distribution.
- For the deviation features used here, the implementation is reverse-coded as `average - current_value`. That means a larger deviation value corresponds to a current value that is lower than its recent/typical average, while a smaller deviation value corresponds to a current value that is higher than its average.
- For the current extractor, some columns labeled `_percentile_rank` are actually centered deviations from the mean rather than true percentile ranks. These should be interpreted as relative-to-average features rather than direct movement intensity measures.

## 2. What the SHAP analysis has told us so far

### Overall pattern
- The most important features tend to come from movement-related and geometry-related families.
- The strongest signals are not only the raw movement features, but also smoothed and derived versions of them.

### Appetitive == 1 subset
- Restricting the analysis to Appetitive == 1 makes the interpretation cleaner than using the full pooled dataset.
- In this subset, the strongest positive SHAP patterns were associated with movement-derived and smoothed features, and with some posture or geometry features.

### Quartile-based interpretation
- A quartile-based analysis is more informative than a simple positive/negative SHAP split because it shows how SHAP changes across the raw feature distribution.
- In the Appetitive == 1 subset, features such as `Total_movement_M1_mean_10_percentile_rank`, `Total_movement_M1_mean_7.5_percentile_rank`, `Total_movement_M1_mean_6_percentile_rank`, `Total_movement_M1_mean_5_percentile_rank`, and `Total_movement_all_bodyparts_deviation` showed an increasing trend from low to high raw-feature quartiles.
- Because these are reverse-coded relative-to-average features, the upward trend should be read as: when the current movement value is lower than the recent or typical average, the corresponding feature becomes larger and the model shifts toward appetitive.
- Other features, especially some raw cumulative movement and tortuosity variables, showed the opposite trend, meaning that higher raw sustained movement values are associated with lower appetitive scores.

## 3. Important caveats

- SHAP importance does not by itself tell us whether higher or lower values are driving the prediction; it must be interpreted together with the raw feature values.
- Some derived features are not direct movement measures. They are better interpreted as relative deviations or rank-based contrasts.
- The interpretation is clearer when we condition on the class of interest (for example, Appetitive == 1) rather than mixing classes.
- A larger raw value does not always mean “more movement” for every feature family. For geometry features it can mean more body spread; for deviation features it can mean “more above average”; and for rolling-window sums it can mean more accumulated signal over time.

## 4. Working interpretation summary

At the current stage, the most defensible summary is:
- posture and geometry features are important,
- movement and smoothed movement summaries are also important,
- the model is most clearly associated with lower-than-average movement in the reverse-coded relative-to-average features, while higher raw sustained movement values are associated with lower appetitive scores,
- and the exact direction of effect depends on the feature family and on whether the analysis is conditioned on a specific behavioral class.

## 5. Open questions

- Should we treat deviation and relative-position features as distinct from raw movement features in the final interpretation?
- Should we emphasize the Appetitive == 1 subset as the main scientific result?
- Would it be useful to produce a compact figure showing the top features from each family?
