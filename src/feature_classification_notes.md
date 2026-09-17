# Feature classification logic (src/assemble_shap_dfs.py)

This documents the rules used to classify each SHAP feature into a semantic
`group`, `bodypart`, `timewindow`, and `reverse_coded` flag, so they can be
reviewed and revised. After you mark up changes, tell me what to change and
I'll update `assemble_shap_dfs.py` to match.

## 0. Duplicate-column removal (`_find_duplicate_percentile_rank_columns`)

**Resolved:** 15 geometry `_percentile_rank` columns were found to be exact
duplicates of their `_deviation` counterpart (both use the same
`mean - current` formula, per `hybrid_feature_extractor.py`), so they're now
dropped from `df_shap`/`df_shap_raw` before any summary is built. This
removed 15 of the original 168 features (168 → 153), and reduced the
`geometry, reverse coded` count from 33 → 18. `Sum_probabilities_deviation`
vs `Sum_probabilities_percentile_rank` were checked too and are **not**
duplicates (kept as-is).

## 1. Reverse-coding suffix (`is_reverse_coded`, `_strip_reverse_suffix`)

A feature is considered **reverse coded** if its name ends in `_deviation` or
`_percentile_rank`:

```python
def is_reverse_coded(feature: str) -> bool:
    return feature.endswith(("_deviation", "_percentile_rank"))
```

Before classifying a feature into `movement` / `geometry` / `other`, this
suffix (and `_deviation_percentile_rank` combos) is stripped off first, so
the semantic group reflects what the feature *measures*, not whether it has
been reverse-coded:

```python
def _strip_reverse_suffix(feature: str) -> str:
    for suffix in ("_percentile_rank", "_deviation"):
        if feature.endswith(suffix):
            feature = feature[: -len(suffix)]
    return feature
```

- Example: `Mouse1_mean_euclid_distances_mean_10_percentile_rank` → strips to
  `Mouse1_mean_euclid_distances_mean_10` → classified as `geometry`.
- Example: `Total_movement_M1_mean_10_percentile_rank` → strips to
  `Total_movement_M1_mean_10` → classified as `movement`.

**Known edge case:** `Sum_probabilities_deviation_percentile_rank` strips
both suffixes down to `Sum_probabilities`, which doesn't match either
pattern, so it lands in `other`.

## 2. Semantic group (`assign_group`)

Applied to the *stripped* feature name (see above), in this order:

**`tortuosity`** — matches if the stripped name starts with `Tortuosity_`.
Split out into its own group per your answer (previously fell into `other`).

**`movement`** — matches if the stripped name starts with any of:
- `Movement_mouse_nose`, `Movement_mouse_tail_base`, `Movement_mouse_left_ear`,
  `Movement_mouse_right_ear`, `Movement_mouse_head_base`
- `Total_movement_all_bodyparts` (with or without a trailing `_M1`)
- `Total_movement_M1_`
- `Tail_base_movement_M1_`
- `Head_base_movement_M1_`
- `Nose_movement_M1_`

**`geometry`** — else matches if the stripped name starts with any of:
- `Mouse_nose_to_tail`, `Mouse_head_to_tail`, `Mouse_Ear_distance`
- `M1_` (any hull-distance feature, e.g. `M1_mean_euclidean_distance_hull`)
- `Mouse1_smallest_euclid_distances_`, `Mouse1_largest_euclid_distances_`,
  `Mouse1_mean_euclid_distances_`

**`other`** — anything that matches neither pattern.

Current counts (153 features total): `movement` = 74, `geometry` = 70,
`tortuosity` = 5, `other` = 4.

**Resolved:** `Total_movement_all_bodyparts_deviation` now correctly matches
the movement pattern (the `_M1` suffix requirement was made optional) and is
classified as `movement`, consistent with its non-reverse-coded sibling
`Total_movement_all_bodyparts_M1`.

**`other` group now contains only these 4 features (pending your decision):**
```
Sum_probabilities
Sum_probabilities_deviation
Sum_probabilities_percentile_rank
Sum_probabilities_deviation_percentile_rank
```
You didn't mark an answer for whether these should get their own
"detection confidence" group — let me know and I'll split them out the same
way `tortuosity` was split out.

**Note:** the `M1_*_euclidean_distance_hull` and `M1_*_euclid_distances_hull_deviation`
features are *not* in `other` — they're already correctly classified as
`geometry` via the `M1_` prefix rule. Only their `bodypart` label is `other`
(see section 3), since "hull" isn't one of the recognised body-part keywords.

## 3. Body part (`assign_bodypart`)

Applied to the lowercased full feature name (including any reverse-coding
suffix), first match wins:

1. Contains `"all_bodyparts"` → `"whole mouse"`
2. Contains one of (checked in this order): `"nose"`, `"tail_base"`,
   `"head_base"`, `"ear"` → corresponding label (`"nose"`, `"tail base"`,
   `"head base"`, `"ears"`)
3. Starts with `"mouse_"` / `"mouse1_"`, or contains `"total"` →
   `"whole mouse"`
4. Otherwise → `"other"`

**Resolved:** `"left ear"` and `"right ear"` are merged into a single
`"ears"` label, and `Mouse_Ear_distance` (previously `"whole mouse"`) is now
also classified as `"ears"`.

Current counts before combining: `whole mouse` = 81, `nose` = 18, `head base` = 17,
`tail base` = 16, `other` = 16, `ears` = 3, `all body parts` = 2. After combining,
`whole mouse` is expected to contain 83 features.

**`other` bodypart currently contains:** the 6 `M1_*_hull*` / `M1_*_hull_deviation`
distance features, the 5 `Tortuosity_Mouse1_*` features, and the
5 `Sum_probabilities*` features. None of these mention a body part by name,
so they fall through to `"other"`.

## 4. Time window (`get_time_window`)

Extracts a numeric suffix from the end of the feature name via
`re.search(r"_(\d+(?:\.\d+)?)(?:_deviation|_percentile_rank)*$", feature)`,
allowing the number to be followed by zero or more reverse-coding suffixes.
If no numeric suffix is found, the feature is labelled `"none"`.

**Resolved:** a reverse-coded feature like `..._mean_10_percentile_rank` now
resolves to window `"10"` instead of `"none"`, since the number only needs
to be at the end once any `_deviation`/`_percentile_rank` suffixes are
accounted for. This stays anchored to the end of the string (rather than
searching anywhere in the name), so it can't accidentally match the `1` in
`M1_`/`Mouse1_` prefixes.

Current counts: `none` = 23, and 26 features each for windows `2`, `5`, `6`,
`7.5`, `10`.

## 5. Display group (`display_group` column)

Simple concatenation for plotting: `f"{group}, reverse coded"` if
`reverse_coded` else just `group`. Current values: `movement`, `geometry`,
`geometry, reverse coded`, `movement, reverse coded`, `other`,
`other, reverse coded`.

---

## Please mark up / annotate

- [x] Should `Total_movement_all_bodyparts_deviation` be `movement` instead of `other`? — **yes, done.**
- [x] Should `Tortuosity_*` be its own group rather than `other`? — **yes, done** (new `tortuosity` group).
- [ ] Should `Sum_probabilities*` be its own group (detection confidence) rather than `other`? — **not yet answered**, still in `other`. Let me know and I'll split it out.
- [x] Should `get_time_window` run on the reverse-suffix-stripped name so reverse-coded features get a real window value instead of `"none"`? — **done**, using an end-anchored regex that allows the reverse-coding suffix after the number (safe against `M1_`/`Mouse1_` false matches).
- [x] Any bodypart labels you want split/merged (e.g. `Mouse_Ear_distance` → dedicated `"ears"` label)? — **done** (merged in an earlier chat turn).
- [ ] Any other feature prefixes you know of that aren't covered by the current regexes? — none noted.
