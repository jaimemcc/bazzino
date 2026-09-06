from pathlib import Path
import pandas as pd
import numpy as np
import re
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

DATAFOLDER = Path('data/shap_values')
shap_frames = []
raw_frames = []
for folder in sorted([f for f in DATAFOLDER.iterdir() if f.is_dir()]):
    shap_frames.append(pd.read_csv(folder / 'SHAP_values_Appetitive.csv'))
    raw_frames.append(pd.read_csv(folder / 'RAW_SHAP_feature_values_Appetitive.csv'))

df_shap = pd.concat(shap_frames, axis=0, ignore_index=True)
df_shap_raw = pd.concat(raw_frames, axis=0, ignore_index=True)
df_shap = df_shap.drop(columns=['Unnamed: 0','Prediction_probability','Sum','Expected_value'], errors='ignore')
df_shap_raw = df_shap_raw.drop(columns=['Unnamed: 0'], errors='ignore')

subset = df_shap.query('Appetitive == 1').copy()
raw_subset = df_shap_raw.loc[subset.index].copy()
feature_cols = [c for c in subset.columns if c not in {'Appetitive'}]


def assign_group(col):
    if re.search(r'^(Mouse_nose_to_tail|Mouse_head_to_tail|Mouse_Ear_distance)$', col):
        return 'body_geometry'
    if re.search(r'^Movement_mouse_(nose|tail_base|left_ear|right_ear|head_base)$', col):
        return 'body_part_movement'
    if re.search(r'^(Total_movement_all_bodyparts_M1|Total_movement_M1_)', col):
        return 'aggregate_movement'
    if re.search(r'^(M1_|Mouse1_(smallest|largest|mean)_euclid_distances_)', col):
        return 'hull_geometry'
    if re.search(r'^(Tail_base_movement_M1_|Head_base_movement_M1_|Nose_movement_M1_)', col):
        return 'rolling_window_summary'
    if re.search(r'(_deviation|_percentile_rank)$', col):
        return 'derived_scores'
    return 'other'

rows = []
for feat in feature_cols:
    if feat not in raw_subset.columns:
        continue
    x = pd.to_numeric(raw_subset[feat], errors='coerce')
    y = pd.to_numeric(subset[feat], errors='coerce')
    valid = ~(x.isna() | y.isna())
    x = x[valid]
    y = y[valid]
    if len(x) < 10:
        continue
    corr = np.corrcoef(x, y)[0, 1]
    q = pd.qcut(x, q=4, labels=['Q1','Q2','Q3','Q4'], duplicates='drop')
    qmeans = y.groupby(q).mean()
    trend = qmeans.iloc[-1] - qmeans.iloc[0] if len(qmeans) >= 2 else np.nan
    rows.append({'feature': feat, 'group': assign_group(feat), 'mean_abs_shap': abs(y).mean(), 'mean_shap': y.mean(), 'corr_raw_vs_shap': corr, 'trend': trend})
metrics = pd.DataFrame(rows)

summary_group = metrics.groupby('group').agg(
    n_features=('feature','count'),
    total_importance=('mean_abs_shap','sum'),
    mean_importance=('mean_abs_shap','mean'),
    mean_corr=('corr_raw_vs_shap','mean'),
    positive_corr_count=('corr_raw_vs_shap', lambda s: (s > 0).sum()),
    negative_corr_count=('corr_raw_vs_shap', lambda s: (s < 0).sum()),
).sort_values('total_importance', ascending=False)

print('GROUP SUMMARY')
print(summary_group.to_string())
print('\nTOP FEATURES BY IMPORTANCE')
print(metrics.sort_values('mean_abs_shap', ascending=False).head(20)[['feature','group','mean_abs_shap','mean_shap','corr_raw_vs_shap','trend']].to_string(index=False))
print('\nTOP FEATURES WITH POSITIVE RAW-SHAP CORR')
print(metrics.sort_values('corr_raw_vs_shap', ascending=False).head(20)[['feature','group','corr_raw_vs_shap','mean_abs_shap','trend']].to_string(index=False))
print('\nTOP FEATURES WITH NEGATIVE RAW-SHAP CORR')
print(metrics.sort_values('corr_raw_vs_shap', ascending=True).head(20)[['feature','group','corr_raw_vs_shap','mean_abs_shap','trend']].to_string(index=False))

X = subset[feature_cols].astype(float).fillna(0)
X_std = StandardScaler().fit_transform(X)
pca = PCA(n_components=2, random_state=42)
pca.fit(X_std)
loadings = pd.DataFrame(pca.components_.T, index=feature_cols, columns=['PC1','PC2'])
loadings['group'] = loadings.index.map(assign_group)
loadings['abs_PC1'] = loadings['PC1'].abs()
loadings['abs_PC2'] = loadings['PC2'].abs()
print('\nTOP PC1 LOADINGS')
print(loadings.sort_values('abs_PC1', ascending=False).head(20)[['group','PC1','abs_PC1']].to_string())
print('\nTOP PC2 LOADINGS')
print(loadings.sort_values('abs_PC2', ascending=False).head(20)[['group','PC2','abs_PC2']].to_string())
print('\nGROUP PC LOADINGS')
print(loadings.groupby('group').agg(mean_abs_pc1=('abs_PC1','mean'), mean_abs_pc2=('abs_PC2','mean')).sort_values('mean_abs_pc1', ascending=False).to_string())
