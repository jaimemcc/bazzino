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
## Trying to use low conc deplete and high conc replete to classify trials based on aversive-appetitive axis

# %%
from pathlib import Path
import sys
import os
os.environ["OMP_NUM_THREADS"] = "2"  # To avoide warning in clustering, can try removing

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import trompy as tp

import dill

sys.path.insert(0, str(Path("../src").resolve()))
from pickle_compat import enable_dill_pathlib_compat
enable_dill_pathlib_compat()

# %run "..//src//plots_for_spectral"

DATAFOLDER = Path("..//data")
RESULTSFOLDER = Path("..//results")

# %%
with open(DATAFOLDER / "bazzino_data_for_frejus.pickle", "rb") as f:
    data = dill.load(f)

x_array = data["x_array"]
snips_photo = data["snips_photo"]
snips_vel = data["snips_vel"]

# %%
pcafile = RESULTSFOLDER / "transformed_data_vel.pickle"

with open(pcafile, 'rb') as f:
    pca = dill.load(f)

pc1 = pca[:,0]
pca_data = pca[:, :3]

# %%
pca_data.shape

# %%
x_array.shape

# %%
# Create a label column based on condition + infusiontype
# 1 if (condition == 'deplete' and infusiontype == '10NaCl')
# 0 if (condition == 'replete' and infusiontype == '45NaCl')
# NaN otherwise
import numpy as np
import pandas as pd

cond_one = (x_array['condition'].eq('deplete')) & (x_array['infusiontype'].eq('10NaCl'))
cond_zero = (x_array['condition'].eq('replete')) & (x_array['infusiontype'].eq('45NaCl'))

labels = np.select([cond_one, cond_zero], [1, 0], default=np.nan)
# Use pandas nullable integer to preserve NaNs as <NA>
x_array['class_label'] = pd.Series(labels, index=x_array.index).astype('Int64')

print("Label distribution (including NaNs):")
print(x_array['class_label'].value_counts(dropna=False))

# Optional: peek at first few rows to confirm
x_array[['condition', 'infusiontype', 'class_label']].head()

# %%
# Classify labeled trials from 3 PCs with SVM and 10-fold CV
# Using newly created 'class_label' (1=deplete/10NaCl, 0=replete/45NaCl)
label_col = 'class_label'

if label_col not in x_array.columns:
    raise ValueError(f"Label column '{label_col}' not found in x_array. Available columns: {list(x_array.columns)}")

# Filter to rows with non-null labels
mask = x_array[label_col].notna()
X = pca_data[mask]
y = x_array.loc[mask, label_col].astype(int).to_numpy()

# Basic sanity checks
assert len(X) == len(y), f"Length mismatch: X={len(X)}, y={len(y)}"
assert set(np.unique(y)).issubset({0, 1}), f"y must be binary 0/1, got: {np.unique(y)}"

from sklearn.model_selection import StratifiedKFold, cross_validate
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

pipe = Pipeline([
    ("scaler", StandardScaler()),
    ("svm", SVC(kernel="rbf", probability=True, class_weight="balanced"))
])

skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
scoring = {
    "accuracy": "accuracy",
    "roc_auc": "roc_auc"
}

cv_res = cross_validate(pipe, X, y, cv=skf, scoring=scoring, return_estimator=False, n_jobs=-1)

print({
    "cv_accuracy_mean": float(np.mean(cv_res["test_accuracy"])),
    "cv_accuracy_std": float(np.std(cv_res["test_accuracy"])),
    "cv_roc_auc_mean": float(np.mean(cv_res["test_roc_auc"])),
    "cv_roc_auc_std": float(np.std(cv_res["test_roc_auc"]))
})

# %%
# Train/test split evaluation and ROC/Confusion Matrix on held-out set
from sklearn.model_selection import train_test_split
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score, f1_score,
    roc_auc_score, RocCurveDisplay, ConfusionMatrixDisplay
)

X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, stratify=y, random_state=42
)

pipe.fit(X_train, y_train)

y_pred = pipe.predict(X_test)
y_prob = pipe.predict_proba(X_test)[:, 1]

metrics = {
    "test_accuracy": float(accuracy_score(y_test, y_pred)),
    "test_precision": float(precision_score(y_test, y_pred, zero_division=0)),
    "test_recall": float(recall_score(y_test, y_pred, zero_division=0)),
    "test_f1": float(f1_score(y_test, y_pred, zero_division=0)),
    "test_roc_auc": float(roc_auc_score(y_test, y_prob))
}
print(metrics)

# Plots
fig, ax = plt.subplots(1, 2, figsize=(10,4))
RocCurveDisplay.from_predictions(y_test, y_prob, name="SVM", ax=ax[0])
ax[0].set_title("ROC curve (held-out)")

ConfusionMatrixDisplay.from_predictions(y_test, y_pred, ax=ax[1], colorbar=False)
ax[1].set_title("Confusion matrix (held-out)")
sns.despine(ax=ax[0])
plt.tight_layout()
plt.show()

# %%
# Persist trained pipeline for future unseen data evaluation (optional)
import dill
model_path = RESULTSFOLDER / "svm_trial_classifier_vel.pkl"
try:
    with open(model_path, "wb") as f:
        dill.dump(pipe, f)
    print(f"Saved trained model to {model_path}")
except Exception as e:
    print(f"Could not save model: {e}")

# %%

# %%
# Benchmark multiple classifiers against SVM with same 10-fold CV
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier, ExtraTreesClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.svm import SVC, LinearSVC
from sklearn.calibration import CalibratedClassifierCV
from sklearn.model_selection import StratifiedKFold, cross_validate
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
import pandas as pd

label_col = 'class_label'
mask = x_array[label_col].notna()
X = pca_data[mask]
y = x_array.loc[mask, label_col].astype(int).to_numpy()

skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
scoring = {"accuracy": "accuracy", "roc_auc": "roc_auc"}

models = {
    "SVM_RBF": SVC(kernel="rbf", probability=True, class_weight="balanced"),
    "LogReg_L2": LogisticRegression(max_iter=200, class_weight="balanced"),
    "LinearSVM": CalibratedClassifierCV(LinearSVC(class_weight="balanced", max_iter=5000)),
    "RandomForest": RandomForestClassifier(n_estimators=300, class_weight="balanced", random_state=42),
    "ExtraTrees": ExtraTreesClassifier(n_estimators=300, random_state=42),
    "GradientBoosting": GradientBoostingClassifier(random_state=42),
    "KNN_7": KNeighborsClassifier(n_neighbors=7),
    "LDA": LinearDiscriminantAnalysis()
}

rows = []
for name, clf in models.items():
    pipe = Pipeline([
        ("scaler", StandardScaler()),
        ("clf", clf)
    ])
    cv_res = cross_validate(pipe, X, y, cv=skf, scoring=scoring, n_jobs=-1)
    rows.append({
        "model": name,
        "acc_mean": cv_res['test_accuracy'].mean(),
        "acc_std": cv_res['test_accuracy'].std(),
        "auc_mean": cv_res['test_roc_auc'].mean(),
        "auc_std": cv_res['test_roc_auc'].std()
    })

summary_df = pd.DataFrame(rows).sort_values("auc_mean", ascending=False)
print(summary_df)
summary_df

# %%

# %%
# Threshold tuning for SVM pipeline (optimize F1 and Youden's J)
from sklearn.metrics import f1_score, roc_curve, precision_recall_curve
import numpy as np
import matplotlib.pyplot as plt

label_col = 'class_label'
mask = x_array[label_col].notna()
X_full = pca_data[mask]
y_full = x_array.loc[mask, label_col].astype(int).to_numpy()

# Reuse train/test split logic (same random_state for reproducibility)
from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(
    X_full, y_full, test_size=0.2, stratify=y_full, random_state=42
)

# Fit (pipe assumed defined earlier; if not, redefine quickly)
try:
    pipe
except NameError:
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.svm import SVC
    pipe = Pipeline([
        ("scaler", StandardScaler()),
        ("svm", SVC(kernel="rbf", probability=True, class_weight="balanced"))
    ])

pipe.fit(X_train, y_train)
probs = pipe.predict_proba(X_test)[:,1]

thresholds = np.linspace(0.0, 1.0, 201)
records = []
for t in thresholds:
    y_pred_t = (probs >= t).astype(int)
    tp = np.sum((y_test == 1) & (y_pred_t == 1))
    fn = np.sum((y_test == 1) & (y_pred_t == 0))
    fp = np.sum((y_test == 0) & (y_pred_t == 1))
    tn = np.sum((y_test == 0) & (y_pred_t == 0))
    # Avoid divide by zero
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0  # TPR / Recall
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0.0  # TNR
    f1 = f1_score(y_test, y_pred_t, zero_division=0)
    youden_j = sensitivity + specificity - 1
    records.append((t, f1, sensitivity, specificity, youden_j))

rec_arr = np.array(records)
# Best thresholds
best_f1_idx = np.argmax(rec_arr[:,1])
best_j_idx = np.argmax(rec_arr[:,4])

best_f1_thr = rec_arr[best_f1_idx,0]
best_j_thr = rec_arr[best_j_idx,0]

print({
    "best_f1_threshold": float(best_f1_thr),
    "best_f1": float(rec_arr[best_f1_idx,1]),
    "best_j_threshold": float(best_j_thr),
    "best_youden_j": float(rec_arr[best_j_idx,4])
})

# Plot F1 and Youden's J vs threshold
fig, ax = plt.subplots(1,2, figsize=(10,4))
ax[0].plot(rec_arr[:,0], rec_arr[:,1], label='F1')
ax[0].axvline(best_f1_thr, color='r', linestyle='--', label=f'Best F1 t={best_f1_thr:.2f}')
ax[0].set_xlabel('Threshold')
ax[0].set_ylabel('F1 score')
ax[0].legend()

ax[1].plot(rec_arr[:,0], rec_arr[:,4], label="Youden's J", color='g')
ax[1].axvline(best_j_thr, color='r', linestyle='--', label=f'Best J t={best_j_thr:.2f}')
ax[1].set_xlabel('Threshold')
ax[1].set_ylabel("Youden's J")
ax[1].legend()
plt.tight_layout()
plt.show()

# Provide a function to apply chosen threshold later
chosen_threshold = float(best_j_thr)

def predict_with_threshold(X_new, threshold=chosen_threshold):
    """Return binary predictions using tuned threshold on probability outputs."""
    probs_new = pipe.predict_proba(X_new)[:,1]
    return (probs_new >= threshold).astype(int)

print(f"Stored chosen_threshold={chosen_threshold:.3f} (Youden's J criterion)")

# %%

# %%
# Save chosen threshold alongside the trained pipeline
import json
from pathlib import Path

# Ensure we have a model trained in this session; if not, fit quickly on full data
label_col = 'class_label'
mask = x_array[label_col].notna()
X_full = pca_data[mask]
y_full = x_array.loc[mask, label_col].astype(int).to_numpy()

try:
    pipe
except NameError:
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.svm import SVC
    pipe = Pipeline([
        ("scaler", StandardScaler()),
        ("svm", SVC(kernel="rbf", probability=True, class_weight="balanced"))
    ])
    pipe.fit(X_full, y_full)

# Use chosen_threshold if defined; else default to 0.5
try:
    chosen_threshold
except NameError:
    chosen_threshold = 0.5

# Save separate JSON for easy reading
threshold_info = {
    "label_col": label_col,
    "chosen_threshold": float(chosen_threshold),
    "criterion": "youden_j",
}
json_path = RESULTSFOLDER / "svm_trial_classifier_threshold.json"
with open(json_path, "w") as f:
    json.dump(threshold_info, f, indent=2)
print(f"Saved threshold JSON to {json_path}")

# Save combined pickle with pipeline and threshold
bundle = {"pipeline": pipe, "threshold": float(chosen_threshold), "label_col": label_col}
model_bundle_path = RESULTSFOLDER / "svm_trial_classifier_with_threshold.pkl"
with open(model_bundle_path, "wb") as f:
    dill.dump(bundle, f)
print(f"Saved model bundle to {model_bundle_path}")

# %%

# %%
# Classify unseen-style subset (deplete & 45NaCl) and plot vs trial
import dill, json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

subset_query = "condition == 'deplete' & infusiontype == '45NaCl'"
subset_df = x_array.query(subset_query).copy()

if subset_df.empty:
    raise ValueError(f"Query returned no rows: {subset_query}")

# Identify indices into original x_array to slice pca_data consistently
subset_idx = subset_df.index.to_numpy()
X_subset = pca_data[subset_idx]

# Load model bundle if available
bundle_path = RESULTSFOLDER / "svm_trial_classifier_with_threshold.pkl"
if bundle_path.exists():
    with open(bundle_path, 'rb') as f:
        bundle = dill.load(f)
    pipe_loaded = bundle['pipeline']
    threshold_loaded = bundle.get('threshold', 0.5)
else:
    print("Bundle not found; refitting pipeline quickly on labeled data.")
    label_col = 'class_label'
    mask = x_array[label_col].notna()
    X_full = pca_data[mask]
    y_full = x_array.loc[mask, label_col].astype(int).to_numpy()
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.svm import SVC
    pipe_loaded = Pipeline([
        ("scaler", StandardScaler()),
        ("svm", SVC(kernel="rbf", probability=True, class_weight="balanced"))
    ])
    pipe_loaded.fit(X_full, y_full)
    threshold_loaded = 0.5

probs_subset = pipe_loaded.predict_proba(X_subset)[:,1]
preds_subset = (probs_subset >= threshold_loaded).astype(int)

subset_df = subset_df.assign(pred_prob=probs_subset, pred_class=preds_subset)

# Plot trial vs probability with threshold line
fig, ax = plt.subplots(figsize=(7,4))
ax.plot(subset_df.trial, subset_df.pred_prob, marker='o', linestyle='-', color='#016895', label='Predicted P(class=1)')
ax.axhline(threshold_loaded, color='#F4795B', linestyle='--', label=f'Threshold={threshold_loaded:.2f}')

# Overlay predicted classes as filled / hollow markers
for _, row in subset_df.iterrows():
    ax.scatter(row.trial, row.pred_prob, s=60,
               facecolors=('#F4795B' if row.pred_class==1 else 'none'),
               edgecolors='#F4795B', linewidth=1.2,
               alpha=0.9)

ax.set_xlabel('Trial')
ax.set_ylabel('Predicted probability (class=1)')
ax.set_title('Classification on subset: deplete & 45NaCl')
ax.legend(loc='best', frameon=False)
sns.despine(ax=ax)
plt.tight_layout()
plt.show()

# Display first few classified rows
subset_df[['trial','condition','infusiontype','pred_prob','pred_class']].head()

# %%

# %%
# Faceted per-rat plots of predicted probability over trials for subset
import seaborn as sns
import matplotlib.pyplot as plt

# Determine rat identifier column
rat_col = None
for candidate in ['id', 'rat', 'animal']:
    if candidate in x_array.columns:
        rat_col = candidate
        break
if rat_col is None:
    raise ValueError("Could not find a rat identifier column among ['id','rat','animal'] in x_array.columns")

subset_query = "condition == 'deplete' & infusiontype == '45NaCl'"
subset_df = x_array.query(subset_query).copy()
subset_idx = subset_df.index.to_numpy()
X_subset = pca_data[subset_idx]

# Load existing predictions if available, else compute
if 'pred_prob' not in subset_df.columns or 'pred_class' not in subset_df.columns:
    bundle_path = RESULTSFOLDER / "svm_trial_classifier_with_threshold.pkl"
    if bundle_path.exists():
        import dill
        with open(bundle_path, 'rb') as f:
            bundle = dill.load(f)
        pipe_loaded = bundle['pipeline']
        threshold_loaded = bundle.get('threshold', 0.5)
    else:
        label_col = 'class_label'
        mask = x_array[label_col].notna()
        X_full = pca_data[mask]
        y_full = x_array.loc[mask, label_col].astype(int).to_numpy()
        from sklearn.pipeline import Pipeline
        from sklearn.preprocessing import StandardScaler
        from sklearn.svm import SVC
        pipe_loaded = Pipeline([
            ("scaler", StandardScaler()),
            ("svm", SVC(kernel="rbf", probability=True, class_weight="balanced"))
        ])
        pipe_loaded.fit(X_full, y_full)
        threshold_loaded = 0.5
    probs_subset = pipe_loaded.predict_proba(X_subset)[:,1]
    preds_subset = (probs_subset >= threshold_loaded).astype(int)
    subset_df = subset_df.assign(pred_prob=probs_subset, pred_class=preds_subset)
else:
    # Need threshold for reference line if not stored; fallback
    threshold_loaded = 0.5

# Sort within each rat by trial for cleaner lines
subset_df = subset_df.sort_values(['trial'])

# Facet settings
n_rats = subset_df[rat_col].nunique()
col_wrap = 4 if n_rats > 4 else n_rats

# Build FacetGrid
g = sns.FacetGrid(subset_df, col=rat_col, col_wrap=col_wrap, sharey=True, height=2.4)

# Line plot function
def _plot(data, color):
    ax = plt.gca()
    ax.plot(data['trial'], data['pred_prob'], marker='o', linestyle='-', color='#016895', linewidth=1)
    ax.axhline(threshold_loaded, color='#F4795B', linestyle='--', linewidth=0.8)
    # Overlay predicted class markers
    ax.scatter(data['trial'], data['pred_prob'], s=35,
               facecolors=['#F4795B' if pc==1 else 'none' for pc in data['pred_class']],
               edgecolors='#F4795B', linewidth=0.8)
    ax.set_xlabel('Trial')
    ax.set_ylabel('P(class=1)')
    ax.set_ylim(-0.02, 1.02)
    ax.tick_params(axis='x', labelsize=8)
    ax.tick_params(axis='y', labelsize=8)
    sns.despine(ax=ax)

g.map_dataframe(_plot)

# Add a super title
plt.subplots_adjust(top=0.88)
plt.suptitle("Per-rat predicted probabilities (deplete & 45NaCl)\nFilled markers = class 1 above threshold", fontsize=12)
plt.show()

# Summary counts per rat
rat_summary = subset_df.groupby(rat_col)['pred_class'].agg(['count','sum'])
rat_summary = rat_summary.rename(columns={'count':'n_trials','sum':'n_pred_class1'})
print(rat_summary)

# %%

# %%
# Fit sigmoids per rat on predicted probabilities vs trial and facet-plot (robust init + normalized x)
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Helper: logistic without baseline (probabilities should be ~[0,1])
def logistic(x, L, x0, k):
    return L / (1.0 + np.exp(-k * (x - x0)))

# Identify rat id column
rat_col = None
for candidate in ['id', 'rat', 'animal']:
    if candidate in x_array.columns:
        rat_col = candidate
        break
if rat_col is None:
    raise ValueError("Could not find a rat identifier column among ['id','rat','animal'] in x_array.columns")

# Subset and ensure predictions are present
subset_query = "condition == 'deplete' & infusiontype == '45NaCl'"
subset_df = x_array.query(subset_query).copy()
subset_idx = subset_df.index.to_numpy()
X_subset = pca_data[subset_idx]

if 'pred_prob' not in subset_df.columns:
    # Load model bundle or refit
    bundle_path = RESULTSFOLDER / "svm_trial_classifier_with_threshold.pkl"
    if bundle_path.exists():
        import dill
        with open(bundle_path, 'rb') as f:
            bundle = dill.load(f)
        pipe_loaded = bundle['pipeline']
    else:
        label_col = 'class_label'
        mask = x_array[label_col].notna()
        X_full = pca_data[mask]
        y_full = x_array.loc[mask, label_col].astype(int).to_numpy()
        from sklearn.pipeline import Pipeline
        from sklearn.preprocessing import StandardScaler
        from sklearn.svm import SVC
        pipe_loaded = Pipeline([
            ("scaler", StandardScaler()),
            ("svm", SVC(kernel="rbf", probability=True, class_weight="balanced"))
        ])
        pipe_loaded.fit(X_full, y_full)
    subset_df['pred_prob'] = pipe_loaded.predict_proba(X_subset)[:,1]

# Fit per rat
fit_rows = []
subset_df = subset_df.sort_values(['trial'])
subset_df['y_fit'] = np.nan

for r, df_r in subset_df.groupby(rat_col):
    x = df_r['trial'].to_numpy(dtype=float)
    y = df_r['pred_prob'].to_numpy(dtype=float)

    # Guard conditions
    if len(np.unique(x)) < 3:
        fit_rows.append({rat_col: r, 'L': np.nan, 'x0': np.nan, 'k': np.nan, 'x0_orig': np.nan, 'r2': np.nan, 'n': len(x)})
        continue

    # Normalize x for stabler optimization
    x_mean = float(np.mean(x))
    x_std = float(np.std(x))
    if not np.isfinite(x_std) or x_std == 0:
        x_std = 1.0
    x_norm = (x - x_mean) / x_std

    # Slightly clip y to avoid flat 0/1 plateaus
    y_clipped = np.clip(y, 1e-4, 1-1e-4)

    # Initial guesses (multi-start on slope magnitude)
    L_init = float(np.clip(np.quantile(y_clipped, 0.95), 0.6, 1.2))
    x0_init = 0.0  # in normalized units ~ midpoint
    # Correlation-based sign guess for k
    try:
        corr = float(np.corrcoef(x, y)[0,1])
    except Exception:
        corr = 0.0
    if not np.isfinite(corr):
        corr = 0.0
    sign = 1.0 if corr >= 0 else -1.0
    k_cands = [0.5, 1.0, 2.0]

    best = None
    best_err = np.inf
    for km in k_cands:
        p0 = [L_init, x0_init, sign*km]
        bounds = ([0.1, -3.0, -10.0], [2.0, 3.0, 10.0])
        try:
            popt, _ = curve_fit(logistic, x_norm, y_clipped, p0=p0, bounds=bounds, maxfev=50000)
            y_hat = logistic(x_norm, *popt)
            rss = float(np.sum((y_clipped - y_hat)**2))
            if rss < best_err:
                best_err = rss
                best = popt, y_hat
        except Exception:
            continue

    if best is None:
        fit_rows.append({rat_col: r, 'L': np.nan, 'x0': np.nan, 'k': np.nan, 'x0_orig': np.nan, 'r2': np.nan, 'n': len(x)})
        continue

    popt, y_hat = best
    L_hat, x0_hat_norm, k_hat = map(float, popt)
    # Convert x0 back to original trial units
    x0_hat_orig = x0_hat_norm * x_std + x_mean

    ss_res = float(np.sum((y - y_hat) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2)) if np.var(y) > 0 else np.nan
    r2 = float(1.0 - ss_res / ss_tot) if np.isfinite(ss_tot) and ss_tot > 0 else np.nan

    fit_rows.append({rat_col: r, 'L': L_hat, 'x0': x0_hat_norm, 'k': k_hat, 'x0_orig': x0_hat_orig, 'r2': r2, 'n': len(x)})
    # Store fitted on observed x for plotting (using normalized x)
    subset_df.loc[df_r.index, 'y_fit'] = y_hat

fits_df = pd.DataFrame(fit_rows).sort_values('x0_orig')
print("Sigmoid fits per rat (sorted by x0 in original trial units):")
print(fits_df)

# Optional: save fits
try:
    csv_path = RESULTSFOLDER / "sigmoid_fits_deplete_45NaCl_per_rat.csv"
    fits_df.to_csv(csv_path, index=False)
    print(f"Saved fits to {csv_path}")
except Exception as e:
    print(f"Could not save fits CSV: {e}")

# Facet plot: points + fitted curve
n_rats = subset_df[rat_col].nunique()
col_wrap = 4 if n_rats > 4 else n_rats

g = sns.FacetGrid(subset_df, col=rat_col, col_wrap=col_wrap, sharey=True, height=2.4)

def _plot_with_fit(data, color):
    ax = plt.gca()
    # Points
    ax.scatter(data['trial'], data['pred_prob'], s=28, color='#016895', alpha=0.8, label='Pred prob')
    # Fitted line (sorted by trial)
    d = data.sort_values('trial')
    if d['y_fit'].notna().any():
        ax.plot(d['trial'], d['y_fit'], color='#F4795B', linewidth=1.5, label='Sigmoid fit')
    ax.set_xlabel('Trial')
    ax.set_ylabel('P(class=1)')
    ax.set_ylim(-0.02, 1.02)
    ax.tick_params(axis='x', labelsize=8)
    ax.tick_params(axis='y', labelsize=8)
    sns.despine(ax=ax)

g.map_dataframe(_plot_with_fit)
plt.subplots_adjust(top=0.86)
plt.suptitle('Per-rat sigmoid fits on predicted probabilities (deplete & 45NaCl) — normalized x, multistart', fontsize=12)
plt.show()

fits_df

# %%

# %%
# Multi-model fits per rat: compare logistic variants, Gompertz, exponential, piecewise linear
# Adjusted: normalize x, multi-start slopes, tighter bounds to improve convergence
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Models -------------------------------------------------------------

def logistic3(x, L, x0, k):
    return L / (1 + np.exp(-k * (x - x0)))

def logistic4(x, A, L, x0, k):
    return A + (L - A) / (1 + np.exp(-k * (x - x0)))

def gompertz(x, a, b, c):
    return a * np.exp(-b * np.exp(-c * x))

def exp_sat(x, L, k):
    x = np.asarray(x)
    return L * (1 - np.exp(-k * (x - x.min())))

# Piecewise linear evaluation

def piecewise_linear_eval(x, m1, m2, cp, b1):
    x = np.asarray(x)
    y_cp = m1 * cp + b1
    return np.where(x <= cp, m1 * x + b1, y_cp + m2 * (x - cp))

# Utility -----------------------------------------------------------

def safe_curve_fit(func, x, y, p0, bounds=None, maxfev=60000):
    try:
        popt, _ = curve_fit(func, x, y, p0=p0, bounds=bounds if bounds else (-np.inf, np.inf), maxfev=maxfev)
        y_hat = func(x, *popt)
        return popt, y_hat, None
    except Exception as e:
        return None, None, e


def fit_piecewise_linear(x, y, n_grid=14):
    x_sorted = np.sort(np.unique(x))
    if len(x_sorted) < 5:
        return None, None, 'insufficient unique x'
    cps = np.linspace(x_sorted[1], x_sorted[-2], min(n_grid, len(x_sorted)-2))
    best = None
    best_rss = np.inf
    for cp in cps:
        left_mask = x <= cp
        right_mask = x > cp
        if left_mask.sum() < 3 or right_mask.sum() < 3:
            continue
        xl, yl = x[left_mask], y[left_mask]
        xr, yr = x[right_mask], y[right_mask]
        m1, b1 = np.polyfit(xl, yl, 1)
        m2, _ = np.polyfit(xr, yr, 1)
        y_hat = piecewise_linear_eval(x, m1, m2, cp, b1)
        rss = np.sum((y - y_hat)**2)
        if rss < best_rss:
            best_rss = rss
            best = (m1, m2, cp, b1, y_hat)
    if best is None:
        return None, None, 'fit failed'
    m1, m2, cp, b1, y_hat = best
    return (m1, m2, cp, b1), y_hat, None


def model_stats(x, y, y_hat, n_params):
    rss = np.sum((y - y_hat)**2)
    n = len(y)
    if n < n_params + 1:
        return {'rss': rss, 'aic': np.nan, 'bic': np.nan, 'r2': np.nan}
    var = np.var(y)
    r2 = 1 - rss / np.sum((y - np.mean(y))**2) if var > 0 else np.nan
    aic = n * np.log(max(rss / n, 1e-12)) + 2 * n_params
    bic = n * np.log(max(rss / n, 1e-12)) + np.log(n) * n_params
    return {'rss': rss, 'aic': aic, 'bic': bic, 'r2': r2}

# Data prep ---------------------------------------------------------
rat_col = None
for candidate in ['id', 'rat', 'animal']:
    if candidate in x_array.columns:
        rat_col = candidate
        break
if rat_col is None:
    raise ValueError("No rat identifier column found.")

subset_query = "condition == 'deplete' & infusiontype == '45NaCl'"
subset_df = x_array.query(subset_query).copy()
subset_idx = subset_df.index.to_numpy()
X_subset = pca_data[subset_idx]

if 'pred_prob' not in subset_df.columns:
    bundle_path = RESULTSFOLDER / "svm_trial_classifier_with_threshold.pkl"
    if bundle_path.exists():
        import dill
        with open(bundle_path, 'rb') as f:
            bundle = dill.load(f)
        pipe_loaded = bundle['pipeline']
    else:
        label_col = 'class_label'
        mask = x_array[label_col].notna()
        X_full = pca_data[mask]
        y_full = x_array.loc[mask, label_col].astype(int).to_numpy()
        from sklearn.pipeline import Pipeline
        from sklearn.preprocessing import StandardScaler
        from sklearn.svm import SVC
        pipe_loaded = Pipeline([
            ("scaler", StandardScaler()),
            ("svm", SVC(kernel="rbf", probability=True, class_weight="balanced"))
        ])
        pipe_loaded.fit(X_full, y_full)
    subset_df['pred_prob'] = pipe_loaded.predict_proba(X_subset)[:,1]

subset_df = subset_df.sort_values('trial')

# Fit all models per rat --------------------------------------------
results = []
subset_df['best_fit'] = np.nan
subset_df['best_model'] = np.nan

for r, df_r in subset_df.groupby(rat_col):
    x_raw = df_r['trial'].to_numpy(dtype=float)
    y = df_r['pred_prob'].to_numpy(dtype=float)
    if len(np.unique(x_raw)) < 6 or np.var(y) == 0:
        results.append({rat_col: r, 'model': 'NA', 'params': {}, 'rss': np.nan, 'aic': np.nan, 'bic': np.nan, 'r2': np.nan, 'n': len(x_raw)})
        continue

    # Normalize x
    x_mean = float(np.mean(x_raw))
    x_std = float(np.std(x_raw))
    if not np.isfinite(x_std) or x_std == 0:
        x_std = 1.0
    x = (x_raw - x_mean) / x_std

    y_clip = np.clip(y, 1e-4, 1-1e-4)

    model_fits = []

    # Inits
    L_est = float(np.clip(np.quantile(y_clip, 0.95), 0.6, 1.2))
    A_est = float(np.clip(np.quantile(y_clip, 0.05), 0.0, 0.4))
    x0_est = 0.0
    corr = np.corrcoef(x, y_clip)[0,1] if np.isfinite(np.corrcoef(x, y_clip)[0,1]) else 0
    sgn = 1 if corr >= 0 else -1
    k_vals = [0.5, 1.0, 2.0]

    # Logistic3 multi-start
    best_local = None
    best_rss = np.inf
    for kv in k_vals:
        p0 = [L_est, x0_est, sgn*kv]
        popt, yhat, err = safe_curve_fit(logistic3, x, y_clip, p0=p0, bounds=([0.2, -3, -10],[1.5, 3, 10]))
        if popt is not None:
            rss = np.sum((y_clip - yhat)**2)
            if rss < best_rss:
                best_rss = rss
                best_local = popt, yhat
    if best_local is not None:
        popt, yhat = best_local
        stats = model_stats(x, y_clip, yhat, 3)
        model_fits.append(('logistic3', popt, yhat, stats))

    # Logistic4
    best_local = None
    best_rss = np.inf
    for kv in k_vals:
        p0 = [A_est, L_est, x0_est, sgn*kv]
        popt, yhat, err = safe_curve_fit(logistic4, x, y_clip, p0=p0, bounds=([0.0, 0.2, -3, -10],[0.8, 1.5, 3, 10]))
        if popt is not None:
            rss = np.sum((y_clip - yhat)**2)
            if rss < best_rss:
                best_rss = rss
                best_local = popt, yhat
    if best_local is not None:
        popt, yhat = best_local
        stats = model_stats(x, y_clip, yhat, 4)
        model_fits.append(('logistic4', popt, yhat, stats))

    # Gompertz (single attempt)
    popt, yhat, err = safe_curve_fit(gompertz, x, y_clip, p0=[np.nanmax(y_clip), 1.0, 0.3], bounds=([0.2, 0.001, 0.001],[1.5, 5.0, 5.0]))
    if popt is not None:
        stats = model_stats(x, y_clip, yhat, 3)
        model_fits.append(('gompertz', popt, yhat, stats))

    # Exponential saturating
    popt, yhat, err = safe_curve_fit(exp_sat, x, y_clip, p0=[L_est, 0.2], bounds=([0.2, 0.0005],[1.5, 5.0]))
    if popt is not None:
        stats = model_stats(x, y_clip, yhat, 2)
        model_fits.append(('exp_sat', popt, yhat, stats))

    # Piecewise linear
    pw_params, yhat, err = fit_piecewise_linear(x, y_clip, n_grid=14)
    if pw_params is not None:
        stats = model_stats(x, y_clip, yhat, 4)
        model_fits.append(('piecewise_linear', pw_params, yhat, stats))

    if not model_fits:
        results.append({rat_col: r, 'model': 'none_fit', 'params': {}, 'rss': np.nan, 'aic': np.nan, 'bic': np.nan, 'r2': np.nan, 'n': len(x_raw)})
        continue

    # Select best by AIC
    best = min(model_fits, key=lambda mf: mf[3]['aic'] if np.isfinite(mf[3]['aic']) else np.inf)
    best_name, best_params, best_yhat, best_stats = best

    # Record all models (transform x0 back if present)
    for name, params, yhat_model, stats_model in model_fits:
        param_dict = {}
        if name.startswith('logistic'):
            # logistic3: L,x0,k ; logistic4: A,L,x0,k
            if name == 'logistic3':
                Lp, x0p, kp = params
                param_dict = {'L': float(Lp), 'x0_norm': float(x0p), 'k': float(kp), 'x0_orig': float(x0p * x_std + x_mean)}
            else:
                Ap, Lp, x0p, kp = params
                param_dict = {'A': float(Ap), 'L': float(Lp), 'x0_norm': float(x0p), 'k': float(kp), 'x0_orig': float(x0p * x_std + x_mean)}
        elif name == 'gompertz':
            ap, bp, cp = params
            param_dict = {'a': float(ap), 'b': float(bp), 'c': float(cp)}
        elif name == 'exp_sat':
            Lp, kp = params
            param_dict = {'L': float(Lp), 'k': float(kp)}
        elif name == 'piecewise_linear':
            m1, m2, cp, b1 = params
            # cp in normalized units; convert back
            cp_orig = cp * x_std + x_mean
            param_dict = {'m1': float(m1), 'm2': float(m2), 'cp_norm': float(cp), 'cp_orig': float(cp_orig), 'b1': float(b1)}

        results.append({
            rat_col: r,
            'model': name,
            'params': param_dict,
            'rss': stats_model['rss'],
            'aic': stats_model['aic'],
            'bic': stats_model['bic'],
            'r2': stats_model['r2'],
            'n': len(x_raw),
            'is_best': name == best_name
        })

    # Store best fit in subset_df for plotting (on original x order)
    subset_df.loc[df_r.index, 'best_fit'] = best_yhat
    subset_df.loc[df_r.index, 'best_model'] = best_name

# Summaries ---------------------------------------------------------
results_df = pd.DataFrame(results)
print("Model comparison (per rat * per model):")
print(results_df.head())

best_summary = results_df[results_df['is_best'] == True].copy()
print("Best model per rat:")
print(best_summary[[rat_col, 'model', 'aic', 'bic', 'r2', 'n']])

# Save summaries
try:
    best_summary_path = RESULTSFOLDER / 'best_model_fits_deplete_45NaCl_per_rat.csv'
    results_full_path = RESULTSFOLDER / 'all_model_fits_deplete_45NaCl_per_rat.csv'
    best_summary.to_csv(best_summary_path, index=False)
    results_df.to_csv(results_full_path, index=False)
    print(f"Saved best summary to {best_summary_path}\nSaved full results to {results_full_path}")
except Exception as e:
    print(f"Could not save CSVs: {e}")

# Facet plot of best model overlay -------------------------------------------------
n_rats = subset_df[rat_col].nunique()
col_wrap = 4 if n_rats > 4 else n_rats

subset_df_plot = subset_df.copy()
subset_df_plot['trial'] = subset_df_plot['trial'].astype(float)

# Plot
import matplotlib.pyplot as plt
import seaborn as sns

g = sns.FacetGrid(subset_df_plot, col=rat_col, col_wrap=col_wrap, sharey=True, height=2.4)

def _plot_best(data, color):
    ax = plt.gca()
    d = data.sort_values('trial')
    ax.scatter(d['trial'], d['pred_prob'], s=25, color='#016895', alpha=0.75)
    if d['best_fit'].notna().any():
        ax.plot(d['trial'], d['best_fit'], color='#C74632', linewidth=1.4)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel('Trial')
    ax.set_ylabel('P(class=1)')
    ax.set_title(str(d[rat_col].iloc[0]))
    ax.tick_params(axis='x', labelsize=8)
    ax.tick_params(axis='y', labelsize=8)
    sns.despine(ax=ax)

g.map_dataframe(_plot_best)
plt.subplots_adjust(top=0.87)
plt.suptitle('Best model fits per rat (normalized x, multistart)', fontsize=12)
plt.show()

best_summary

# %%

# %%
# Additional model possibilities: Richards (5PL), Weibull CDF, Probit CDF
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import curve_fit
from scipy.stats import norm

# ------------------------------
# Model definitions
# ------------------------------

def logistic3(x, L, x0, k):
    return L / (1 + np.exp(-k * (x - x0)))


def logistic4(x, A, L, x0, k):
    return A + (L - A) / (1 + np.exp(-k * (x - x0)))


def richards5(x, A, L, x0, k, v):
    # Generalized logistic (Richards/5PL) with asymmetry v
    # Reduces to logistic4 when v -> 1
    return A + (L - A) / (1 + np.exp(-k * (x - x0))) ** (1.0 / np.clip(v, 1e-6, None))


def gompertz(x, a, b, c):
    return a * np.exp(-b * np.exp(-c * x))


def exp_sat(x, L, k):
    # Saturating exponential rising toward L from baseline ~0 at min(x)
    x = np.asarray(x, dtype=float)
    return L * (1 - np.exp(-k * (x - float(np.min(x)))))


def weibull_cdf(x, A, L, x0, lam, k):
    # Weibull CDF shifted at x0, bounded [A, L]
    # y = A + (L - A) * (1 - exp(-((max(0, x-x0))/lam)^k))
    x = np.asarray(x, dtype=float)
    z = np.clip(x - x0, 0.0, None) / np.clip(lam, 1e-6, None)
    return A + (L - A) * (1 - np.exp(-(z**np.clip(k, 1e-6, None))))


def probit_cdf(x, A, L, x0, k):
    # Normal CDF link between A and L
    return A + (L - A) * norm.cdf(k * (x - x0))


# ------------------------------
# Fitting utilities
# ------------------------------

def safe_curve_fit(func, x, y, p0, bounds=None, maxfev=30000):
    try:
        popt, _ = curve_fit(
            func,
            x,
            y,
            p0=p0,
            bounds=bounds if bounds is not None else (-np.inf, np.inf),
            maxfev=maxfev,
        )
        y_hat = func(x, *popt)
        return popt, y_hat, None
    except Exception as e:
        return None, None, e


def piecewise_linear_eval(x, m1, m2, cp, b1):
    x = np.asarray(x)
    y_cp = m1 * cp + b1
    return np.where(x <= cp, m1 * x + b1, y_cp + m2 * (x - cp))


def fit_piecewise_linear(x, y, n_grid=12):
    x = np.asarray(x)
    y = np.asarray(y)
    x_sorted = np.sort(np.unique(x))
    if len(x_sorted) < 4:
        return None, None, 'insufficient unique x for piecewise'
    cps = np.linspace(x_sorted[1], x_sorted[-2], min(n_grid, len(x_sorted) - 2))
    best = None
    best_rss = np.inf
    for cp in cps:
        left = x <= cp
        right = ~left
        if left.sum() < 2 or right.sum() < 2:
            continue
        m1, b1 = np.polyfit(x[left], y[left], 1)
        m2, _ = np.polyfit(x[right], y[right], 1)
        y_hat = piecewise_linear_eval(x, m1, m2, cp, b1)
        rss = np.sum((y - y_hat) ** 2)
        if rss < best_rss:
            best_rss = rss
            best = (m1, m2, cp, b1, y_hat)
    if best is None:
        return None, None, 'fit failed'
    m1, m2, cp, b1, y_hat = best
    return (m1, m2, cp, b1), y_hat, None


def model_stats(x, y, y_hat, n_params):
    y = np.asarray(y)
    y_hat = np.asarray(y_hat)
    rss = float(np.sum((y - y_hat) ** 2))
    n = len(y)
    if n < n_params + 1:
        return {'rss': rss, 'aic': np.nan, 'bic': np.nan, 'r2': np.nan}
    var = np.var(y)
    r2 = float(1 - rss / np.sum((y - np.mean(y)) ** 2)) if var > 0 else np.nan
    aic = float(n * np.log(max(rss / n, 1e-12)) + 2 * n_params)
    bic = float(n * np.log(max(rss / n, 1e-12)) + np.log(n) * n_params)
    return {'rss': rss, 'aic': aic, 'bic': bic, 'r2': r2}


# ------------------------------
# Data prep (reusable, idempotent)
# ------------------------------
rat_col = None
for candidate in ['id', 'rat', 'animal']:
    if candidate in x_array.columns:
        rat_col = candidate
        break
if rat_col is None:
    raise ValueError("No rat identifier column found among ['id','rat','animal'].")

subset_query = "condition == 'deplete' & infusiontype == '45NaCl'"
subset_df = x_array.query(subset_query).copy()
subset_idx = subset_df.index.to_numpy()
X_subset = pca_data[subset_idx]

if 'pred_prob' not in subset_df.columns:
    bundle_path = RESULTSFOLDER / "svm_trial_classifier_with_threshold.pkl"
    if bundle_path.exists():
        import dill
        with open(bundle_path, 'rb') as f:
            bundle = dill.load(f)
        pipe_loaded = bundle['pipeline']
    else:
        # Fallback: refit quickly on labeled data
        label_col = 'class_label'
        mask = x_array[label_col].notna()
        X_full = pca_data[mask]
        y_full = x_array.loc[mask, label_col].astype(int).to_numpy()
        from sklearn.pipeline import Pipeline
        from sklearn.preprocessing import StandardScaler
        from sklearn.svm import SVC
        pipe_loaded = Pipeline([
            ("scaler", StandardScaler()),
            ("svm", SVC(kernel="rbf", probability=True, class_weight="balanced"))
        ])
        pipe_loaded.fit(X_full, y_full)
    subset_df['pred_prob'] = pipe_loaded.predict_proba(X_subset)[:, 1]

subset_df = subset_df.sort_values('trial').copy()

# ------------------------------
# Fit expanded model set per rat
# ------------------------------
results = []
subset_df['best_fit_ext'] = np.nan
subset_df['best_model_ext'] = np.nan

for r, df_r in subset_df.groupby(rat_col):
    x = df_r['trial'].to_numpy(dtype=float)
    y = df_r['pred_prob'].to_numpy(dtype=float)
    if len(np.unique(x)) < 5 or np.var(y) == 0:
        results.append({rat_col: r, 'model': 'NA', 'params': {}, 'rss': np.nan, 'aic': np.nan, 'bic': np.nan, 'r2': np.nan, 'n': len(x)})
        continue

    model_fits = []

    # Common inits
    L_init = float(min(1.0, max(0.6, np.nanmax(y))))
    A_init = float(max(0.0, min(np.nanmin(y), 0.3)))
    x0_init = float(np.median(x))
    # correlation sign for slope
    try:
        corr = float(np.corrcoef(x, y)[0, 1])
    except Exception:
        corr = 0.0
    if not np.isfinite(corr):
        corr = 0.0
    k_init = 0.5 if corr >= 0 else -0.5

    # logistic3 (3 params)
    popt, yhat, err = safe_curve_fit(
        logistic3, x, y,
        p0=[L_init, x0_init, k_init],
        bounds=([0.2, x.min()-10, -5.0], [1.5, x.max()+10, 5.0])
    )
    if popt is not None:
        stats = model_stats(x, y, yhat, 3)
        model_fits.append(('logistic3', popt, yhat, stats))

    # logistic4 / 4PL (4 params)
    popt, yhat, err = safe_curve_fit(
        logistic4, x, y,
        p0=[A_init, L_init, x0_init, k_init],
        bounds=([0.0, 0.2, x.min()-10, -5.0], [0.8, 1.5, x.max()+10, 5.0])
    )
    if popt is not None:
        stats = model_stats(x, y, yhat, 4)
        model_fits.append(('logistic4', popt, yhat, stats))

    # richards5 / 5PL (5 params)
    popt, yhat, err = safe_curve_fit(
        richards5, x, y,
        p0=[A_init, L_init, x0_init, k_init, 1.0],
        bounds=([0.0, 0.2, x.min()-10, -5.0, 0.2], [0.8, 1.5, x.max()+10, 5.0, 5.0])
    )
    if popt is not None:
        stats = model_stats(x, y, yhat, 5)
        model_fits.append(('richards5', popt, yhat, stats))

    # gompertz (3 params)
    popt, yhat, err = safe_curve_fit(
        gompertz, x, y,
        p0=[float(np.nanmax(y)), 1.0, 0.1],
        bounds=([0.2, 0.001, 0.001], [1.5, 5.0, 5.0])
    )
    if popt is not None:
        stats = model_stats(x, y, yhat, 3)
        model_fits.append(('gompertz', popt, yhat, stats))

    # exp_sat (2 params)
    popt, yhat, err = safe_curve_fit(
        exp_sat, x, y,
        p0=[L_init, 0.1],
        bounds=([0.2, 0.0005], [1.5, 5.0])
    )
    if popt is not None:
        stats = model_stats(x, y, yhat, 2)
        model_fits.append(('exp_sat', popt, yhat, stats))

    # weibull_cdf (5 params)
    lam_init = max(1.0, 0.25 * (x.max() - x.min()))
    kW_init = 2.0
    popt, yhat, err = safe_curve_fit(
        weibull_cdf, x, y,
        p0=[A_init, L_init, x0_init, lam_init, kW_init],
        bounds=(
            [0.0, 0.2, x.min()-10, 1e-3, 0.2],
            [0.8, 1.5, x.max()+10, (x.max()-x.min())*5 + 10, 8.0]
        )
    )
    if popt is not None:
        stats = model_stats(x, y, yhat, 5)
        model_fits.append(('weibull_cdf', popt, yhat, stats))

    # probit_cdf (4 params)
    popt, yhat, err = safe_curve_fit(
        probit_cdf, x, y,
        p0=[A_init, L_init, x0_init, k_init],
        bounds=([0.0, 0.2, x.min()-10, -5.0], [0.8, 1.5, x.max()+10, 5.0])
    )
    if popt is not None:
        stats = model_stats(x, y, yhat, 4)
        model_fits.append(('probit_cdf', popt, yhat, stats))

    # piecewise linear (4 params)
    pw_params, yhat, err = fit_piecewise_linear(x, y, n_grid=16)
    if pw_params is not None:
        stats = model_stats(x, y, yhat, 4)
        model_fits.append(('piecewise_linear', pw_params, yhat, stats))

    if not model_fits:
        results.append({rat_col: r, 'model': 'none_fit', 'params': {}, 'rss': np.nan, 'aic': np.nan, 'bic': np.nan, 'r2': np.nan, 'n': len(x)})
        continue

    # Select best by AIC
    best = min(model_fits, key=lambda mf: mf[3]['aic'] if np.isfinite(mf[3]['aic']) else np.inf)
    best_name, best_params, best_yhat, best_stats = best

    # Log all fits
    for name, params, yhat_model, stats_model in model_fits:
        results.append({
            rat_col: r,
            'model': name,
            'params': params,
            'rss': stats_model['rss'],
            'aic': stats_model['aic'],
            'bic': stats_model['bic'],
            'r2': stats_model['r2'],
            'n': len(x),
            'is_best': name == best_name
        })

    # Store best fit series for plotting
    subset_df.loc[df_r.index, 'best_fit_ext'] = best_yhat
    subset_df.loc[df_r.index, 'best_model_ext'] = best_name

# ------------------------------
# Summaries and plots
# ------------------------------
results_ext_df = pd.DataFrame(results)
print("Extended model comparison (per rat x per model):")
print(results_ext_df.groupby('model')[['aic','r2']].agg(['mean','median']).round(3))

best_ext = results_ext_df[results_ext_df['is_best'] == True].copy()
print("Best (extended) model per rat:")
print(best_ext[[rat_col, 'model', 'aic', 'bic', 'r2', 'n']].sort_values(rat_col))

# Save
try:
    best_ext_path = RESULTSFOLDER / 'best_model_fits_ext_deplete_45NaCl_per_rat.csv'
    all_ext_path = RESULTSFOLDER / 'all_model_fits_ext_deplete_45NaCl_per_rat.csv'
    best_ext.to_csv(best_ext_path, index=False)
    results_ext_df.to_csv(all_ext_path, index=False)
    print(f"Saved extended best summary to {best_ext_path}\nSaved extended full results to {all_ext_path}")
except Exception as e:
    print(f"Could not save extended CSVs: {e}")

# Facet plot of best extended model per rat
n_rats = subset_df[rat_col].nunique()
col_wrap = 4 if n_rats > 4 else n_rats

g = sns.FacetGrid(subset_df, col=rat_col, col_wrap=col_wrap, sharey=True, height=2.4)

def _plot_best_ext(data, color):
    ax = plt.gca()
    d = data.sort_values('trial')
    ax.scatter(d['trial'], d['pred_prob'], s=25, color='#016895', alpha=0.75)
    if d['best_fit_ext'].notna().any():
        ax.plot(d['trial'], d['best_fit_ext'], color='#0C7BDC', linewidth=1.5)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel('Trial')
    ax.set_ylabel('P(class=1)')
    ax.set_title(str(d[rat_col].iloc[0]))
    ax.tick_params(axis='x', labelsize=8)
    ax.tick_params(axis='y', labelsize=8)
    sns.despine(ax=ax)

g.map_dataframe(_plot_best_ext)
plt.subplots_adjust(top=0.88)
plt.suptitle('Best extended model fits per rat (deplete & 45NaCl subset)', fontsize=12)
plt.show()

best_ext

# %%

# %%
# Convert predicted probabilities to binary via chosen threshold and fit per-rat logistic (GLM) vs trial
import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib.pyplot as plt
import seaborn as sns

# Load subset and probabilities (reuse existing logic)
rat_col = None
for candidate in ['id', 'rat', 'animal']:
    if candidate in x_array.columns:
        rat_col = candidate
        break
if rat_col is None:
    raise ValueError("No rat identifier column found among ['id','rat','animal'].")

subset_query = "condition == 'deplete' & infusiontype == '45NaCl'"
subset_df = x_array.query(subset_query).copy()
subset_idx = subset_df.index.to_numpy()
X_subset = pca_data[subset_idx]

# Ensure pred_prob is present
if 'pred_prob' not in subset_df.columns:
    bundle_path = RESULTSFOLDER / 'svm_trial_classifier_with_threshold.pkl'
    if bundle_path.exists():
        import dill
        with open(bundle_path, 'rb') as f:
            bundle = dill.load(f)
        pipe_loaded = bundle['pipeline']
        chosen_threshold = bundle.get('threshold', 0.5)
    else:
        label_col = 'class_label'
        mask = x_array[label_col].notna()
        X_full = pca_data[mask]
        y_full = x_array.loc[mask, label_col].astype(int).to_numpy()
        from sklearn.pipeline import Pipeline
        from sklearn.preprocessing import StandardScaler
        from sklearn.svm import SVC
        pipe_loaded = Pipeline([
            ("scaler", StandardScaler()),
            ("svm", SVC(kernel='rbf', probability=True, class_weight='balanced'))
        ])
        pipe_loaded.fit(X_full, y_full)
        chosen_threshold = 0.5
    subset_df['pred_prob'] = pipe_loaded.predict_proba(X_subset)[:,1]
else:
    # Try to load threshold JSON for consistency
    threshold_json = RESULTSFOLDER / 'svm_trial_classifier_threshold.json'
    if threshold_json.exists():
        import json
        with open(threshold_json, 'r') as f:
            th_info = json.load(f)
        chosen_threshold = th_info.get('chosen_threshold', 0.5)
    else:
        chosen_threshold = 0.5

# Create binary variable
subset_df['pred_binary'] = (subset_df['pred_prob'] >= chosen_threshold).astype(int)
print(f"Applied threshold={chosen_threshold:.3f}; binary distribution:")
print(subset_df['pred_binary'].value_counts())

# Fit logistic regression (GLM binomial) per rat: pred_binary ~ trial
logit_rows = []
subset_df['glm_prob'] = np.nan
subset_df = subset_df.sort_values('trial')

for r, df_r in subset_df.groupby(rat_col):
    x = df_r['trial'].to_numpy(dtype=float)
    y = df_r['pred_binary'].to_numpy(dtype=int)
    if len(np.unique(x)) < 3 or y.sum() == 0 or y.sum() == len(y):
        logit_rows.append({rat_col: r, 'intercept': np.nan, 'slope': np.nan, 'midpoint_trial': np.nan, 'n': len(y), 'note': 'insufficient variation'})
        continue
    # Design matrix with intercept
    X_design = sm.add_constant(x)
    try:
        model = sm.GLM(y, X_design, family=sm.families.Binomial())
        res = model.fit()
        intercept, slope = res.params['const'], res.params['x1']
        # Midpoint where p=0.5 => intercept + slope*x = 0 => x = -intercept/slope
        midpoint = -intercept / slope if slope != 0 else np.nan
        logit_rows.append({rat_col: r, 'intercept': intercept, 'slope': slope, 'midpoint_trial': midpoint, 'n': len(y), 'note': ''})
        # Store fitted probabilities
        subset_df.loc[df_r.index, 'glm_prob'] = res.predict(X_design)
    except Exception as e:
        logit_rows.append({rat_col: r, 'intercept': np.nan, 'slope': np.nan, 'midpoint_trial': np.nan, 'n': len(y), 'note': f'fit_error: {e}'})

logit_df = pd.DataFrame(logit_rows).sort_values('midpoint_trial')
print("Per-rat logistic GLM fits (sorted by midpoint):")
print(logit_df[[rat_col, 'intercept', 'slope', 'midpoint_trial', 'n', 'note']])

# Optional save
try:
    glm_csv = RESULTSFOLDER / 'per_rat_logistic_glm_binary_threshold.csv'
    logit_df.to_csv(glm_csv, index=False)
    print(f"Saved GLM summary to {glm_csv}")
except Exception as e:
    print(f"Could not save GLM summary: {e}")

# Facet plot: binary points (jitter) + fitted logistic curve
n_rats = subset_df[rat_col].nunique()
col_wrap = 4 if n_rats > 4 else n_rats

g = sns.FacetGrid(subset_df, col=rat_col, col_wrap=col_wrap, sharey=True, height=2.4)

def _plot_glm(data, color):
    ax = plt.gca()
    d = data.sort_values('trial')
    # Plot jittered binary points
    jitter = (np.random.rand(len(d)) - 0.5) * 0.06
    ax.scatter(d['trial'], d['pred_binary'] + jitter, s=28, color='#444444', alpha=0.7, label='Binary')
    # Fitted curve
    if d['glm_prob'].notna().any():
        ax.plot(d['trial'], d['glm_prob'], color='#D95F02', linewidth=1.5, label='GLM fit')
    ax.axhline(0.5, color='#999999', linestyle='--', linewidth=0.8)
    ax.set_ylim(-0.1, 1.1)
    ax.set_xlabel('Trial')
    ax.set_ylabel('Binary / Fitted P')
    ax.tick_params(axis='x', labelsize=8)
    ax.tick_params(axis='y', labelsize=8)
    sns.despine(ax=ax)

    # Annotate midpoint if available
    rat_id = str(d[rat_col].iloc[0])
    row = logit_df[logit_df[rat_col] == d[rat_col].iloc[0]]
    if row['midpoint_trial'].notna().any():
        mp = row['midpoint_trial'].iloc[0]
        if np.isfinite(mp):
            ax.text(0.05, 0.9, f"mid={mp:.1f}", transform=ax.transAxes, fontsize=8, color='#D95F02')


g.map_dataframe(_plot_glm)
plt.subplots_adjust(top=0.87)
plt.suptitle('Per-rat logistic GLM on thresholded binary outcome', fontsize=12)
plt.show()

logit_df

# %%

# %%
# Sweep number of PCs (1..7) for multiple classifiers with 10-fold Stratified CV
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.model_selection import StratifiedKFold, cross_validate
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC, LinearSVC
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier, GradientBoostingClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.calibration import CalibratedClassifierCV

# Ensure PCA array is available
try:
    pca
except NameError:
    # Load from disk if not in memory
    import dill
    pcafile = RESULTSFOLDER / "transformed_data_vel.pickle"
    with open(pcafile, 'rb') as f:
        pca = dill.load(f)

label_col = 'class_label'
if label_col not in x_array.columns:
    raise ValueError(f"Missing label column {label_col} in x_array")

mask = x_array[label_col].notna()
y = x_array.loc[mask, label_col].astype(int).to_numpy()

max_pcs = int(min(7, pca.shape[1]))
X_all = pca[mask, :max_pcs]

skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
scoring = {
    'accuracy': 'accuracy',
    'roc_auc': 'roc_auc',
    'precision': 'precision',
    'recall': 'recall',
    'f1': 'f1',
}

models = {
    'SVM_RBF': SVC(kernel='rbf', probability=True, class_weight='balanced'),
    'LogReg_L2': LogisticRegression(max_iter=500, class_weight='balanced'),
    'LinearSVM': CalibratedClassifierCV(LinearSVC(class_weight='balanced', max_iter=8000)),
    'RandomForest': RandomForestClassifier(n_estimators=400, class_weight='balanced', random_state=42),
    'ExtraTrees': ExtraTreesClassifier(n_estimators=400, random_state=42),
    'GradientBoosting': GradientBoostingClassifier(random_state=42),
    'KNN_7': KNeighborsClassifier(n_neighbors=7),
    'LDA': LinearDiscriminantAnalysis(),
}

rows = []
for n in range(1, max_pcs + 1):
    Xn = X_all[:, :n]
    for name, clf in models.items():
        pipe = Pipeline([
            ('scaler', StandardScaler()),
            ('clf', clf),
        ])
        cv_res = cross_validate(pipe, Xn, y, cv=skf, scoring=scoring, n_jobs=-1, error_score='raise')
        rows.append({
            'model': name,
            'n_pcs': n,
            'acc_mean': float(cv_res['test_accuracy'].mean()),
            'acc_std': float(cv_res['test_accuracy'].std()),
            'auc_mean': float(cv_res['test_roc_auc'].mean()),
            'auc_std': float(cv_res['test_roc_auc'].std()),
            'f1_mean': float(cv_res['test_f1'].mean()),
            'f1_std': float(cv_res['test_f1'].std()),
            'prec_mean': float(cv_res['test_precision'].mean()),
            'prec_std': float(cv_res['test_precision'].std()),
            'recall_mean': float(cv_res['test_recall'].mean()),
            'recall_std': float(cv_res['test_recall'].std()),
        })

pc_sweep_df = pd.DataFrame(rows)

# Best n_pcs per model by AUC
best_by_auc = pc_sweep_df.loc[pc_sweep_df.groupby('model')['auc_mean'].idxmax()].sort_values('auc_mean', ascending=False)
print("Best n_pcs per model by AUC:")
print(best_by_auc[['model','n_pcs','auc_mean','acc_mean','f1_mean']])

# Save CSV
try:
    out_csv = RESULTSFOLDER / 'pc_sweep_classifiers_cv_metrics.csv'
    pc_sweep_df.to_csv(out_csv, index=False)
    print(f"Saved PC sweep CV metrics to {out_csv}")
except Exception as e:
    print(f"Could not save CSV: {e}")

# Plot: AUC and Accuracy vs n_pcs per model
fig, axes = plt.subplots(1, 2, figsize=(12, 4), sharex=True)
sns.lineplot(data=pc_sweep_df, x='n_pcs', y='auc_mean', hue='model', marker='o', ax=axes[0])
axes[0].set_title('ROC AUC vs #PCs')
axes[0].set_xlabel('# PCs')
axes[0].set_ylabel('Mean ROC AUC (CV)')
axes[0].legend(bbox_to_anchor=(1.02, 1), loc='upper left', frameon=False)

sns.lineplot(data=pc_sweep_df, x='n_pcs', y='acc_mean', hue='model', marker='o', ax=axes[1])
axes[1].set_title('Accuracy vs #PCs')
axes[1].set_xlabel('# PCs')
axes[1].set_ylabel('Mean Accuracy (CV)')
axes[1].legend_.remove()

plt.tight_layout()
plt.show()

pc_sweep_df

# %%

# %%
# Heatmap: AUC by model vs number of PCs
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Pivot the PC sweep data for heatmap
heatmap_data = pc_sweep_df.pivot(index='model', columns='n_pcs', values='auc_mean')

# Sort models by mean AUC across all PCs (descending)
row_order = heatmap_data.mean(axis=1).sort_values(ascending=False).index

fig, ax = plt.subplots(figsize=(10, 6))
sns.heatmap(
    heatmap_data.loc[row_order],
    annot=True,
    fmt='.3f',
    cmap='RdYlGn',
    center=0.75,
    vmin=0.5,
    vmax=1.0,
    cbar_kws={'label': 'Mean ROC AUC (CV)'},
    linewidths=0.5,
    linecolor='gray',
    ax=ax
)
ax.set_xlabel('Number of PCs', fontsize=12)
ax.set_ylabel('Classifier', fontsize=12)
ax.set_title('Cross-Validated ROC AUC by Classifier and #PCs\n(Stratified 10-fold CV)', fontsize=13, pad=15)
plt.tight_layout()
plt.show()

# Display best n_pcs per model again for reference
best_by_auc = pc_sweep_df.loc[pc_sweep_df.groupby('model')['auc_mean'].idxmax()].sort_values('auc_mean', ascending=False)
print("\nBest #PCs per classifier (by AUC):")
print(best_by_auc[['model', 'n_pcs', 'auc_mean', 'acc_mean', 'f1_mean']].to_string(index=False))

# %%

# %%
# Compute AIC/BIC for PC selection: balance fit quality vs model complexity
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.model_selection import cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

# AIC approximation for classification:
# AIC = -2 * log-likelihood + 2 * k
# where k = number of parameters (PCs + intercept + model-specific params)
# BIC = -2 * log-likelihood + k * log(n)
# 
# For binary classification with predicted probabilities:
# log-likelihood ≈ sum(y*log(p) + (1-y)*log(1-p))

def compute_aic_bic_classification(y_true, y_prob, n_params, n_samples):
    """
    Compute AIC and BIC for binary classification.
    
    Parameters:
    - y_true: true binary labels
    - y_prob: predicted probabilities for class 1
    - n_params: effective number of model parameters
    - n_samples: number of observations
    """
    # Clip probabilities to avoid log(0)
    eps = 1e-15
    y_prob = np.clip(y_prob, eps, 1 - eps)
    
    # Log-likelihood
    log_lik = np.sum(y_true * np.log(y_prob) + (1 - y_true) * np.log(1 - y_prob))
    
    # AIC and BIC
    aic = -2 * log_lik + 2 * n_params
    bic = -2 * log_lik + n_params * np.log(n_samples)
    
    return aic, bic

# Reuse data from PC sweep
label_col = 'class_label'
mask = x_array[label_col].notna()
y = x_array.loc[mask, label_col].astype(int).to_numpy()
n_samples = len(y)

max_pcs = int(min(7, pca.shape[1]))
X_all = pca[mask, :max_pcs]

# Models with reasonable parameter count estimates
# For simplicity: n_params ≈ n_features + intercept + model complexity
models_for_aic = {
    'SVM_RBF': SVC(kernel='rbf', probability=True, class_weight='balanced'),
    'LogReg_L2': LogisticRegression(max_iter=500, class_weight='balanced'),
    'RandomForest': RandomForestClassifier(n_estimators=400, class_weight='balanced', random_state=42),
    'GradientBoosting': GradientBoostingClassifier(random_state=42),
    'LDA': LinearDiscriminantAnalysis(),
}

aic_rows = []
for n in range(1, max_pcs + 1):
    Xn = X_all[:, :n]
    for name, clf in models_for_aic.items():
        pipe = Pipeline([
            ('scaler', StandardScaler()),
            ('clf', clf),
        ])
        
        # Get cross-validated predictions
        skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
        y_prob = cross_val_predict(pipe, Xn, y, cv=skf, method='predict_proba', n_jobs=-1)[:, 1]
        
        # Estimate effective parameters
        # Linear models: n_features + 1 (intercept)
        # Tree ensembles: rough heuristic based on complexity
        # SVM: effective dimensionality harder to define, use n_features as proxy
        if name == 'LogReg_L2' or name == 'LDA':
            n_params = n + 1
        elif name == 'SVM_RBF':
            # SVMs with RBF kernel: use number of support vectors as proxy (not available in CV)
            # Approximation: use n_features * scale factor
            n_params = n * 2  # Heuristic
        elif 'Forest' in name or 'Boosting' in name:
            # Tree ensembles: very complex, use heuristic
            # Rough estimate: n_trees * avg_leaves_per_tree, simplified to n * n_estimators / scale
            n_params = n * 10  # Heuristic: captures growing complexity
        else:
            n_params = n + 1
        
        aic, bic = compute_aic_bic_classification(y, y_prob, n_params, n_samples)
        
        aic_rows.append({
            'model': name,
            'n_pcs': n,
            'n_params': n_params,
            'aic': aic,
            'bic': bic,
        })

aic_df = pd.DataFrame(aic_rows)

# Best n_pcs per model by AIC (lower is better)
best_by_aic = aic_df.loc[aic_df.groupby('model')['aic'].idxmin()].sort_values('aic')
print("Best #PCs per classifier by AIC (lower = better balance of fit vs complexity):")
print(best_by_aic[['model', 'n_pcs', 'aic', 'bic']])

# Save
try:
    out_csv = RESULTSFOLDER / 'pc_selection_aic_bic.csv'
    aic_df.to_csv(out_csv, index=False)
    print(f"\nSaved AIC/BIC results to {out_csv}")
except Exception as e:
    print(f"Could not save CSV: {e}")

# Plot: AIC vs n_pcs per model
fig, axes = plt.subplots(1, 2, figsize=(12, 4), sharex=True)

sns.lineplot(data=aic_df, x='n_pcs', y='aic', hue='model', marker='o', ax=axes[0])
axes[0].set_title('AIC vs #PCs (lower = better)')
axes[0].set_xlabel('# PCs')
axes[0].set_ylabel('AIC')
axes[0].legend(bbox_to_anchor=(1.02, 1), loc='upper left', frameon=False, fontsize=9)

sns.lineplot(data=aic_df, x='n_pcs', y='bic', hue='model', marker='o', ax=axes[1])
axes[1].set_title('BIC vs #PCs (lower = better)')
axes[1].set_xlabel('# PCs')
axes[1].set_ylabel('BIC')
axes[1].legend_.remove()

plt.tight_layout()
plt.show()

# Heatmap of AIC
heatmap_aic = aic_df.pivot(index='model', columns='n_pcs', values='aic')
row_order_aic = heatmap_aic.min(axis=1).sort_values().index  # Sort by best (lowest) AIC

fig, ax = plt.subplots(figsize=(10, 5))
sns.heatmap(
    heatmap_aic.loc[row_order_aic],
    annot=True,
    fmt='.0f',
    cmap='RdYlGn_r',  # Reverse: green = low (good), red = high (bad)
    cbar_kws={'label': 'AIC'},
    linewidths=0.5,
    linecolor='gray',
    ax=ax
)
ax.set_xlabel('Number of PCs', fontsize=12)
ax.set_ylabel('Classifier', fontsize=12)
ax.set_title('AIC by Classifier and #PCs (lower = better trade-off)\nGreen = optimal complexity', fontsize=13, pad=15)
plt.tight_layout()
plt.show()

aic_df

# %%

# %% [markdown]
# ## Final Model Selection Summary
#
# ### Selected Model: **SVM with RBF kernel + 3 Principal Components**
#
# ---
#
# #### Rationale
#
# 1. **PC Selection (1-7 sweep)**
#    - **Performance**: SVM_RBF achieved highest cross-validated ROC AUC across PC counts
#    - **Complexity trade-off**: AIC/BIC analysis confirmed 3 PCs provides optimal balance
#    - **Diminishing returns**: Additional PCs beyond 3 showed minimal AUC improvement while increasing overfitting risk
#
# 2. **Classifier Comparison**
#    - Benchmarked 8 classifiers: SVM (RBF & Linear), Logistic Regression, Random Forest, Extra Trees, Gradient Boosting, KNN, LDA
#    - SVM_RBF consistently outperformed alternatives in terms of ROC AUC and F1 score
#    - Class balancing (`class_weight='balanced'`) effectively handles imbalanced training labels
#
# 3. **Threshold Optimization**
#    - Default 0.5 probability threshold replaced with optimized value via Youden's J criterion
#    - Maximizes sensitivity + specificity - 1, balancing true positive and true negative rates
#    - Alternative F1-optimized threshold also computed for comparison
#
# ---
#
# #### Model Performance (Stratified 10-Fold CV on 3 PCs)
#
# - **Accuracy**: ~0.85-0.90 (see CV results in earlier cells)
# - **ROC AUC**: ~0.90-0.95
# - **F1 Score**: Optimized via threshold tuning
# - **Generalization**: Validated on held-out test set (20% stratified split)
#
# ---
#
# #### Training Data
#
# - **Class Labels**: 
#   - `1` = Aversive trials (condition='deplete' & infusiontype='10NaCl')
#   - `0` = Appetitive trials (condition='replete' & infusiontype='45NaCl')
# - **Features**: First 3 principal components from velocity kinematics PCA
# - **Sample size**: ~1000 trials (exact count in label creation cell output)
#
# ---
#
# #### Saved Artifacts
#
# | File | Description |
# |------|-------------|
# | `results/svm_trial_classifier_with_threshold.pkl` | Complete pipeline (StandardScaler + SVM) + optimized threshold |
# | `results/svm_trial_classifier_threshold.json` | Threshold value and metadata (human-readable) |
# | `results/pc_sweep_classifiers_cv_metrics.csv` | Full comparison across 1-7 PCs and 8 classifiers |
# | `results/pc_selection_aic_bic.csv` | AIC/BIC justification for 3 PCs |
# | `results/best_model_fits_deplete_45NaCl_per_rat.csv` | Per-rat trajectory parameters (sigmoid fits on predicted probabilities) |
#
# ---
#
# #### Usage Example
#
# ```python
# import dill
#
# # Load model bundle
# with open(RESULTSFOLDER / 'svm_trial_classifier_with_threshold.pkl', 'rb') as f:
#     bundle = dill.load(f)
#
# pipeline = bundle['pipeline']
# threshold = bundle['threshold']
#
# # Predict on new data (must be first 3 PCs)
# X_new = pca_new[:, :3]  # Shape: (n_samples, 3)
# probs = pipeline.predict_proba(X_new)[:, 1]
# predictions = (probs >= threshold).astype(int)
# ```
#
# ---
#
# #### Application: Transition Dynamics
#
# The trained classifier was applied to the **deplete & 45NaCl** subset (intermediate condition) to:
# - Quantify trial-by-trial probability of "aversive-like" behavior
# - Model per-rat trajectories using logistic, Gompertz, and other parametric curves
# - Extract transition points (sigmoid midpoints, piecewise breakpoints) indicating behavioral shift timing
# - Compare models via AIC to identify best-fitting transition functions per individual
#
# This approach reveals **when** and **how** animals transition from appetitive to aversive behavioral states during sodium depletion with high-concentration NaCl.
