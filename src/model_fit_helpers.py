"""Reusable helpers for trial-wise model fitting and model selection."""

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import pearsonr


def _linear_model(x, m, b):
    return m * x + b


def _exp_model(x, a, b, c):
    return a * np.exp(b * x) + c


def _sigmoid_model(x, L, k, x0, b):
    z = np.clip(-k * (x - x0), -60, 60)
    return L / (1 + np.exp(z)) + b


def _logistic4_model(x, A, L, x0, k):
    z = np.clip(-k * (x - x0), -60, 60)
    return A + (L - A) / (1 + np.exp(z))


def _logistic3_model(x, L, x0, k):
    z = np.clip(-k * (x - x0), -60, 60)
    return L / (1 + np.exp(z))


def binarize_series_threshold(values, low=-0.7, high=0.7):
    """Binarize a series with a middle exclusion zone.

    Values >= high map to 1, values <= low map to 0, and middle values become NaN.
    """
    arr = np.asarray(values, dtype=float)
    out = np.full(arr.shape, np.nan, dtype=float)
    out[arr >= high] = 1.0
    out[arr <= low] = 0.0
    return out


def quantize_series_threshold(values, low=-0.7, high=0.7):
    """Quantize a series to {-1, 0, +1} using low/high thresholds."""
    arr = np.asarray(values, dtype=float)
    out = np.zeros(arr.shape, dtype=float)
    out[arr >= high] = 1.0
    out[arr <= low] = -1.0
    return out


def _safe_pearson(y_true, y_pred):
    if np.allclose(np.std(y_true), 0) or np.allclose(np.std(y_pred), 0):
        return np.nan, np.nan
    return pearsonr(y_true, y_pred)


def _model_metrics(y_true, y_pred, n_params):
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    n = len(y_true)
    residuals = y_true - y_pred
    sse = np.nansum(residuals ** 2)
    sse = max(float(sse), 1e-12)
    rmse = np.sqrt(sse / n) if n > 0 else np.nan
    aic = n * np.log(sse / n) + 2 * n_params if n > 0 else np.nan
    bic = n * np.log(sse / n) + n_params * np.log(n) if n > 1 else np.nan
    aicc = aic + (2 * n_params * (n_params + 1)) / (n - n_params - 1) if n > (n_params + 1) else np.nan
    return rmse, aic, aicc, bic


def _fit_linear(x, y):
    try:
        m, b = np.polyfit(x, y, 1)
        yhat = _linear_model(x, m, b)
        r, p = _safe_pearson(y, yhat)
        rmse, aic, aicc, bic = _model_metrics(y, yhat, n_params=2)
        return {"r": r, "p": p, "rmse": rmse, "aic": aic, "aicc": aicc, "bic": bic, "params": np.array([m, b]), "pcov": None}
    except Exception:
        return {"r": np.nan, "p": np.nan, "rmse": np.nan, "aic": np.nan, "aicc": np.nan, "bic": np.nan, "params": np.array([np.nan, np.nan]), "pcov": None}


def _fit_curve(model_name, model_func, x, y):
    try:
        if model_name == "exponential":
            p0 = [np.ptp(y) if np.ptp(y) != 0 else 1.0, 0.05, np.min(y)]
            bounds = ([-np.inf, -np.inf, -np.inf], [np.inf, np.inf, np.inf])
            n_params = 3
        elif model_name == "sigmoidal":
            p0 = [np.ptp(y) if np.ptp(y) != 0 else 1.0, 0.2, np.median(x), np.min(y)]
            bounds = ([-np.inf, -np.inf, np.min(x), -np.inf], [np.inf, np.inf, np.max(x), np.inf])
            n_params = 4
        else:
            raise ValueError(f"Unknown model: {model_name}")

        params, pcov = curve_fit(model_func, x, y, p0=p0, bounds=bounds, maxfev=30000)
        yhat = model_func(x, *params)
        r, p = _safe_pearson(y, yhat)
        rmse, aic, aicc, bic = _model_metrics(y, yhat, n_params=n_params)
        return {"r": r, "p": p, "rmse": rmse, "aic": aic, "aicc": aicc, "bic": bic, "params": params, "pcov": pcov}
    except Exception:
        n_params = 3 if model_name == "exponential" else 4
        return {"r": np.nan, "p": np.nan, "rmse": np.nan, "aic": np.nan, "aicc": np.nan, "bic": np.nan, "params": np.full(n_params, np.nan), "pcov": None}


def fit_curve_series(x, y, model_name="sigmoidal", p0=None, bounds=None, maxfev=30000):
    """Fit one series with linear, exponential, or sigmoidal model.

    Returns a dictionary with fitted params, predictions, and correlation stats.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if model_name == "linear":
        res = _fit_linear(x, y)
        y_hat = _linear_model(x, *res["params"])
        return {
            "success": np.all(np.isfinite(res["params"])),
            "model": "linear",
            "params": res["params"],
            "pcov": None,
            "y_hat": y_hat,
            "r": res["r"],
            "p": res["p"],
        }

    if model_name not in {"exponential", "sigmoidal"}:
        raise ValueError("model_name must be one of: linear, exponential, sigmoidal")

    model_func = _exp_model if model_name == "exponential" else _sigmoid_model
    if p0 is None or bounds is None:
        seeded = _fit_curve(model_name, model_func, x, y)
        params = seeded["params"]
        if np.any(~np.isfinite(params)):
            return {
                "success": False,
                "model": model_name,
                "params": params,
                "pcov": None,
                "y_hat": None,
                "r": np.nan,
                "p": np.nan,
            }
        y_hat = model_func(x, *params)
        return {
            "success": True,
            "model": model_name,
            "params": params,
            "pcov": seeded["pcov"],
            "y_hat": y_hat,
            "r": seeded["r"],
            "p": seeded["p"],
        }

    try:
        params, pcov = curve_fit(model_func, x, y, p0=p0, bounds=bounds, maxfev=maxfev)
        y_hat = model_func(x, *params)
        r, p = _safe_pearson(y, y_hat)
        return {
            "success": True,
            "model": model_name,
            "params": params,
            "pcov": pcov,
            "y_hat": y_hat,
            "r": r,
            "p": p,
        }
    except Exception:
        n_params = 3 if model_name == "exponential" else 4
        return {
            "success": False,
            "model": model_name,
            "params": np.full(n_params, np.nan),
            "pcov": None,
            "y_hat": None,
            "r": np.nan,
            "p": np.nan,
        }


def _sigmoid_quality_checks(x, params, pcov):
    if params is None or len(params) != 4 or np.any(~np.isfinite(params)):
        return False, "fit_failed", {"x0_interior": False, "k_plausible": False, "ci_finite": False, "asymptotes_covered": False}

    L, k, x0, b = params
    x_min, x_max = float(np.min(x)), float(np.max(x))
    x_range = max(x_max - x_min, 1.0)
    edge_margin = 0.15 * x_range
    x0_interior = (x_min + edge_margin) <= x0 <= (x_max - edge_margin)
    # Keep a minimum steepness floor but allow very steep transitions.
    k_plausible = np.isfinite(k) and (0.02 <= abs(k))

    ci_finite = False
    if pcov is not None:
        diag = np.diag(pcov)
        if np.all(np.isfinite(diag)) and np.all(diag >= 0):
            se = np.sqrt(diag)
            ci_finite = np.all(np.isfinite(se)) and np.all(se > 0)

    amplitude = abs(L)
    if amplitude < 1e-8:
        asymptotes_covered = False
    else:
        lower_asym = min(b, b + L)
        upper_asym = max(b, b + L)
        y_start = _sigmoid_model(np.array([x_min]), L, k, x0, b)[0]
        y_end = _sigmoid_model(np.array([x_max]), L, k, x0, b)[0]

        tol = 0.2
        start_low = abs(y_start - lower_asym) / amplitude <= tol
        start_high = abs(y_start - upper_asym) / amplitude <= tol
        end_low = abs(y_end - lower_asym) / amplitude <= tol
        end_high = abs(y_end - upper_asym) / amplitude <= tol

        start_side = "low" if start_low else ("high" if start_high else None)
        end_side = "low" if end_low else ("high" if end_high else None)
        asymptotes_covered = (start_side is not None) and (end_side is not None) and (start_side != end_side)

    checks = {
        "x0_interior": bool(x0_interior),
        "k_plausible": bool(k_plausible),
        "ci_finite": bool(ci_finite),
        "asymptotes_covered": bool(asymptotes_covered),
    }

    failed = [name for name, ok in checks.items() if not ok]
    is_valid = len(failed) == 0
    reasons = "ok" if is_valid else ";".join(failed)
    return is_valid, reasons, checks


def fit_logistic_per_series(y, x=None, prefer_4p=True, direction=None, maxfev=60000):
    """Fit logistic curve to binary/near-binary data for one series.

    Parameters
    ----------
    y : array-like
        Response values, typically binary or bounded in [0, 1].
    x : array-like, optional
        Trial index values. If None, uses ``np.arange(len(y))``.
    prefer_4p : bool
        If True, attempt 4-parameter logistic first, then fallback to 3-parameter.
    direction : {"increasing", "decreasing"}, optional
        Optional expected direction to seed slope sign.
    maxfev : int
        Maximum function evaluations for ``curve_fit``.

    Returns
    -------
    dict
        Fit summary with model name, params, fitted trace, transition point, and flags.
    """
    if x is None:
        x = np.arange(len(y), dtype=float)
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    finite_mask = np.isfinite(y) & np.isfinite(x)
    if finite_mask.sum() < 4:
        return {
            "model": None,
            "x0_orig": np.nan,
            "k": np.nan,
            "r_squared": np.nan,
            "params": {},
            "y_hat": None,
            "x_fit": None,
            "y_fit": None,
            "success": False,
            "note": "too_few_points",
        }

    x = x[finite_mask]
    y = y[finite_mask]
    if np.unique(np.round(y, 6)).size < 2:
        return {
            "model": None,
            "x0_orig": np.nan,
            "k": np.nan,
            "r_squared": np.nan,
            "params": {},
            "y_hat": None,
            "x_fit": None,
            "y_fit": None,
            "success": False,
            "note": "no_variation",
        }

    x_mean, x_std = float(np.mean(x)), float(np.std(x))
    if not np.isfinite(x_std) or x_std == 0:
        x_std = 1.0
    x_norm = (x - x_mean) / x_std
    y_clip = np.clip(y, 1e-4, 1 - 1e-4)

    y_min, y_max = float(np.min(y_clip)), float(np.max(y_clip))
    A_init, L_init, x0_init = y_min, y_max, 0.0

    if direction is None:
        try:
            c = float(np.corrcoef(x, y_clip)[0, 1])
        except Exception:
            c = 0.0
        if not np.isfinite(c):
            c = 0.0
        sign = 1.0 if c >= 0 else -1.0
    else:
        sign = 1.0 if direction == "increasing" else -1.0

    k_mags = [0.5, 1.0, 2.0]
    trend_tol = 1e-6

    def direction_ok(y_hat):
        if direction is None:
            return True
        if direction == "decreasing":
            return bool(y_hat[0] > y_hat[-1] + trend_tol)
        if direction == "increasing":
            return bool(y_hat[-1] > y_hat[0] + trend_tol)
        return True

    def try_fit(func, p0_list, bounds):
        best, best_rss = None, np.inf
        for p0 in p0_list:
            try:
                popt, _ = curve_fit(func, x_norm, y_clip, p0=p0, bounds=bounds, maxfev=maxfev)
                y_hat = func(x_norm, *popt)
                rss = float(np.sum((y_clip - y_hat) ** 2))
                if rss < best_rss:
                    best_rss, best = rss, (popt, y_hat)
            except Exception:
                continue
        return best

    res4 = None
    if prefer_4p:
        p0s_4 = [[A_init, L_init, x0_init, sign * km] for km in k_mags]
        bnds_4 = ([-0.1, 0.4, -3.0, -10.0], [0.6, 1.6, 3.0, 10.0])
        res4 = try_fit(_logistic4_model, p0s_4, bnds_4)

    p0s_3 = [[L_init, x0_init, sign * km] for km in k_mags]
    bnds_3 = ([0.4, -3.0, -10.0], [1.6, 3.0, 10.0])
    res3 = try_fit(_logistic3_model, p0s_3, bnds_3)

    if res4 is not None:
        popt, y_hat = res4
        if not direction_ok(y_hat):
            return {
                "model": None,
                "x0_orig": np.nan,
                "k": np.nan,
                "r_squared": np.nan,
                "params": {},
                "y_hat": None,
                "x_fit": x,
                "y_fit": y,
                "success": False,
                "note": "direction_mismatch",
            }
        A, L, x0n, k = map(float, popt)
        x0_orig = x0n * x_std + x_mean
        ss_res = float(np.sum((y_clip - y_hat) ** 2))
        ss_tot = float(np.sum((y_clip - np.mean(y_clip)) ** 2))
        r_squared = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 0.0
        return {
            "model": "logistic4",
            "x0_orig": x0_orig,
            "k": k,
            "r_squared": r_squared,
            "params": {"A": A, "L": L, "x0_norm": x0n, "x0_orig": x0_orig, "k": k},
            "y_hat": y_hat,
            "x_fit": x,
            "y_fit": y,
            "success": True,
            "note": "",
        }

    if res3 is not None:
        popt, y_hat = res3
        if not direction_ok(y_hat):
            return {
                "model": None,
                "x0_orig": np.nan,
                "k": np.nan,
                "r_squared": np.nan,
                "params": {},
                "y_hat": None,
                "x_fit": x,
                "y_fit": y,
                "success": False,
                "note": "direction_mismatch",
            }
        L, x0n, k = map(float, popt)
        x0_orig = x0n * x_std + x_mean
        ss_res = float(np.sum((y_clip - y_hat) ** 2))
        ss_tot = float(np.sum((y_clip - np.mean(y_clip)) ** 2))
        r_squared = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 0.0
        A = 0.0
        return {
            "model": "logistic3",
            "x0_orig": x0_orig,
            "k": k,
            "r_squared": r_squared,
            "params": {"A": A, "L": L, "x0_norm": x0n, "x0_orig": x0_orig, "k": k},
            "y_hat": y_hat,
            "x_fit": x,
            "y_fit": y,
            "success": True,
            "note": "4p failed; used 3p",
        }

    return {
        "model": None,
        "x0_orig": np.nan,
        "k": np.nan,
        "r_squared": np.nan,
        "params": {},
        "y_hat": None,
        "x_fit": None,
        "y_fit": None,
        "success": False,
        "note": "fit failed",
    }


def fit_continuous_sigmoid_per_series(
    y,
    x=None,
    maxfev=30000,
    require_negative_k=False,
):
    """Fit 4-parameter sigmoid to continuous data for one series.

    The response is normalized to [0, 1] before fitting. Returned parameters are
    in normalized-y space, while ``x0_orig`` is in the original x-space.
    """
    if x is None:
        x = np.arange(1, len(y) + 1, dtype=float)

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    finite_mask = np.isfinite(y)

    if finite_mask.sum() < 4:
        return {
            "success": False,
            "is_valid": False,
            "note": "too_few_points",
            "params": {},
            "x0_orig": np.nan,
            "checks": {
                "x0_interior": False,
                "k_plausible": True,
                "k_negative": False,
                "trend_decreasing": False,
                "ci_finite": False,
                "asymptotes_covered": False,
            },
            "fit_trace": None,
        }

    x_fit = x[finite_mask]
    y_fit_raw = y[finite_mask]
    y_min, y_max = np.nanmin(y_fit_raw), np.nanmax(y_fit_raw)
    y_range = y_max - y_min

    if y_range < 1e-8:
        return {
            "success": False,
            "is_valid": False,
            "note": "no_variation",
            "params": {},
            "x0_orig": np.nan,
            "checks": {
                "x0_interior": False,
                "k_plausible": True,
                "k_negative": False,
                "trend_decreasing": False,
                "ci_finite": False,
                "asymptotes_covered": False,
            },
            "fit_trace": None,
        }

    y_norm = (y_fit_raw - y_min) / y_range
    p0 = [1.0, -1.0 if require_negative_k else 1.0, float(np.median(x_fit)), 0.0]
    bounds = (
        [-np.inf, -np.inf, np.min(x_fit), -np.inf],
        [np.inf, np.inf, np.max(x_fit), np.inf],
    )

    try:
        params_fit, pcov = curve_fit(
            _sigmoid_model,
            x_fit,
            y_norm,
            p0=p0,
            bounds=bounds,
            maxfev=maxfev,
        )
        L, k, x0, b = [float(v) for v in params_fit]
        y_hat_norm = _sigmoid_model(x_fit, *params_fit)
        y_hat_raw = y_hat_norm * y_range + y_min
        trend_decreasing = bool(y_hat_norm[0] > y_hat_norm[-1] + 1e-6)

        _, _, checks = _sigmoid_quality_checks(x_fit, params_fit, pcov)
        checks = dict(checks)
        checks["k_negative"] = bool(np.isfinite(k) and (k < 0))
        checks["trend_decreasing"] = trend_decreasing
        if require_negative_k and (not checks["k_negative"]):
            checks["k_plausible"] = False
        if require_negative_k and (not checks["trend_decreasing"]):
            checks["trend_decreasing"] = False

        failed = [name for name, ok in checks.items() if not ok]
        is_valid = len(failed) == 0
        note = "ok" if is_valid else ";".join(failed)

        return {
            "success": True,
            "is_valid": bool(is_valid),
            "note": note,
            "params": {"L": L, "k": k, "x0_orig": x0, "b": b},
            "x0_orig": x0,
            "checks": checks,
            "fit_trace": {
                "x_fit": x_fit,
                "y_fit_raw": y_fit_raw,
                "y_hat_raw": y_hat_raw,
                "y_norm": y_norm,
                "y_hat_norm": y_hat_norm,
            },
        }
    except Exception:
        return {
            "success": False,
            "is_valid": False,
            "note": "fit_failed",
            "params": {},
            "x0_orig": np.nan,
            "checks": {
                "x0_interior": False,
                "k_plausible": True,
                "k_negative": False,
                "trend_decreasing": False,
                "ci_finite": False,
                "asymptotes_covered": False,
            },
            "fit_trace": None,
        }


def fit_logistic_transitions_by_id(
    df,
    signal_col,
    id_col="id",
    trial_col="trial",
    value_transform=None,
    prefer_4p=True,
    direction="decreasing",
    maxfev=60000,
    min_x0=0,
    max_x0=None,
):
    """Fit logistic transitions for each subject in a dataframe."""
    all_fits = []
    fit_traces = {}

    for subject_id in sorted(df[id_col].dropna().unique()):
        sdf = df.loc[df[id_col] == subject_id].sort_values(trial_col)
        y = sdf[signal_col].to_numpy(dtype=float)
        if value_transform is not None:
            y = np.asarray(value_transform(y), dtype=float)
        x = np.arange(len(y), dtype=float)

        fit = fit_logistic_per_series(
            y,
            x=x,
            prefer_4p=prefer_4p,
            direction=direction,
            maxfev=maxfev,
        )
        all_fits.append({
            id_col: subject_id,
            **fit["params"],
            "model": fit["model"],
            "x0_orig": fit["x0_orig"],
            "k": fit["k"],
            "r_squared": fit["r_squared"],
            "success": fit["success"],
            "note": fit["note"],
        })

        if fit["success"]:
            fit_traces[subject_id] = {
                "x": fit.get("x_fit", x),
                "y": fit.get("y_fit", y),
                "y_hat": fit["y_hat"],
                "model": fit["model"],
                "x0_orig": fit["x0_orig"],
            }

    fits_df = pd.DataFrame(all_fits)
    if fits_df.empty:
        return fits_df, fit_traces

    mask = fits_df["success"] == True
    if min_x0 is not None:
        mask &= fits_df["x0_orig"] > min_x0
    if max_x0 is not None:
        mask &= fits_df["x0_orig"] < max_x0

    fits_df = fits_df.loc[mask].copy()
    valid_ids = set(fits_df[id_col].tolist())
    fit_traces = {k: v for k, v in fit_traces.items() if k in valid_ids}
    return fits_df, fit_traces


def fit_continuous_sigmoid_transitions_by_id(
    df,
    signal_col,
    id_col="id",
    trial_col="trial",
    maxfev=30000,
    require_negative_k=True,
    min_x0=0,
    max_x0=None,
):
    """Fit continuous sigmoid transitions for each subject in a dataframe."""
    all_fits = []
    fit_traces = {}

    for subject_id in sorted(df[id_col].dropna().unique()):
        sdf = df.loc[df[id_col] == subject_id].sort_values(trial_col)
        x = np.arange(1, len(sdf) + 1, dtype=float)
        y = sdf[signal_col].to_numpy(dtype=float)

        fit = fit_continuous_sigmoid_per_series(
            y,
            x=x,
            maxfev=maxfev,
            require_negative_k=require_negative_k,
        )

        params = fit.get("params", {})
        checks = fit.get("checks", {})
        all_fits.append({
            id_col: subject_id,
            "L": params.get("L", np.nan),
            "k": params.get("k", np.nan),
            "x0_orig": params.get("x0_orig", np.nan),
            "b": params.get("b", np.nan),
            "success": fit["success"],
            "is_valid": fit["is_valid"],
            "note": fit["note"],
            "x0_interior": checks.get("x0_interior", False),
            "k_plausible": checks.get("k_plausible", True),
            "k_negative": checks.get("k_negative", False),
            "trend_decreasing": checks.get("trend_decreasing", False),
            "ci_finite": checks.get("ci_finite", False),
            "asymptotes_covered": checks.get("asymptotes_covered", False),
        })

        if fit["success"] and (fit.get("fit_trace") is not None):
            fit_traces[subject_id] = fit["fit_trace"]

    fits_df = pd.DataFrame(all_fits)
    if fits_df.empty:
        return fits_df, fit_traces

    mask = (fits_df["success"] == True) & (fits_df["is_valid"] == True)
    if min_x0 is not None:
        mask &= fits_df["x0_orig"] > min_x0
    if max_x0 is not None:
        mask &= fits_df["x0_orig"] < max_x0

    fits_df = fits_df.loc[mask].copy()
    valid_ids = set(fits_df[id_col].tolist())
    fit_traces = {k: v for k, v in fit_traces.items() if k in valid_ids}
    return fits_df, fit_traces


def build_realigned_trials(df, fits_df, output_col, id_col="id", trial_col="trial"):
    """Create a trial-aligned column from a transition table containing x0_orig."""
    if fits_df is None or fits_df.empty:
        return pd.Series(np.nan, index=df.index, name=output_col, dtype=float)

    transitions = (
        fits_df[[id_col, "x0_orig"]]
        .drop_duplicates(subset=[id_col])
        .assign(_transition=lambda frame: frame["x0_orig"].round().astype(int))
        .set_index(id_col)["_transition"]
    )

    aligned = []
    for _, row in df.iterrows():
        subject_id = row[id_col]
        if subject_id not in transitions.index:
            aligned.append(np.nan)
        else:
            aligned.append(row[trial_col] - int(transitions.loc[subject_id]))

    return pd.Series(aligned, index=df.index, name=output_col, dtype=float)


def run_model_fit_comparison(fit_inputs, round_digits=4):
    """Fit linear/exponential/sigmoidal models for each condition.

    Parameters
    ----------
    fit_inputs : dict[str, np.ndarray]
        Mapping of condition label to 1D response array ordered by trial.
    round_digits : int
        Number of decimals to report in result tables.

    Returns
    -------
    tuple[pd.DataFrame, pd.DataFrame]
        fit_results table and sigmoid_diagnostics table.
    """
    rows = []
    sigmoid_diagnostics = []

    for condition_label, y in fit_inputs.items():
        y = np.asarray(y, dtype=float)
        x = np.arange(1, len(y) + 1, dtype=float)

        lin = _fit_linear(x, y)
        rows.append({
            "Condition": condition_label,
            "Model": "Linear",
            "r": lin["r"],
            "p": lin["p"],
            "RMSE": lin["rmse"],
            "AICc": lin["aicc"],
            "BIC": lin["bic"],
            "Sigmoid_valid": np.nan,
            "Sigmoid_flags": "",
        })

        exp_fit = _fit_curve("exponential", _exp_model, x, y)
        rows.append({
            "Condition": condition_label,
            "Model": "Exponential",
            "r": exp_fit["r"],
            "p": exp_fit["p"],
            "RMSE": exp_fit["rmse"],
            "AICc": exp_fit["aicc"],
            "BIC": exp_fit["bic"],
            "Sigmoid_valid": np.nan,
            "Sigmoid_flags": "",
        })

        sig_fit = _fit_curve("sigmoidal", _sigmoid_model, x, y)
        is_valid, reasons, checks = _sigmoid_quality_checks(x, sig_fit["params"], sig_fit["pcov"])
        rows.append({
            "Condition": condition_label,
            "Model": "Sigmoidal",
            "r": sig_fit["r"],
            "p": sig_fit["p"],
            "RMSE": sig_fit["rmse"],
            "AICc": sig_fit["aicc"],
            "BIC": sig_fit["bic"],
            "Sigmoid_valid": is_valid,
            "Sigmoid_flags": reasons,
        })

        L, k, x0, b = sig_fit["params"] if len(sig_fit["params"]) == 4 else [np.nan, np.nan, np.nan, np.nan]
        sigmoid_diagnostics.append({
            "Condition": condition_label,
            "L": L,
            "k": k,
            "x0": x0,
            "b": b,
            "Sigmoid_valid": is_valid,
            "Sigmoid_flags": reasons,
            "x0_interior": checks["x0_interior"],
            "k_plausible": checks["k_plausible"],
            "ci_finite": checks["ci_finite"],
            "asymptotes_covered": checks["asymptotes_covered"],
        })

    fit_results = pd.DataFrame(rows)
    fit_results["AICc_rank"] = fit_results.groupby("Condition")["AICc"].rank(method="dense")
    fit_results["Best_by_AICc"] = fit_results["AICc_rank"] == 1

    for col in ["r", "p", "RMSE", "AICc", "BIC", "AICc_rank"]:
        fit_results[col] = fit_results[col].round(round_digits)

    sigmoid_diagnostics_table = pd.DataFrame(sigmoid_diagnostics)
    if not sigmoid_diagnostics_table.empty:
        sigmoid_diagnostics_table[["L", "k", "x0", "b"]] = sigmoid_diagnostics_table[["L", "k", "x0", "b"]].round(round_digits)

    return fit_results, sigmoid_diagnostics_table


def summarize_best_valid_model(fit_results, round_digits=4):
    """Select best model by condition with conservative nonlinear gating.

    Default behavior is linear unless a nonlinear model passes all gates:
    - p-value below threshold
    - AICc improvement versus linear of at least ``min_delta_aicc``
    - (sigmoidal only) passes sigmoid quality checks
    """
    p_threshold = 0.05
    min_delta_aicc = 2.0
    summary_rows = []

    for condition_label, group in fit_results.groupby("Condition"):
        group = group.copy()
        linear_rows = group[group["Model"] == "Linear"]
        linear_aicc = linear_rows["AICc"].iloc[0] if len(linear_rows) > 0 else np.nan

        group["passes_sigmoid_quality"] = np.where(
            group["Model"] == "Sigmoidal",
            group["Sigmoid_valid"].fillna(False),
            True,
        )
        group["passes_p_gate"] = group["p"].notna() & (group["p"] < p_threshold)
        group["delta_aicc_vs_linear"] = linear_aicc - group["AICc"]
        group["passes_delta_aicc_gate"] = (
            group["AICc"].notna()
            & np.isfinite(linear_aicc)
            & (group["delta_aicc_vs_linear"] >= min_delta_aicc)
        )

        is_linear = group["Model"] == "Linear"
        is_nonlinear = group["Model"].isin(["Exponential", "Sigmoidal"])
        group["Is_valid_for_selection"] = is_linear | (
            is_nonlinear
            & group["passes_sigmoid_quality"]
            & group["passes_p_gate"]
            & group["passes_delta_aicc_gate"]
        )

        valid_group = group[group["Is_valid_for_selection"] == True]
        if len(valid_group) > 0 and valid_group["AICc"].notna().any():
            best = valid_group.sort_values("AICc", ascending=True).iloc[0]
            if best["Model"] == "Linear":
                selection_reason = "linear_default_or_nonlinear_failed_gate"
            else:
                selection_reason = "nonlinear_passed_gate_and_best_AICc"
        elif len(linear_rows) > 0:
            best = linear_rows.iloc[0]
            selection_reason = "fallback_linear_missing_AICc"
        else:
            best = group.sort_values("AICc", ascending=True).iloc[0]
            selection_reason = "fallback_no_linear_available"

        summary_rows.append({
            "Condition": condition_label,
            "Selected_Model": best["Model"],
            "Selection_Reason": selection_reason,
            "r": best["r"],
            "p": best["p"],
            "RMSE": best["RMSE"],
            "AICc": best["AICc"],
            "BIC": best["BIC"],
            "Sigmoid_valid": best["Sigmoid_valid"],
            "Sigmoid_flags": best["Sigmoid_flags"],
        })

    best_model_summary = pd.DataFrame(summary_rows)
    for col in ["r", "p", "RMSE", "AICc", "BIC"]:
        best_model_summary[col] = best_model_summary[col].round(round_digits)

    return best_model_summary
