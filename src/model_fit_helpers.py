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


def _sigmoid_quality_checks(x, params, pcov):
    if params is None or len(params) != 4 or np.any(~np.isfinite(params)):
        return False, "fit_failed", {"x0_interior": False, "k_plausible": False, "ci_finite": False, "asymptotes_covered": False}

    L, k, x0, b = params
    x_min, x_max = float(np.min(x)), float(np.max(x))
    x_range = max(x_max - x_min, 1.0)
    edge_margin = 0.15 * x_range
    x0_interior = (x_min + edge_margin) <= x0 <= (x_max - edge_margin)
    k_plausible = np.isfinite(k) and (0.02 <= abs(k) <= 2.5)

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
