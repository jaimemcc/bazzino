"""Reusable statistical helpers for figure notebooks."""

from itertools import combinations

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import mannwhitneyu, wilcoxon
from statsmodels.stats.multitest import multipletests


def wilcoxon_vs_zero(values):
    """Run one-sample Wilcoxon signed-rank test against zero with safeguards."""
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    n = len(vals)

    out = {
        "n": n,
        "Mean": np.nan if n == 0 else float(np.mean(vals)),
        "Median": np.nan if n == 0 else float(np.median(vals)),
        "SD": np.nan if n < 2 else float(np.std(vals, ddof=1)),
        "Test": "wilcoxon signed-rank",
        "Statistic": np.nan,
        "p-value": np.nan,
        "Note": "",
    }

    if n == 0:
        out["Note"] = "no finite observations"
        return out

    if np.allclose(vals, 0):
        out["Statistic"] = 0.0
        out["p-value"] = 1.0
        out["Note"] = "all values are zero"
        return out

    try:
        w_res = wilcoxon(vals, zero_method="wilcox", alternative="two-sided")
        out["Statistic"] = float(w_res.statistic)
        out["p-value"] = float(w_res.pvalue)
    except ValueError as e:
        out["Note"] = str(e)

    return out


def one_sample_wilcoxon_table(group_data, round_digits=4):
    """Return one-sample Wilcoxon summary table for a dict of named arrays."""
    rows = []
    for label, vals in group_data.items():
        rows.append({"Group": label, **wilcoxon_vs_zero(vals)})

    out = pd.DataFrame(rows)
    for col in ["Mean", "Median", "SD", "Statistic", "p-value"]:
        out[col] = out[col].astype(float).round(round_digits)
    return out


def prepare_mixedlm_auc_df(x_array, value_col="auc_snips"):
    """Aggregate trial-level values to subject x condition x infusiontype x sex.

    Parameters
    ----------
    x_array : pd.DataFrame
        Source dataframe containing id, condition, infusiontype, sex, and value_col.
    value_col : str, default "auc_snips"
        Column to aggregate as the response variable.
    """
    model_df = (
        x_array[["id", "condition", "infusiontype", "sex", value_col]]
        .dropna(subset=["id", "condition", "infusiontype", "sex", value_col])
        .groupby(["id", "condition", "infusiontype", "sex"], as_index=False)[value_col]
        .mean()
        .rename(columns={value_col: "auc"})
    )

    n_cond_by_id = model_df.groupby("id")["condition"].nunique()
    valid_ids = n_cond_by_id[n_cond_by_id > 1].index
    model_df = model_df[model_df["id"].isin(valid_ids)].copy()

    model_df["condition"] = pd.Categorical(
        model_df["condition"],
        categories=["replete", "deplete"],
        ordered=True,
    )
    model_df["infusiontype"] = pd.Categorical(model_df["infusiontype"])
    model_df["sex"] = pd.Categorical(model_df["sex"])

    return model_df


def fit_mixedlm_with_fallback(model_df, candidate_formulas=None, group_col="id"):
    """Fit a mixed model, falling back to simpler formulas if needed."""
    if candidate_formulas is None:
        candidate_formulas = [
            "auc ~ C(condition) * C(infusiontype) * C(sex)",
            "auc ~ C(condition) * C(infusiontype) + C(sex)",
            "auc ~ C(condition) * C(infusiontype)",
            "auc ~ C(condition) + C(infusiontype) + C(sex)",
            "auc ~ C(condition) + C(infusiontype)",
        ]

    _methods = ["lbfgs", "nm", "powell"]
    fit_errors = []
    for formula in candidate_formulas:
        for method in _methods:
            try:
                model = smf.mixedlm(formula, data=model_df, groups=model_df[group_col])
                res = model.fit(reml=False, method=method)
                return res, formula, fit_errors
            except Exception as e:
                fit_errors.append((formula, method, str(e)))

    return None, None, fit_errors


def omnibus_wald_table(mixed_result):
    """Return Wald table for fixed effects from a fitted MixedLM result."""
    wald = mixed_result.wald_test_terms()
    return wald.table.reset_index().rename(columns={"index": "Term"})


def posthoc_pairwise_auc(x_array, value_col="auc_snips", round_digits=4):
    """Run pairwise posthoc comparisons across condition x infusiontype groups.

    Uses paired Wilcoxon within infusion contrasts and Mann-Whitney U for
    across-infusion contrasts; applies Holm correction to valid p-values.

    Parameters
    ----------
    x_array : pd.DataFrame
        Source dataframe containing id, condition, infusiontype, and value_col.
    value_col : str, default "auc_snips"
        Column to aggregate as the response variable.
    round_digits : int, default 4
        Decimal places for reported statistics.
    """
    posthoc_df = (
        x_array[["id", "condition", "infusiontype", value_col]]
        .dropna(subset=["id", "condition", "infusiontype", value_col])
        .groupby(["id", "condition", "infusiontype"], as_index=False)[value_col]
        .mean()
        .rename(columns={value_col: "auc"})
    )

    posthoc_df["group"] = (
        posthoc_df["condition"].astype(str)
        + " | "
        + posthoc_df["infusiontype"].astype(str)
    )
    groups = sorted(posthoc_df["group"].unique())

    rows = []
    for g1, g2 in combinations(groups, 2):
        d1 = (
            posthoc_df.loc[posthoc_df["group"] == g1, ["id", "auc"]]
            .rename(columns={"auc": "auc_1"})
        )
        d2 = (
            posthoc_df.loc[posthoc_df["group"] == g2, ["id", "auc"]]
            .rename(columns={"auc": "auc_2"})
        )

        inf1 = g1.split(" | ")[1]
        inf2 = g2.split(" | ")[1]

        if inf1 == inf2:
            paired = d1.merge(d2, on="id", how="inner")
            n_pairs = len(paired)
            if n_pairs > 0:
                if np.allclose(paired["auc_1"].to_numpy(), paired["auc_2"].to_numpy()):
                    stat = 0.0
                    p_raw = 1.0
                    note = "all paired differences are zero"
                else:
                    w = wilcoxon(
                        paired["auc_1"],
                        paired["auc_2"],
                        alternative="two-sided",
                        zero_method="wilcox",
                    )
                    stat = float(w.statistic)
                    p_raw = float(w.pvalue)
                    note = ""
                test_name = "Wilcoxon signed-rank (paired)"
            else:
                stat = np.nan
                p_raw = np.nan
                note = "no paired observations"
                test_name = "Wilcoxon signed-rank (paired)"
            n1, n2 = len(d1), len(d2)
        else:
            v1 = d1["auc_1"].to_numpy(dtype=float)
            v2 = d2["auc_2"].to_numpy(dtype=float)
            n1, n2 = len(v1), len(v2)
            n_pairs = np.nan
            if n1 > 0 and n2 > 0:
                u = mannwhitneyu(v1, v2, alternative="two-sided")
                stat = float(u.statistic)
                p_raw = float(u.pvalue)
                note = ""
            else:
                stat = np.nan
                p_raw = np.nan
                note = "one or both groups empty"
            test_name = "Mann-Whitney U (unpaired)"

        rows.append(
            {
                "Comparison": f"{g1} vs {g2}",
                "Test": test_name,
                "n1": n1,
                "n2": n2,
                "n_pairs": n_pairs,
                "Statistic": stat,
                "p_raw": p_raw,
                "Note": note,
            }
        )

    results = pd.DataFrame(rows)
    valid_mask = results["p_raw"].notna()
    results["p_holm"] = np.nan
    results["reject_holm_0.05"] = False
    if valid_mask.any():
        rej, p_adj, _, _ = multipletests(
            results.loc[valid_mask, "p_raw"],
            alpha=0.05,
            method="holm",
        )
        results.loc[valid_mask, "p_holm"] = p_adj
        results.loc[valid_mask, "reject_holm_0.05"] = rej

    for col in ["Statistic", "p_raw", "p_holm"]:
        results[col] = results[col].astype(float).round(round_digits)

    return results, groups


def holm_matrix_from_posthoc(results_df, groups, round_digits=4):
    """Build symmetric Holm-adjusted p-value matrix from posthoc result table."""
    p_matrix = pd.DataFrame(np.nan, index=groups, columns=groups)
    for g in groups:
        p_matrix.loc[g, g] = 0.0

    for _, row in results_df.iterrows():
        if pd.notna(row["p_holm"]):
            g1, g2 = row["Comparison"].split(" vs ")
            p_matrix.loc[g1, g2] = row["p_holm"]
            p_matrix.loc[g2, g1] = row["p_holm"]

    return p_matrix.round(round_digits)
