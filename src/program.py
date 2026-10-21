"""Historical scoring helpers. The manuscript gene set is fixed in data/definitions."""

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


def bh_fdr(pvals):
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    if n == 0:
        return np.array([])
    order = np.argsort(pvals)
    ranked = pvals[order]
    fdr = ranked * n / (np.arange(1, n + 1))
    fdr = np.minimum.accumulate(fdr[::-1])[::-1]
    fdr = np.clip(fdr, 0, 1)
    out = np.empty(n, dtype=float)
    out[order] = fdr
    return out


def safe_mwu(x, y):
    x = pd.to_numeric(pd.Series(x), errors="coerce").dropna().values
    y = pd.to_numeric(pd.Series(y), errors="coerce").dropna().values
    if len(x) < 3 or len(y) < 3:
        return np.nan, np.nan
    try:
        stat, p = mannwhitneyu(x, y, alternative="two-sided")
        return stat, p
    except Exception:
        return np.nan, np.nan


def diff_table_numeric(df, labels, group_a, group_b, min_abs_effect=0.0):
    idx_a = labels[labels == group_a].index
    idx_b = labels[labels == group_b].index

    rows = []
    for col in df.columns:
        xa = pd.to_numeric(df.loc[idx_a, col], errors="coerce")
        xb = pd.to_numeric(df.loc[idx_b, col], errors="coerce")

        _, p = safe_mwu(xa, xb)
        mean_a = xa.mean()
        mean_b = xb.mean()
        med_a = xa.median()
        med_b = xb.median()
        effect = mean_a - mean_b
        rows.append(
            {
                "feature": col,
                f"mean_{group_a}": mean_a,
                f"mean_{group_b}": mean_b,
                f"median_{group_a}": med_a,
                f"median_{group_b}": med_b,
                "effect_mean_diff": effect,
                "abs_effect": abs(effect),
                "pvalue": p,
            }
        )

    res = pd.DataFrame(rows)
    res = res.dropna(subset=["pvalue"]).copy()
    if len(res) == 0:
        return res
    res["fdr"] = bh_fdr(res["pvalue"].values)
    res = res.sort_values(["fdr", "abs_effect"], ascending=[True, False]).reset_index(
        drop=True
    )

    if min_abs_effect > 0:
        res = res[res["abs_effect"] >= min_abs_effect].copy()

    return res


def zscore_by_dataset(df):
    df = df.copy()
    mu = df.mean(axis=0)
    sd = df.std(axis=0, ddof=0)
    sd[sd == 0] = 1.0
    sd = sd.fillna(1.0)
    out = (df - mu) / sd
    return out, mu, sd


def build_driver_score(expr_z, up_genes, down_genes):
    up_genes = [g for g in up_genes if g in expr_z.columns]
    down_genes = [g for g in down_genes if g in expr_z.columns]

    if len(up_genes) == 0 and len(down_genes) == 0:
        return pd.Series(np.nan, index=expr_z.index)

    up_score = expr_z[up_genes].mean(axis=1) if len(up_genes) > 0 else 0.0
    down_score = expr_z[down_genes].mean(axis=1) if len(down_genes) > 0 else 0.0
    score = up_score - down_score
    return score
