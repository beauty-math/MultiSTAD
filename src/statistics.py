import numpy as np, pandas as pd
from scipy import stats

SEED = 20261011
B = 5000


def bh(p):
    p = np.asarray(p, float)
    q = np.full(len(p), np.nan)
    ix = np.flatnonzero(np.isfinite(p))
    s = ix[np.argsort(p[ix])]
    q[s] = np.minimum(
        1, np.minimum.accumulate((p[s] * len(s) / np.arange(1, len(s) + 1))[::-1])[::-1]
    )
    return q


def fdr(df, groups, p="p_value", q="q_BH"):
    df[q] = np.nan
    for _, g in df.groupby(groups, dropna=False):
        df.loc[g.index, q] = bh(g[p])
    return df


def corr_matrix(a):
    r = stats.rankdata(a, axis=0)
    r -= r.mean(axis=0)
    return (r.T @ r) / np.sqrt(np.outer((r * r).sum(axis=0), (r * r).sum(axis=0)))


def boot_matrix(a, seed=SEED):
    rng = np.random.default_rng(seed)
    n, k = a.shape
    out = np.empty((B, k, k))
    for j in range(0, B, 100):
        ix = rng.integers(n, size=(min(100, B - j), n))
        r = stats.rankdata(a[ix], axis=1)
        r -= r.mean(axis=1, keepdims=True)
        cov = np.einsum("bnk,bnl->bkl", r, r)
        sd = np.sqrt(np.einsum("bnk,bnk->bk", r, r))
        out[j : j + len(ix)] = cov / (sd[:, :, None] * sd[:, None, :])
    return out


def ci(x):
    return np.nanquantile(x, [0.025, 0.975]).tolist()


def diff_p(draw, effect):
    return (1 + np.sum(np.abs(draw - effect) >= abs(effect))) / (1 + len(draw))


def partial(cm, x, y, c):
    return (cm[..., x, y] - cm[..., x, c] * cm[..., y, c]) / np.sqrt(
        (1 - cm[..., x, c] ** 2) * (1 - cm[..., y, c] ** 2)
    )


def paired(d):
    d = np.asarray(d)
    d = d[np.isfinite(d)]
    n = len(d)
    if not n:
        return dict(
            n_patients=0,
            median_paired_difference=np.nan,
            ci_low=np.nan,
            ci_high=np.nan,
            p_value=np.nan,
            status="NOT_EVALUABLE",
        )
    draw = np.median(d[np.random.default_rng(SEED).integers(n, size=(B, n))], axis=1)
    p = (
        1.0
        if np.all(d == 0)
        else stats.wilcoxon(d, alternative="two-sided", method="auto").pvalue
    )
    return dict(
        n_patients=n,
        median_paired_difference=np.median(d),
        ci_low=ci(draw)[0],
        ci_high=ci(draw)[1],
        p_value=p,
        n_positive=int((d > 0).sum()),
        status="EXPLORATORY_SMALL_N" if n < 5 else "EVALUABLE",
    )
