"""Exploratory immunotherapy association using the frozen score and five comparators.

The DeLong influence calculation is retained from the current evidence analysis.
Bootstrap intervals are read from the archived comparison table, never invented.
"""

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

COMPARATORS = {
    "PKM": "PKM",
    "IFNG": "IFNG",
    "T_cell_7gene": "Tcell_inflamed",
    "CD8": "CD8",
    "Merck18": "Merck18",
}


def evaluate_tumors(scores_path, archived_comparisons, output):
    x = pd.read_csv(scores_path)
    if len(x) != 45 or x.patient_id.nunique() != 45:
        raise ValueError("Expected 45 unique tumor patients")
    y = x.response_binary.to_numpy()
    base = x.ecosystem_score.to_numpy()
    if set(y) != {0, 1} or y.sum() != 12:
        raise ValueError("Expected 12 responders and 33 nonresponders")
    if not np.isfinite(x[["ecosystem_score", *COMPARATORS.values()]].to_numpy()).all():
        raise ValueError("Missing/nonfinite frozen scores")
    if (
        "normal_suffix" in x
        and x.normal_suffix.astype(str).str.lower().isin(["true", "1"]).any()
    ):
        raise ValueError("Normal tissue records must be excluded")

    def influence(score):
        a = score[y == 1]
        b = score[y == 0]
        h = (a[:, None] > b).astype(float) + 0.5 * (a[:, None] == b)
        return h.mean(), h.mean(axis=1), h.mean(axis=0)

    auc, va, vb = influence(base)
    rows = []
    archived = pd.read_csv(archived_comparisons) if archived_comparisons else None
    for name, column in COMPARATORS.items():
        other, vc, vd = influence(x[column].to_numpy())
        delta = auc - other
        variance = np.var(va - vc, ddof=1) / len(va) + np.var(vb - vd, ddof=1) / len(vb)
        p = 2 * stats.norm.sf(abs(delta) / np.sqrt(variance))
        lo = hi = np.nan
        if archived is not None:
            hit = archived[archived.program == "Ecosystem_minus_" + name]
            if (
                len(hit) != 1
                or not np.isclose(hit.estimate.iloc[0], delta, atol=1e-10, rtol=0)
                or not np.isclose(hit.p_value.iloc[0], p, atol=1e-10, rtol=0)
            ):
                raise ValueError(
                    "Archived comparison conflicts with frozen scores: " + name
                )
            lo, hi = hit.ci_low.iloc[0], hit.ci_high.iloc[0]
        rows.append(
            dict(
                comparator=name,
                n=45,
                AUC_program=auc,
                AUC_comparator=other,
                delta_AUC=delta,
                ci_low=lo,
                ci_high=hi,
                paired_DeLong_P=p,
                CI_source="archived paired bootstrap" if archived is not None else "NA",
            )
        )
    comparisons = pd.DataFrame(rows)
    comparisons["Holm_q"] = multipletests(comparisons.paired_DeLong_P, method="holm")[1]
    z = (base - base.mean()) / base.std(ddof=0)
    fit = sm.Logit(y, sm.add_constant(z)).fit(disp=False)
    interval = np.exp(fit.conf_int()[1])
    result = dict(
        n=45,
        responders=12,
        nonresponders=33,
        AUC=auc,
        OR_per_SD=np.exp(fit.params[1]),
        OR_CI_low=interval[0],
        OR_CI_high=interval[1],
        OR_P=fit.pvalues[1],
    )
    output.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([result]).to_csv(output / "association.csv", index=False)
    comparisons.to_csv(output / "paired_comparisons.csv", index=False)
    return result, comparisons
