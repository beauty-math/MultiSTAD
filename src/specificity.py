"""Prespecified paired comparison to the existing no-PKM 7-gene glycolysis comparator."""

from pathlib import Path
import pandas as pd, numpy as np
from scipy import stats
from .statistics import *
from .settings import OUTPUT

T = OUTPUT / "bulk"
rows = []
for num, cohort in enumerate(["GSE62254", "GSE15459", "GSE84433", "GSE26253", "CPTAC"]):
    s = pd.read_csv(T / (cohort + "_phase2_patient_scores.csv"), index_col=0)
    for target in [
        "Hallmark_Hypoxia",
        "IFNG_signature",
        "T_cell_inflamed",
        "CAF_composite",
    ]:
        sub = s[["A20", "glycolysis7_no_PKM", target + "__nonoverlap"]].dropna()
        cm = corr_matrix(sub.to_numpy())
        bm = boot_matrix(sub.to_numpy(), SEED + num)
        d = cm[0, 2] - cm[1, 2]
        b = bm[:, 0, 2] - bm[:, 1, 2]
        rows.append(
            dict(
                cohort=cohort,
                target=target,
                n=len(sub),
                rho_program=cm[0, 2],
                program_ci_low=ci(bm[:, 0, 2])[0],
                program_ci_high=ci(bm[:, 0, 2])[1],
                rho_simple_glycolysis=cm[1, 2],
                simple_ci_low=ci(bm[:, 1, 2])[0],
                simple_ci_high=ci(bm[:, 1, 2])[1],
                delta_rho_program_minus_simple=d,
                ci_low=ci(b)[0],
                ci_high=ci(b)[1],
                p_value=diff_p(b, d),
                method="patient paired bootstrap B5000; centered null approximation",
            )
        )
fdr(pd.DataFrame(rows), ["cohort"]).to_csv(
    T / "01_simple_glycolysis_comparison.csv", index=False
)
