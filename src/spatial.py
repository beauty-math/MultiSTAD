"""Conditional-gate minimal review of existing section effects only. No new spot model."""

from pathlib import Path
import json
import pandas as pd, numpy as np
from .statistics import *
from .settings import DATA, OUTPUT
from .metadata import soft

O = OUTPUT
T = O / "spatial"
T.mkdir(parents=True, exist_ok=True)
D = DATA / "spatial"
gate = json.loads((DATA / "spatial/gate.json").read_text())
assert gate["decision"] in ["GO", "CONDITIONAL GO"]
text, sm = soft("GSE251950")
assert "9 primary gastric cancers" in text
raw = pd.read_csv(D / "sensitivity/PKM_spatial_sensitivity_sample_effects.csv")
raw = raw[raw["sample"] != "GC6-PM"].copy()
assert raw["sample"].nunique() == 9
qc = pd.read_csv(D / "GSE251950_sample_QC_summary.csv")
q = qc[qc.is_primary].copy()
summaries = []
for sample in q["sample"]:
    x = pd.read_csv(D / "samples" / f"{sample}_spot_scores.csv.gz")
    assert len(x) == int(q.loc[q["sample"] == sample, "n_spots"].iloc[0])
    high = x.PKM_high_CAF.astype(bool)
    row = dict(
        sample=sample,
        n_spots=len(x),
        n_patients=1,
        patient_identity="GC1-GC9 distinct primary cancers per GEO; one paired metastasis GC6-PM excluded",
        n_high=int(high.sum()),
        n_reference=int((~high).sum()),
        CAF_mapped_mean=float(x.Fibroblast_fraction.mean()),
        CAF_mapped_high=float(x.loc[high, "Fibroblast_fraction"].mean()),
        CAF_mapped_reference=float(x.loc[~high, "Fibroblast_fraction"].mean()),
        PKM_delta_z=float(
            x.loc[high, "PKM_expression_z"].mean()
            - x.loc[~high, "PKM_expression_z"].mean()
        ),
        median_total_counts=float(x.total_counts.median()),
        median_selected_genes_detected=float(x.n_genes_detected.median()),
        section_ci_low=np.nan,
        section_ci_high=np.nan,
        section_inference="one specimen per patient; biological CI not identifiable; spot inference not substituted",
    )
    for outcome in [
        "Glycolysis_score_no_PKM",
        "Hypoxia_score",
        "Tcell_inflamed",
        "IFNG_score",
        "Immune_fraction",
    ]:
        row[outcome + "_delta_z"] = float(
            raw[
                (raw["sample"] == sample)
                & (raw.qc_mode == "all_tissue_spots")
                & (raw.definition == "composite_tangram_all_spots")
                & (raw.analysis == "direct")
                & (raw.outcome == outcome)
            ].delta_mean_z.iloc[0]
        )
    summaries.append(row)
pd.DataFrame(summaries).to_csv(T / "03_spatial_section_effects.csv", index=False)
rows = []
for keys, g in raw.groupby(["qc_mode", "definition", "analysis", "outcome"]):
    result = paired(g.delta_mean_z.to_numpy())
    old = pd.read_csv(D / "sensitivity/PKM_spatial_sensitivity_summary.csv")
    match = old[
        (old.qc_mode == keys[0])
        & (old.definition == keys[1])
        & (old.analysis == keys[2])
        & (old.outcome == keys[3])
    ].iloc[0]
    rows.append(
        dict(
            qc_mode=keys[0],
            definition=keys[1],
            analysis=keys[2],
            outcome=keys[3],
            mean_section_effect=g.delta_mean_z.mean(),
            n_spots_after_qc_total=int(g.n_spots_after_qc.sum()),
            historical_FDR=match.BH_FDR_within_definition,
            historical_p_value=match.wilcoxon_p,
            original_effect_reused=True,
            **result,
        )
    )
res = fdr(pd.DataFrame(rows), ["qc_mode"])
res.to_csv(T / "03_spatial_specificity.csv", index=False)
raw["ci_low"] = np.nan
raw["ci_high"] = np.nan
raw["CI_status"] = (
    "not identifiable for one biological specimen; descriptive section effect"
)
raw.to_csv(T / "03_existing_all_section_effects.csv", index=False)
main = res[
    (res.qc_mode == "all_tissue_spots")
    & (res.analysis == "direct")
    & (
        res.definition.isin(
            ["composite_tangram_all_spots", "PKM_within_tangram_fibroblast_enriched"]
        )
    )
]
