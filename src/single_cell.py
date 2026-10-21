"""Frozen labels and original score units; inference only over independent patients."""

from pathlib import Path
import re, json, hashlib, sys
import numpy as np, pandas as pd
from scipy import stats
from .statistics import *
from .settings import DATA, OUTPUT, DEFINITIONS
from .metadata import soft

O = OUTPUT
D = DATA / "single_cell"
T = O / "single_cell"
T.mkdir(parents=True, exist_ok=True)
stage = sys.argv[1] if len(sys.argv) > 1 else "main"
p = pd.read_csv(DEFINITIONS / "program.csv")
pos = p.query("direction==1").gene.tolist()
neg = p.query("direction==-1").gene.tolist()
cp = p.query("direction==1 and common_core").gene.tolist()
cn = p.query("direction==-1 and common_core").gene.tolist()
variants = {
    "A20": (pos, neg),
    "B19_no_PKM": ([g for g in pos if g != "PKM"], neg),
    "C_PKM": (["PKM"], []),
    "D_core11": (cp, cn),
}
rows = []
sens = []
annotation = []
sampleagg = []
counts = []
mappingrows = []
coverage = []
qccor = []
for ds in ["GSE183904", "GSE163558"]:
    x = pd.read_csv(D / (ds + "_panel_expression.csv"), index_col=0)
    m = pd.read_csv(D / (ds + "_coarse_celltype_metadata.csv"), index_col=0).reindex(
        x.index
    )
    s = pd.read_csv(D / (ds + "_cell_scores.csv"), index_col=0).reindex(x.index)
    text, sm = soft(ds)
    mapdict = {}
    primary = set()
    if ds == "GSE183904":
        for line in text.splitlines():
            if line.startswith("!Series_summary = NGC"):
                v = line.split(" = ", 1)[1].split("\t")
                patient = v[0].strip()
                for j, token in enumerate(v[1:], 1):
                    if token != "-":
                        sn = token.strip().removesuffix(".csv")
                        mapdict[sn] = patient
                        if j == 1:
                            primary.add(sn)
        m["short_sample"] = m.sample_name.str.extract(r"(sample\d+)")[0]
        m["patient_id"] = m.short_sample.map(mapdict)
        m["is_primary"] = m.short_sample.isin(primary)
    else:
        # Original authors explicitly document donor mapping in PMC8858624, Results 3.1.
        donor = {
            "PT1": "Patient1",
            "Li1": "Patient1",
            "PT2": "Patient2",
            "NT1": "Patient2",
            "O1": "Patient3",
            "PT3": "Patient4",
            "Li2": "Patient4",
            "LN1": "Patient4",
            "LN2": "Patient5",
            "P1": "Patient6",
        }
        m["short_sample"] = m.sample_name.str.split("_").str[-1]
        m["patient_id"] = m.short_sample.map(donor)
        m["is_primary"] = m.short_sample.isin(["PT1", "PT2", "PT3"])
    for sample, mm in m.groupby("sample_name"):
        mappingrows.append(
            dict(
                cohort=ds,
                sample_name=sample,
                patient_id=mm.patient_id.iloc[0],
                is_primary=mm.is_primary.iloc[0],
                source=(
                    "official GEO Manuscript_ID table"
                    if ds == "GSE183904"
                    else "original author donor mapping https://pmc.ncbi.nlm.nih.gov/articles/PMC8858624/ Results 3.1"
                ),
            )
        )
    for v, (pp, nn) in variants.items():
        ap = [g for g in pp if g in x]
        an = [g for g in nn if g in x]
        missing = [g for g in pp + nn if g not in x]
        valid = bool(ap) and (bool(an) or not nn) and (v != "D_core11" or not missing)
        s[v] = (
            x[ap].mean(axis=1) - (x[an].mean(axis=1) if nn else 0) if valid else np.nan
        )
        coverage.append(
            dict(
                cohort=ds,
                program=v,
                requested=len(pp + nn),
                available=len(ap + an),
                positive_denominator=len(ap),
                negative_denominator=len(an),
                missing_genes="|".join(missing),
                score_units="original log1p raw counts module means; no UMI normalization or gene z",
            )
        )
    assert np.allclose(s.A20, s.driver_score_sc, equal_nan=True, atol=1e-12)
    s["no_RPL36AL"] = x[[g for g in pos if g in x]].mean(axis=1) - x[
        [g for g in neg if g in x and g != "RPL36AL"]
    ].mean(axis=1)
    s["proliferation"] = x[
        [g for g in ["MKI67", "TOP2A", "UBE2C", "STMN1"] if g in x]
    ].mean(axis=1)
    if stage == "sensitivity":
        qc = pd.read_csv(D / (ds + "_full_QC.csv.gz"), index_col=0).reindex(x.index)
        assert qc.total_counts.notna().all()
        s["logUMI"] = np.log1p(qc.total_counts)
        s["logGenes"] = np.log1p(qc.n_genes_detected)
        for cov in ["logUMI", "logGenes", "proliferation"]:
            for v in variants:
                s[v + "__" + cov] = np.nan
                for sample, mm in m.groupby("sample_name"):
                    idx = mm.index
                    cv = s.loc[idx, cov].to_numpy()
                    design = np.column_stack([np.ones(len(idx)), cv])
                    y = s.loc[idx, v].to_numpy()
                    s.loc[idx, v + "__" + cov] = (
                        y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
                    )
        for marker in ["COL1A1", "COL1A2", "DCN", "EPCAM", "TAGLN"]:
            s[marker + "_expression"] = x[marker]
            s[marker + "_detected"] = x[marker] > 0
        for marker in [
            "LUM",
            "PTPRC",
            "RGS5",
            "CSPG4",
            "MCAM",
            "NOTCH3",
            "ACTA2",
            "MYH11",
            "CNN1",
            "DES",
        ]:
            if marker + "_raw" in qc:
                s[marker + "_expression"] = np.log1p(qc[marker + "_raw"])
                s[marker + "_detected"] = qc[marker + "_raw"] > 0
        for sample, mm in m.groupby("sample_name"):
            ix = mm.index
            qccor.append(
                dict(
                    cohort=ds,
                    sample_name=sample,
                    n_cells=len(ix),
                    rho_UMI_detected_genes=stats.spearmanr(
                        s.loc[ix, "logUMI"], s.loc[ix, "logGenes"]
                    ).statistic,
                    interpretation="descriptive cells; no cell-level inferential P",
                )
            )
    # Frozen sample-celltype means and minimum cell count.
    merged = s.join(
        m[["sample_name", "patient_id", "is_primary", "predicted_celltype"]]
    )
    agg = merged.groupby(
        ["sample_name", "patient_id", "is_primary", "predicted_celltype"], dropna=False
    ).mean(numeric_only=True)
    agg["n_cells"] = merged.groupby(
        ["sample_name", "patient_id", "is_primary", "predicted_celltype"], dropna=False
    ).size()
    agg = agg.reset_index()
    agg["cohort"] = ds
    sampleagg.append(agg)
    counts.append(
        agg[
            [
                "cohort",
                "sample_name",
                "patient_id",
                "is_primary",
                "predicted_celltype",
                "n_cells",
            ]
        ]
    )

    def comparisons(frame, mode, cols, collection, mincells=20, primary_only=True):
        z = frame[(frame.n_cells >= mincells) & frame.patient_id.notna()]
        if primary_only:
            z = z[z.is_primary]
        for comparator in ["Epithelial", "Myeloid", "T_NK"]:
            for v, col in cols.items():
                subset = z[z.predicted_celltype.isin(["Fibroblast", comparator])]
                wide = (
                    subset.groupby(["patient_id", "predicted_celltype"])[col]
                    .mean()
                    .unstack("predicted_celltype")
                )
                if not set(["Fibroblast", comparator]) <= set(wide.columns):
                    dd = []
                else:
                    dd = (wide.Fibroblast - wide[comparator]).dropna().to_numpy()
                collection.append(
                    dict(
                        cohort=ds,
                        program=v,
                        comparator=comparator,
                        analysis=mode,
                        min_cells_per_sample_type=mincells,
                        scope=(
                            "primary tumor"
                            if primary_only
                            else "all tissues; equal within-patient sample weighting"
                        ),
                        **paired(dd)
                    )
                )

    comparisons(agg, "PRIMARY_FROZEN_SCORE", {v: v for v in variants}, rows)
    if stage == "main":
        comparisons(
            agg,
            "INCLUDE_SMALL_GROUPS_AT_LEAST1",
            {v: v for v in variants},
            sens,
            mincells=1,
        )
        comparisons(
            agg,
            "ALL_TISSUES_PATIENT_AGGREGATED",
            {v: v for v in variants},
            sens,
            primary_only=False,
        )
    if stage == "sensitivity":
        for cov in ["logUMI", "logGenes", "proliferation"]:
            comparisons(
                agg,
                "SINGLE_COVARIATE_RESIDUAL_" + cov,
                {v: v + "__" + cov for v in variants},
                sens,
            )
        comparisons(
            agg,
            "DELETE_RPL36AL_CONTROL",
            {"A20_no_RPL36AL_sensitivity": "no_RPL36AL"},
            sens,
        )
        comparisons(
            agg,
            "QC_MEAN_DIFFERENCES",
            {v: v for v in ["logUMI", "logGenes", "proliferation"]},
            sens,
        )
        for comparator in ["Epithelial", "Myeloid", "T_NK"]:
            balanced = []
            rng = np.random.default_rng(SEED)
            for sample, mm in m[m.is_primary].groupby("sample_name"):
                groups = {
                    ct: mm.index[mm.predicted_celltype == ct].to_numpy()
                    for ct in ["Fibroblast", comparator]
                }
                n = min(map(len, groups.values()))
                if n < 20:
                    continue
                for ct, ids in groups.items():
                    pick = rng.choice(ids, n, replace=False)
                    rec = {v: s.loc[pick, v].mean() for v in variants}
                    rec.update(
                        patient_id=mm.patient_id.iloc[0],
                        sample_name=sample,
                        predicted_celltype=ct,
                        n_cells=n,
                        is_primary=True,
                    )
                    balanced.append(rec)
            bal = pd.DataFrame(balanced)
            if len(bal):
                for v in variants:
                    wide = (
                        bal.groupby(["patient_id", "predicted_celltype"])[v]
                        .mean()
                        .unstack()
                    )
                    dd = (wide.Fibroblast - wide[comparator]).dropna()
                    sens.append(
                        dict(
                            cohort=ds,
                            program=v,
                            comparator=comparator,
                            analysis="EQUAL_CELL_COUNT_SINGLE_FIXED_DRAW",
                            min_cells_per_sample_type=20,
                            scope="primary tumor; patient units",
                            **paired(dd)
                        )
                    )
        # Audit-marker patient mean differences with held-out vs classifier provenance explicitly marked.
        for col in [c for c in s if c.endswith(("_expression", "_detected"))]:
            temp = []
            comparisons(agg, "MARKER_LABEL_AUDIT", {col: col}, temp)
            annotation.extend(temp)
pd.concat(counts).to_csv(T / "02_cell_counts_by_patient_type.csv", index=False)
pd.DataFrame(mappingrows).to_csv(T / "02_sample_patient_mapping.csv", index=False)
pd.DataFrame(coverage).to_csv(T / "02_program_gene_coverage.csv", index=False)
pd.concat(sampleagg).to_csv(
    T / ("02_patient_celltype_summaries_" + stage + ".csv"), index=False
)
if stage == "main":
    fdr(pd.DataFrame(rows), ["cohort"]).to_csv(
        T / "02_sample_level_CAF_comparison.csv", index=False
    )
    fdr(pd.DataFrame(sens), ["cohort", "analysis"]).to_csv(
        T / "02_CAF_initial_sensitivity.csv", index=False
    )
else:
    initial = pd.read_csv(T / "02_CAF_initial_sensitivity.csv")
    a = fdr(pd.DataFrame(sens), ["cohort", "analysis"])
    pd.concat([initial, a], ignore_index=True).to_csv(
        T / "02_CAF_specificity_sensitivity.csv", index=False
    )
    fdr(pd.DataFrame(annotation), ["cohort"]).to_csv(
        T / "02_marker_patient_comparisons.csv", index=False
    )
    pd.DataFrame(qccor).to_csv(T / "02_QC_collinearity.csv", index=False)
print(pd.DataFrame(rows).to_string(index=False))
print("SC_DONE", stage)
