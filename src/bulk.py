"""Frozen Phase 2 bulk specificity analysis. No outcome-based selection."""

from pathlib import Path
import json, hashlib, sys
import numpy as np, pandas as pd
from scipy import stats
from .statistics import *
from .settings import DATA, OUTPUT, DEFINITIONS

O = OUTPUT
D = DATA / "bulk"
T = O / "bulk"
T.mkdir(parents=True, exist_ok=True)
(O / "logs").mkdir(exist_ok=True)
program = pd.read_csv(DEFINITIONS / "program.csv")
pos = program.query("direction==1").gene.tolist()
neg = program.query("direction==-1").gene.tolist()
cp = program.query("direction==1 and common_core").gene.tolist()
cn = program.query("direction==-1 and common_core").gene.tolist()
variants = {
    "A20": (pos, neg),
    "B19_no_PKM": ([g for g in pos if g != "PKM"], neg),
    "C_PKM": (["PKM"], []),
    "D_core11": (cp, cn),
    "E_core10_no_PKM": ([g for g in cp if g != "PKM"], cn),
}
manifest = json.loads((DEFINITIONS / "signatures.json").read_text())
targets = {
    k: manifest["fixed_signatures"][k]
    for k in [
        "Hallmark_Glycolysis",
        "Hallmark_Hypoxia",
        "IFNG_signature",
        "T_cell_inflamed",
        "CAF_composite",
    ]
}
gmtfile = DATA / "external/signature/h.all.v2026.1.Hs.symbols.gmt"
gmt = {
    v[0]: v[2:]
    for line in gmtfile.read_text().splitlines()
    if len(v := line.split("\t")) >= 3
}
data = {
    c: pd.read_csv(D / (c + "_expression.csv.gz"), index_col=0)
    for c in ["GSE62254", "GSE15459", "GSE84433", "GSE26253"]
}
data["CPTAC"] = pd.read_csv(DATA / "bulk/CPTAC_RNA.cct", sep="\t", index_col=0)
meta = pd.read_csv(DATA / "bulk/CPTAC_meta.tsv", sep="\t", index_col=0).drop(
    index="CAT"
)
matching = pd.read_csv(DATA / "bulk/CPTAC_matching.csv").set_index("sample_id")
assert matching.reindex(data["CPTAC"].columns).case_id.nunique() == 151
disc = pd.read_csv(DATA / "discovery/expression.csv", index_col=0)
split = pd.read_csv(DATA / "discovery/split.csv").set_index("case_id")
disc = disc.loc[split.index[split["split"] == "train"]]
discdet = (disc > 0).mean()
# Freeze all 200 random sets per cohort before any target correlation is calculated.
randomsets = []
coverage = []
overlap = []
scoresets = {}
rng = np.random.default_rng(SEED)


def score(z, genes):
    return z.loc[[g for g in genes if g in z.index]].mean(axis=0)


for cohort, x in data.items():
    assert x.columns.is_unique and x.index.is_unique
    x = x.apply(pd.to_numeric, errors="coerce")
    z = x.sub(x.mean(axis=1), axis=0).div(
        x.std(axis=1, ddof=0).replace(0, np.nan), axis=0
    )
    s = pd.DataFrame(index=x.columns)
    for v, (pp, nn) in variants.items():
        ap = [g for g in pp if g in z.index]
        an = [g for g in nn if g in z.index]
        missing = [g for g in pp + nn if g not in z.index]
        valid = (
            len(ap) > 0
            and (len(an) > 0 or not nn)
            and (not v.startswith(("D_", "E_")) or not missing)
        )
        s[v] = score(z, pp) - (score(z, nn) if nn else 0) if valid else np.nan
        coverage.append(
            dict(
                cohort=cohort,
                program=v,
                n_expression_records=x.shape[1],
                n_independent_patients=151 if cohort == "CPTAC" else x.shape[1],
                patient_identity_status=(
                    "verified unique CPTAC Case"
                    if cohort == "CPTAC"
                    else "unique GSM; one-record-per-patient assumption not cross-study deduplicated"
                ),
                requested=len(pp + nn),
                available=len(ap + an),
                positive_denominator=len(ap),
                negative_denominator=len(an),
                available_genes="|".join(ap + an),
                missing_genes="|".join(missing),
                status="EVALUABLE_AVAILABLE_GENE" if valid else "NOT_EVALUABLE",
                score_sd=float(s[v].std(ddof=0)),
                rule="within cohort gene z ddof0, mean available positive minus mean available negative; PKM weight unchanged within module",
            )
        )
    for target, gg in targets.items():
        shared = sorted(set(gg) & set(pos + neg))
        ng = [g for g in gg if g not in shared]
        s[target] = score(z, gg)
        s[target + "__nonoverlap"] = score(z, ng)
        overlap.append(
            dict(
                cohort=cohort,
                target=target,
                shared_genes="|".join(shared),
                requested_before=len(gg),
                requested_after=len(ng),
                available_before=sum(g in z.index for g in gg),
                available_after=sum(g in z.index for g in ng),
                rule="remove every intersection with frozen20 once",
            )
        )
    gly8 = ["HK2", "PFKP", "ALDOA", "GAPDH", "ENO1", "PKM", "LDHA", "SLC2A1"]
    covdefs = {
        "glycolysis7_no_PKM": [g for g in gly8 if g != "PKM"],
        "E2F": gmt["HALLMARK_E2F_TARGETS"],
        "G2M": gmt["HALLMARK_G2M_CHECKPOINT"],
    }
    for c, gg in covdefs.items():
        s[c] = score(z, gg)
    s["global_expression"] = x.mean()
    s["detected_fraction"] = (x > 0).mean() if cohort == "CPTAC" else x.notna().mean()
    if cohort == "CPTAC":
        for c, mc in [
            ("purity_WGS", "Absolute_TumorPurity_WGS"),
            ("stromal_score", "StromalScore"),
            ("CAF_xCell", "Cancer.associated.fibroblast.xCell"),
        ]:
            s[c] = pd.to_numeric(meta[mc], errors="coerce").reindex(s.index)
    scoresets[cohort] = s
    s.to_csv(T / (cohort + "_phase2_patient_scores.csv"))
    pool = x.index[
        (x.std(axis=1, ddof=0) > 0)
        & x.notna().all(axis=1)
        & x.index.isin(discdet.index)
        & ~x.index.isin(pos + neg)
    ].tolist()
    feat = pd.DataFrame(
        {
            "mean": x.mean(axis=1).rank(pct=True),
            "detection": (
                (x > 0).mean(axis=1) if cohort == "CPTAC" else x.notna().mean(axis=1)
            ).rank(pct=True),
            "variance": x.var(axis=1).rank(pct=True),
            "discovery_detection": discdet.reindex(x.index).rank(pct=True),
        }
    ).fillna(0)
    ap = [g for g in pos if g in x.index]
    an = [g for g in neg if g in x.index]
    neighbours = {
        g: np.array(pool)[
            np.argsort(
                ((feat.loc[pool].to_numpy() - feat.loc[g].to_numpy()) ** 2).sum(axis=1)
            )[:50]
        ]
        for g in ap + an
    }
    for rep in range(200):
        used = set()
        donors = []
        for g in ap + an:
            allowed = [h for h in neighbours[g] if h not in used]
            h = str(rng.choice(allowed))
            used.add(h)
            donors.append(h)
        randomsets.append(
            dict(
                cohort=cohort,
                replicate=rep,
                positive="|".join(donors[: len(ap)]),
                negative="|".join(donors[len(ap) :]),
                source_gene_order="|".join(ap + an),
                seed=SEED,
                match_rule="uniform nearest50 by mean/detection/variance percentile plus Discovery detection; excludes fixed20; without replacement within set",
                detection_limit="GEO finite measurement fraction is not biological detection; Discovery nonzero detection added",
            )
        )
frozen = pd.DataFrame(randomsets)
frozen.to_csv(T / "random_sets_frozen_before_associations.csv", index=False)
pd.DataFrame(coverage).to_csv(T / "01_program_ablation.csv", index=False)
rows = []
adjusted = []
randomrows = []
overlaprows = []
diagnostics = []
for count, (cohort, s) in enumerate(scoresets.items()):
    valid = [c for c in s if s[c].notna().all() and s[c].nunique() > 1]
    a = s[valid].to_numpy()
    cm = corr_matrix(a)
    bm = boot_matrix(a, SEED + count)
    ix = {c: j for j, c in enumerate(valid)}
    n = len(s)
    for v in variants:
        for target in targets:
            if v not in ix or target not in ix:
                rows.append(
                    dict(
                        cohort=cohort,
                        program=v,
                        target=target,
                        n=n,
                        status="NOT_EVALUABLE",
                    )
                )
                continue
            i, j = ix[v], ix[target]
            rho = cm[i, j]
            br = bm[:, i, j]
            lo, hi = ci(br)
            p = stats.spearmanr(s[v], s[target]).pvalue
            delta = rho - cm[ix["A20"], j]
            bd = br - bm[:, ix["A20"], j]
            rows.append(
                dict(
                    cohort=cohort,
                    program=v,
                    target=target,
                    n=n,
                    spearman_rho=rho,
                    ci_low=lo,
                    ci_high=hi,
                    p_value=p,
                    delta_rho_vs20=delta,
                    delta_ci_low=ci(bd)[0],
                    delta_ci_high=ci(bd)[1],
                    delta_p_value=diff_p(bd, delta) if v != "A20" else np.nan,
                    status="EVALUABLE",
                    resampling="patient paired bootstrap B5000",
                )
            )
            nt = target + "__nonoverlap"
            k = ix[nt]
            r = cm[i, k]
            bb = bm[:, i, k]
            dd = r - rho
            db = bb - br
            overlaprows.append(
                dict(
                    cohort=cohort,
                    program=v,
                    target=target,
                    n=n,
                    spearman_rho_before=rho,
                    spearman_rho_after=r,
                    ci_low=ci(bb)[0],
                    ci_high=ci(bb)[1],
                    p_value=stats.spearmanr(s[v], s[nt]).pvalue,
                    delta_rho_after_minus_before=dd,
                    delta_ci_low=ci(db)[0],
                    delta_ci_high=ci(db)[1],
                    delta_p_value=(
                        diff_p(db, dd)
                        if k != j and not np.all(s[nt] == s[target])
                        else np.nan
                    ),
                )
            )
            for cov in [
                "glycolysis7_no_PKM",
                "global_expression",
                "detected_fraction",
                "E2F",
                "G2M",
                "purity_WGS",
                "stromal_score",
            ]:
                sub = s[[v, nt, cov]].dropna() if cov in s else pd.DataFrame()
                if (
                    len(sub) < 5
                    or cov not in valid
                    and (len(sub) == 0 or sub[cov].nunique() < 2)
                ):
                    adjusted.append(
                        dict(
                            cohort=cohort,
                            program=v,
                            target=target,
                            covariate=cov,
                            n=len(sub),
                            status="NOT_EVALUABLE",
                            reason="unavailable or invariant covariate",
                        )
                    )
                    continue
                if sub[cov].notna().all() and len(sub) == n and cov in ix:
                    cc = cm
                    bc = bm
                    ii, jj, kk = ix[v], ix[nt], ix[cov]
                else:
                    cc = corr_matrix(sub.to_numpy())
                    bc = boot_matrix(sub.to_numpy(), SEED + count)
                    ii, jj, kk = 0, 1, 2
                rr = partial(cc, ii, jj, kk)
                brp = partial(bc, ii, jj, kk)
                nn = len(sub)
                pval = 2 * stats.t.sf(
                    abs(rr) * np.sqrt((nn - 3) / max(1e-12, 1 - rr**2)), nn - 3
                )
                adjusted.append(
                    dict(
                        cohort=cohort,
                        program=v,
                        target=target,
                        covariate=cov,
                        n=nn,
                        partial_spearman=rr,
                        ci_low=ci(brp)[0],
                        ci_high=ci(brp)[1],
                        p_value=pval,
                        unadjusted_rho=cc[ii, jj],
                        program_covariate_rho=cc[ii, kk],
                        target_covariate_rho=cc[jj, kk],
                        VIF=1 / (1 - cc[ii, kk] ** 2),
                        status="EVALUABLE",
                        inference="single covariate partial rank correlation; approximate t P; patient bootstrap CI",
                    )
                )
    for cov in [
        "glycolysis7_no_PKM",
        "global_expression",
        "E2F",
        "G2M",
        "purity_WGS",
        "stromal_score",
        "CAF_xCell",
    ]:
        if cov in ix:
            i, j = ix["A20"], ix[cov]
            diagnostics.append(
                dict(
                    cohort=cohort,
                    covariate=cov,
                    n=n,
                    rho=cm[i, j],
                    ci_low=ci(bm[:, i, j])[0],
                    ci_high=ci(bm[:, i, j])[1],
                    p_value=stats.spearmanr(s["A20"], s[cov]).pvalue,
                )
            )
    x = data[cohort]
    z = x.sub(x.mean(axis=1), axis=0).div(
        x.std(axis=1, ddof=0).replace(0, np.nan), axis=0
    )
    for target in targets:
        y = s[target + "__nonoverlap"]
        actual = stats.spearmanr(s.A20, y).statistic
        effects = []
        for _, r in frozen[frozen.cohort == cohort].iterrows():
            rr = stats.spearmanr(
                score(z, r.positive.split("|")) - score(z, r.negative.split("|")), y
            ).statistic
            effects.append(rr)
            randomrows.append(
                dict(
                    cohort=cohort,
                    target=target,
                    replicate=r.replicate,
                    n=n,
                    random_rho=rr,
                    actual_rho=actual,
                    seed=SEED,
                )
            )
        effects = np.asarray(effects)
        qlo, qhi = np.quantile(effects, [0.025, 0.975])
        emp = (1 + sum(abs(effects) >= abs(actual))) / 201
        randomrows.append(
            dict(
                cohort=cohort,
                target=target,
                replicate="SUMMARY",
                n=n,
                actual_rho=actual,
                random_median=np.median(effects),
                null_effect_2_5pct=qlo,
                null_effect_97_5pct=qhi,
                empirical_two_sided_p=emp,
                percentile_signed=(sum(effects < actual) + 0.5 * sum(effects == actual))
                / 200,
                seed=SEED,
                interpretation="null distribution percentile; not causal; null interval not patient CI",
            )
        )
    print("DONE", cohort, flush=True)
res = fdr(pd.DataFrame(rows), ["cohort"])
res = fdr(res, ["cohort"], p="delta_p_value", q="delta_q_BH")
res.to_csv(T / "01_cross_cohort_specificity.csv", index=False)
ov = pd.DataFrame(overlaprows).merge(pd.DataFrame(overlap), on=["cohort", "target"])
fdr(ov, ["cohort"]).to_csv(T / "01_gene_overlap_control.csv", index=False)
fdr(pd.DataFrame(adjusted), ["cohort"]).to_csv(
    T / "01_covariate_specificity.csv", index=False
)
fdr(pd.DataFrame(diagnostics), ["cohort"]).to_csv(
    T / "01_covariate_associations.csv", index=False
)
random = pd.DataFrame(randomrows)
mask = random.replicate == "SUMMARY"
random.loc[mask, "empirical_q_BH"] = bh(random.loc[mask, "empirical_two_sided_p"])
random.to_csv(T / "01_matched_random_controls.csv", index=False)
(O / "logs/bulk_run_metadata.json").write_text(
    json.dumps(
        dict(
            seed=SEED,
            B=B,
            random_sets=200,
            random_design_sha256=hashlib.sha256(
                (T / "random_sets_frozen_before_associations.csv").read_bytes()
            ).hexdigest(),
            script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            limitations=[
                "GEO record identity not proven beyond GSM",
                "available20 scores differ in coverage",
                "single-covariate controls do not imply causal independence",
            ],
        ),
        indent=2,
    )
)
