#!/usr/bin/env python3
"""CPU-only cross-cohort transcriptional associations of the fixed host transcriptional program."""

from __future__ import annotations

import csv
import gzip
import json
import math
import re
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats


from .settings import DATA, OUTPUT

BASE = DATA / "external"
OUT = OUTPUT / "cohorts"
TABLES = OUT / "tables"
SCORES = OUT / "per_sample_scores"
FIGURES = OUT / "figures"
for directory in (OUT, TABLES, SCORES, FIGURES):
    directory.mkdir(parents=True, exist_ok=True)

RANDOM_SEED = 2026

MIP1_GENES = [
    "AP2M1",
    "M6PR",
    "RNF167",
    "PKM",
    "BFAR",
    "TNFRSF12A",
    "CCT3",
    "TPM4",
    "RTFDC1",
    "ALAS1",
]
MIP2_GENES = [
    "PLAC8L1",
    "SMN1",
    "ZNF714",
    "PBDC1",
    "DDX46",
    "RPL36AL",
    "RP11-437B10.1",
    "NEURL2",
    "POTEI",
    "CTAGE8",
]

FIXED_SIGNATURES = {
    "CAF_composite": ["FAP", "COL1A1", "COL3A1", "PDGFRA"],
    "IFNG_signature": ["IFNG", "STAT1", "CXCL9", "CXCL10", "IDO1", "HLA-DRA"],
    "CYT_signature": ["GZMA", "PRF1"],
    "T_cell_inflamed": ["CD8A", "CXCL9", "CXCL10", "GZMB", "IFNG", "STAT1", "TBX21"],
}

INDIVIDUAL_TARGETS = ["FAP", "COL1A1", "COL3A1", "PDGFRA"]

COHORTS = {
    "GSE62254": {
        "matrix": BASE / "bulk/GSE62254/GSE62254_series_matrix.txt.gz",
        "annotation": BASE / "annotation/GPL570.annot.gz",
        "platform": "GPL570",
        "filter_excluded_titles": False,
    },
    "GSE15459": {
        "matrix": BASE / "bulk/GSE15459/GSE15459_series_matrix.txt.gz",
        "annotation": BASE / "annotation/GPL570.annot.gz",
        "platform": "GPL570",
        "filter_excluded_titles": True,
    },
    "GSE84433": {
        "matrix": BASE / "bulk/GSE84433/GSE84433_series_matrix.txt.gz",
        "annotation": BASE / "annotation/GPL6947.annot.gz",
        "platform": "GPL6947",
        "filter_excluded_titles": False,
    },
    "GSE26253": {
        "matrix": BASE / "bulk/GSE26253/GSE26253_series_matrix.txt.gz",
        "annotation": BASE
        / "annotation/GPL8432_HUMANREF-8_V3_0_R1_11282963_A_WGDASL.txt.gz",
        "platform": "GPL8432",
        "filter_excluded_titles": False,
    },
}


def clean_field(value: str) -> str:
    return value.strip().strip('"').strip()


def parse_symbol_field(value: str) -> list[str]:
    raw = clean_field(value)
    if not raw or raw in {"---", "NA", "N/A", "null"}:
        return []
    symbols = []
    for part in re.split(r"\s*///\s*|\s*//\s*|\s*;\s*", raw):
        symbol = part.strip().upper()
        if symbol and symbol not in {"---", "NA", "N/A", "NULL"}:
            symbols.append(symbol)
    return list(dict.fromkeys(symbols))


def read_gmt(path: Path) -> dict[str, list[str]]:
    result = {}
    with path.open("rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) >= 3:
                result[fields[0]] = [g.upper() for g in fields[2:] if g]
    return result


def load_probe_map(
    path: Path, target_genes: set[str]
) -> tuple[dict[str, str], dict[str, int]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = None
        for candidate in reader:
            cleaned = [clean_field(x) for x in candidate]
            is_geo = "ID" in cleaned and "Gene symbol" in cleaned
            is_illumina = "Probe_Id" in cleaned and "Symbol" in cleaned
            if is_geo or is_illumina:
                header = cleaned
                break
        if header is None:
            raise ValueError(f"Cannot locate annotation table header in {path}")
        positions = {name: idx for idx, name in enumerate(header)}
        id_idx = positions.get("ID", positions.get("Probe_Id"))
        symbol_idx = positions.get("Gene symbol", positions.get("Symbol"))
        synonym_idx = positions.get("Synonyms")
        if id_idx is None or symbol_idx is None:
            raise ValueError(f"Cannot find probe/symbol columns in {path}")
        mapping = {}
        counts = defaultdict(int)
        aliases = {"PKM2": "PKM"}
        required_width = max(
            id_idx, symbol_idx, synonym_idx if synonym_idx is not None else 0
        )
        for row in reader:
            if len(row) <= required_width:
                continue
            probe = clean_field(row[id_idx])
            symbols = parse_symbol_field(row[symbol_idx])
            if synonym_idx is not None:
                symbols.extend(parse_symbol_field(row[synonym_idx]))
            canonical = [aliases.get(g, g) for g in symbols]
            relevant = list(dict.fromkeys(g for g in canonical if g in target_genes))
            if len(relevant) == 1:
                mapping[probe] = relevant[0]
                counts[relevant[0]] += 1
        return mapping, dict(counts)


def parse_geo_metadata(path: Path) -> dict[str, list[str]]:
    metadata: dict[str, list[str]] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!Sample_"):
                parts = line.rstrip("\r\n").split("\t")
                metadata[parts[0]] = [clean_field(x) for x in parts[1:]]
    return metadata


def float_array(fields: list[str]) -> np.ndarray:
    values = np.empty(len(fields), dtype=float)
    for idx, value in enumerate(fields):
        value = clean_field(value)
        try:
            values[idx] = float(value)
        except (TypeError, ValueError):
            values[idx] = np.nan
    return values


def read_geo_expression(
    matrix_path: Path,
    probe_map: dict[str, str],
    filter_excluded_titles: bool,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, int]]:
    metadata = parse_geo_metadata(matrix_path)
    titles = metadata.get("!Sample_title", [])
    sample_geo = metadata.get("!Sample_geo_accession", [])
    collected: dict[str, list[np.ndarray]] = defaultdict(list)
    sample_ids: list[str] = []

    with gzip.open(matrix_path, "rt", encoding="utf-8", errors="replace") as handle:
        in_table = False
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if not in_table:
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            parts = line.rstrip("\r\n").split("\t")
            if not sample_ids:
                sample_ids = [clean_field(x) for x in parts[1:]]
                continue
            probe = clean_field(parts[0])
            gene = probe_map.get(probe)
            if gene is not None and len(parts) == len(sample_ids) + 1:
                collected[gene].append(float_array(parts[1:]))

    if sample_geo and sample_geo != sample_ids:
        raise ValueError(f"Sample order mismatch in {matrix_path.name}")
    if titles and len(titles) != len(sample_ids):
        raise ValueError(f"Title count mismatch in {matrix_path.name}")

    keep = np.ones(len(sample_ids), dtype=bool)
    exclusion_reason = np.array([""] * len(sample_ids), dtype=object)
    if filter_excluded_titles:
        excluded = np.array(
            ["[EXCLUDED]" in title.upper() for title in titles], dtype=bool
        )
        keep &= ~excluded
        exclusion_reason[excluded] = "GEO title marked [EXCLUDED]"

    sample_qc = pd.DataFrame(
        {
            "sample_id": sample_ids,
            "sample_title": titles if titles else sample_ids,
            "included": keep,
            "exclusion_reason": exclusion_reason,
        }
    )

    aggregated = {}
    probe_counts = {}
    for gene, arrays in collected.items():
        stacked = np.vstack(arrays)
        aggregated[gene] = np.nanmedian(stacked, axis=0)[keep]
        probe_counts[gene] = len(arrays)
    expression = pd.DataFrame(aggregated, index=np.array(sample_ids)[keep]).T
    expression.index.name = "gene"
    expression.columns.name = "sample_id"
    return expression, sample_qc, probe_counts


def zscore_rows(expression: pd.DataFrame) -> pd.DataFrame:
    means = expression.mean(axis=1, skipna=True)
    stds = expression.std(axis=1, skipna=True, ddof=0).replace(0, np.nan)
    return expression.sub(means, axis=0).div(stds, axis=0)


def score_signature(
    z_expression: pd.DataFrame, genes: list[str]
) -> tuple[pd.Series, list[str]]:
    available = [g for g in genes if g in z_expression.index]
    if not available:
        return pd.Series(np.nan, index=z_expression.columns, dtype=float), []
    return z_expression.loc[available].mean(axis=0, skipna=True), available


def fisher_ci(rho: float, n: int, alpha: float = 0.05) -> tuple[float, float]:
    if not np.isfinite(rho) or n <= 3:
        return np.nan, np.nan
    clipped = np.clip(rho, -0.999999, 0.999999)
    z = np.arctanh(clipped)
    delta = stats.norm.ppf(1 - alpha / 2) / math.sqrt(n - 3)
    return float(np.tanh(z - delta)), float(np.tanh(z + delta))


def correlation(x: pd.Series, y: pd.Series) -> dict[str, float]:
    pair = pd.concat([x, y], axis=1).dropna()
    n = len(pair)
    if n < 4 or pair.iloc[:, 0].nunique() < 2 or pair.iloc[:, 1].nunique() < 2:
        return {
            "n": n,
            "spearman_rho": np.nan,
            "p_value": np.nan,
            "ci_low": np.nan,
            "ci_high": np.nan,
        }
    rho, p_value = stats.spearmanr(pair.iloc[:, 0], pair.iloc[:, 1])
    ci_low, ci_high = fisher_ci(float(rho), n)
    return {
        "n": n,
        "spearman_rho": float(rho),
        "p_value": float(p_value),
        "ci_low": ci_low,
        "ci_high": ci_high,
    }


def bh_adjust(values: pd.Series) -> pd.Series:
    p = values.to_numpy(dtype=float)
    result = np.full(len(p), np.nan)
    valid = np.isfinite(p)
    valid_p = p[valid]
    if valid_p.size:
        order = np.argsort(valid_p)
        ranked = valid_p[order]
        adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
        adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
        adjusted = np.clip(adjusted, 0, 1)
        restored = np.empty_like(adjusted)
        restored[order] = adjusted
        result[np.where(valid)[0]] = restored
    return pd.Series(result, index=values.index)


def random_effects_meta(frame: pd.DataFrame) -> dict[str, float]:
    usable = frame.dropna(subset=["spearman_rho", "n"]).copy()
    usable = usable[usable["n"] > 3]
    k = len(usable)
    if not k:
        return {}
    yi = np.arctanh(
        np.clip(usable["spearman_rho"].to_numpy(float), -0.999999, 0.999999)
    )
    vi = 1.0 / (usable["n"].to_numpy(float) - 3.0)
    wi = 1.0 / vi
    fixed = np.sum(wi * yi) / np.sum(wi)
    q = float(np.sum(wi * (yi - fixed) ** 2))
    df = k - 1
    c_value = float(np.sum(wi) - np.sum(wi**2) / np.sum(wi))
    tau2 = max(0.0, (q - df) / c_value) if c_value > 0 else 0.0
    wr = 1.0 / (vi + tau2)
    pooled_z = float(np.sum(wr * yi) / np.sum(wr))
    se = float(math.sqrt(1.0 / np.sum(wr)))
    z_stat = pooled_z / se
    p_value = float(2 * stats.norm.sf(abs(z_stat)))
    low_z = pooled_z - 1.96 * se
    high_z = pooled_z + 1.96 * se
    i2 = max(0.0, (q - df) / q * 100.0) if q > 0 and df > 0 else 0.0
    p_heterogeneity = float(stats.chi2.sf(q, df)) if df > 0 else np.nan
    return {
        "k": k,
        "pooled_spearman_rho": float(np.tanh(pooled_z)),
        "ci_low": float(np.tanh(low_z)),
        "ci_high": float(np.tanh(high_z)),
        "p_value": p_value,
        "tau2_fisher_z": tau2,
        "Q": q,
        "Q_df": df,
        "p_heterogeneity": p_heterogeneity,
        "I2_percent": i2,
    }


def read_cct_targets(path: Path, targets: set[str]) -> pd.DataFrame:
    collected: dict[str, list[np.ndarray]] = defaultdict(list)
    with path.open("rt", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        sample_ids = [clean_field(x) for x in header[1:]]
        for line in handle:
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) != len(sample_ids) + 1:
                continue
            gene = clean_field(parts[0]).upper()
            if gene in targets:
                collected[gene].append(float_array(parts[1:]))
    aggregated = {}
    for gene, values in collected.items():
        stacked = np.vstack(values)
        if np.isfinite(stacked).any():
            aggregated[gene] = np.nanmedian(stacked, axis=0)
    frame = pd.DataFrame(aggregated, index=sample_ids).T
    frame.index.name = "gene"
    frame.columns.name = "sample_id"
    return frame


def hedges_g(high: np.ndarray, low: np.ndarray) -> float:
    high = np.asarray(high, dtype=float)
    low = np.asarray(low, dtype=float)
    high = high[np.isfinite(high)]
    low = low[np.isfinite(low)]
    n1, n0 = len(high), len(low)
    if n1 < 2 or n0 < 2:
        return np.nan
    pooled_var = ((n1 - 1) * np.var(high, ddof=1) + (n0 - 1) * np.var(low, ddof=1)) / (
        n1 + n0 - 2
    )
    if pooled_var <= 0:
        return np.nan
    d_value = (np.mean(high) - np.mean(low)) / math.sqrt(pooled_var)
    correction = 1.0 - 3.0 / (4.0 * (n1 + n0) - 9.0)
    return float(correction * d_value)


def bootstrap_hedges_ci(
    high: np.ndarray, low: np.ndarray, n_boot: int = 10000
) -> tuple[float, float]:
    high = np.asarray(high, dtype=float)
    low = np.asarray(low, dtype=float)
    high = high[np.isfinite(high)]
    low = low[np.isfinite(low)]
    rng = np.random.default_rng(RANDOM_SEED)
    estimates = np.empty(n_boot, dtype=float)
    for idx in range(n_boot):
        high_sample = rng.choice(high, size=len(high), replace=True)
        low_sample = rng.choice(low, size=len(low), replace=True)
        estimates[idx] = hedges_g(high_sample, low_sample)
    estimates = estimates[np.isfinite(estimates)]
    if not estimates.size:
        return np.nan, np.nan
    return tuple(np.quantile(estimates, [0.025, 0.975]).astype(float))


def format_p(value: float) -> str:
    if not np.isfinite(value):
        return "P = NA"
    if value < 0.001:
        return f"P = {value:.2e}"
    return f"P = {value:.3f}"


def save_figure(fig: plt.Figure, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def make_cross_cohort_figure(correlations: pd.DataFrame, meta: pd.DataFrame) -> None:
    display_targets = [
        "FAP",
        "COL1A1",
        "COL3A1",
        "PDGFRA",
        "CAF_composite",
        "Hallmark_Glycolysis",
        "Hallmark_Hypoxia",
        "IFNG_signature",
        "CYT_signature",
        "T_cell_inflamed",
    ]
    labels = {
        "CAF_composite": "CAF composite",
        "Hallmark_Glycolysis": "Hallmark glycolysis",
        "Hallmark_Hypoxia": "Hallmark hypoxia",
        "IFNG_signature": "IFNG signature",
        "CYT_signature": "CYT signature",
        "T_cell_inflamed": "T-cell inflamed",
    }
    cohorts = list(COHORTS)
    matrix = np.full((len(display_targets), len(cohorts)), np.nan)
    fdr = np.full_like(matrix, np.nan)
    for i, target in enumerate(display_targets):
        for j, cohort in enumerate(cohorts):
            row = correlations[
                (correlations["target"] == target) & (correlations["cohort"] == cohort)
            ]
            if len(row):
                matrix[i, j] = row.iloc[0]["spearman_rho"]
                fdr[i, j] = row.iloc[0]["fdr_global"]

    fig = plt.figure(figsize=(13.4, 9.2))
    grid = fig.add_gridspec(1, 2, width_ratios=[0.86, 1.55], wspace=0.34)
    ax_h = fig.add_subplot(grid[0, 0])
    cmap = LinearSegmentedColormap.from_list("assoc", ["#3B78B5", "#F6F6F3", "#D95F4B"])
    image = ax_h.imshow(matrix, cmap=cmap, vmin=-0.75, vmax=0.75, aspect="auto")
    ax_h.set_xticks(np.arange(len(cohorts)), labels=cohorts, rotation=35, ha="right")
    ax_h.set_yticks(
        np.arange(len(display_targets)),
        labels=[labels.get(x, x) for x in display_targets],
    )
    ax_h.set_title(
        "Cross-cohort association pattern", loc="left", fontweight="bold", pad=12
    )
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if np.isfinite(matrix[i, j]):
                star = "*" if np.isfinite(fdr[i, j]) and fdr[i, j] < 0.05 else ""
                color = "white" if abs(matrix[i, j]) > 0.48 else "#202020"
                ax_h.text(
                    j,
                    i,
                    f"{matrix[i, j]:.2f}{star}",
                    ha="center",
                    va="center",
                    color=color,
                    fontsize=7,
                )
    ax_h.set_xlabel("Independent gastric cancer cohort")
    ax_h.set_ylabel("Fixed molecular readout")
    ax_h.text(
        0.0,
        -0.16,
        "Values are Spearman rho; * global BH-FDR < 0.05",
        transform=ax_h.transAxes,
        fontsize=7,
        color="#555555",
    )
    cbar = fig.colorbar(image, ax=ax_h, fraction=0.045, pad=0.04)
    cbar.set_label("Spearman rho")
    ax_h.text(
        -0.20, 1.035, "a", transform=ax_h.transAxes, fontsize=13, fontweight="bold"
    )

    ax_f = fig.add_subplot(grid[0, 1])
    forest_targets = [
        "CAF_composite",
        "Hallmark_Glycolysis",
        "Hallmark_Hypoxia",
        "IFNG_signature",
        "CYT_signature",
        "T_cell_inflamed",
    ]
    y = 0
    y_ticks = []
    y_labels = []
    separators = []
    for target in forest_targets:
        ax_f.text(-0.79, y, labels[target], fontweight="bold", va="center", fontsize=8)
        y += 0.8
        subset = correlations[correlations["target"] == target].set_index("cohort")
        for cohort in cohorts:
            if cohort not in subset.index:
                continue
            row = subset.loc[cohort]
            ax_f.errorbar(
                row["spearman_rho"],
                y,
                xerr=[
                    [row["spearman_rho"] - row["ci_low"]],
                    [row["ci_high"] - row["spearman_rho"]],
                ],
                fmt="o",
                ms=3.4,
                color="#7C8794",
                ecolor="#AAB2BA",
                elinewidth=0.9,
                capsize=1.8,
            )
            y_ticks.append(y)
            y_labels.append(f"{cohort}  (n={int(row['n'])})")
            y += 0.7
        pooled = meta[meta["target"] == target]
        if len(pooled):
            row = pooled.iloc[0]
            ax_f.errorbar(
                row["pooled_spearman_rho"],
                y,
                xerr=[
                    [row["pooled_spearman_rho"] - row["ci_low"]],
                    [row["ci_high"] - row["pooled_spearman_rho"]],
                ],
                fmt="D",
                ms=4.8,
                color="#A62B68",
                ecolor="#A62B68",
                elinewidth=1.4,
                capsize=2.3,
            )
            y_ticks.append(y)
            y_labels.append(f"Random-effects pooled  (I2={row['I2_percent']:.0f}%)")
            ax_f.text(
                0.82,
                y,
                f"{row['pooled_spearman_rho']:.2f}",
                va="center",
                fontsize=7,
                color="#A62B68",
            )
            y += 0.8
        separators.append(y - 0.35)
        y += 0.35

    for sep in separators[:-1]:
        ax_f.axhline(sep, color="#E4E7EA", lw=0.7)
    ax_f.axvline(0, color="#7B838B", lw=0.8, ls="--")
    ax_f.set_xlim(-0.82, 0.92)
    ax_f.set_ylim(y - 0.2, -0.7)
    ax_f.set_yticks(y_ticks, labels=y_labels)
    ax_f.tick_params(axis="y", labelsize=6.4)
    ax_f.set_xlabel("Spearman rho (95% CI)")
    ax_f.set_title(
        "Cohort-specific effects and random-effects meta-analysis",
        loc="left",
        fontweight="bold",
        pad=12,
    )
    ax_f.text(
        -0.10, 1.035, "b", transform=ax_f.transAxes, fontsize=13, fontweight="bold"
    )
    ax_f.spines["left"].set_visible(False)
    ax_f.tick_params(axis="y", length=0)

    fig.suptitle(
        "Cross-cohort transcriptional associations of the fixed host transcriptional program",
        x=0.51,
        y=0.995,
        fontsize=12,
        fontweight="bold",
    )
    save_figure(fig, FIGURES / "Supplementary_Figure_cross_cohort_conservation")


def add_scatter(
    ax: plt.Axes, x: pd.Series, y: pd.Series, title: str, ylabel: str, color: str
) -> dict[str, float]:
    pair = pd.concat([x, y], axis=1).dropna()
    result = correlation(pair.iloc[:, 0], pair.iloc[:, 1])
    ax.scatter(
        pair.iloc[:, 0],
        pair.iloc[:, 1],
        s=16,
        color=color,
        alpha=0.68,
        edgecolor="white",
        linewidth=0.3,
    )
    if len(pair) >= 2:
        slope, intercept = np.polyfit(pair.iloc[:, 0], pair.iloc[:, 1], 1)
        xx = np.linspace(pair.iloc[:, 0].min(), pair.iloc[:, 0].max(), 100)
        ax.plot(xx, slope * xx + intercept, color="#25344F", lw=1.2)
    ax.set_title(title, loc="left", fontweight="bold")
    ax.set_xlabel("RNA-derived fixed host-program score")
    ax.set_ylabel(ylabel)
    ax.text(
        0.04,
        0.95,
        f"Spearman rho = {result['spearman_rho']:.2f}\n{format_p(result['p_value'])}\nn = {int(result['n'])}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=7,
    )
    return result


def make_cptac_figure(scores: pd.DataFrame, comparison: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 7.6))
    add_scatter(
        axes[0, 0],
        scores["RNA_driver_score"],
        scores["PKM_protein_z"],
        "Fixed host-program score vs PKM protein",
        "PKM protein abundance (z)",
        "#C95781",
    )

    ax = axes[0, 1]
    low = (
        scores.loc[scores["driver_group"] == "Low", "PKM_protein_z"].dropna().to_numpy()
    )
    high = (
        scores.loc[scores["driver_group"] == "High", "PKM_protein_z"]
        .dropna()
        .to_numpy()
    )
    box = ax.boxplot(
        [low, high],
        positions=[0, 1],
        widths=0.46,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "#1E2633", "linewidth": 1.2},
        whiskerprops={"color": "#4C5561"},
        capprops={"color": "#4C5561"},
    )
    box["boxes"][0].set(facecolor="#CFE2F3", edgecolor="#3478B9", linewidth=1.0)
    box["boxes"][1].set(facecolor="#F3CDDC", edgecolor="#B52C69", linewidth=1.0)
    rng = np.random.default_rng(RANDOM_SEED)
    ax.scatter(
        rng.normal(0, 0.055, len(low)),
        low,
        s=12,
        color="#3478B9",
        alpha=0.55,
        edgecolor="none",
    )
    ax.scatter(
        rng.normal(1, 0.055, len(high)),
        high,
        s=12,
        color="#B52C69",
        alpha=0.55,
        edgecolor="none",
    )
    row = comparison[comparison["protein_readout"] == "PKM_protein_z"].iloc[0]
    ax.set_xticks([0, 1], [f"Low\n(n={len(low)})", f"High\n(n={len(high)})"])
    ax.set_ylabel("PKM protein abundance (z)")
    ax.set_xlabel("Fixed host-program score group (median split)")
    ax.set_title(
        "PKM protein across fixed host-program score groups",
        loc="left",
        fontweight="bold",
    )
    ax.text(
        0.04,
        0.95,
        f"{format_p(row['mann_whitney_p'])}\nHedges' g = {row['hedges_g']:.2f}\n95% CI {row['g_ci_low']:.2f} to {row['g_ci_high']:.2f}",
        transform=ax.transAxes,
        va="top",
        fontsize=7,
    )

    add_scatter(
        axes[1, 0],
        scores["RNA_driver_score"],
        scores["Protein_Hallmark_Glycolysis"],
        "Fixed host-program score vs protein glycolysis",
        "Protein glycolysis score",
        "#D59638",
    )
    add_scatter(
        axes[1, 1],
        scores["RNA_driver_score"],
        scores["Protein_Hallmark_Hypoxia"],
        "Fixed host-program score vs protein hypoxia",
        "Protein hypoxia score",
        "#4D9B84",
    )
    for label, ax in zip("abcd", axes.flat):
        ax.text(
            -0.16, 1.05, label, transform=ax.transAxes, fontsize=12, fontweight="bold"
        )
    fig.suptitle(
        "RNA–protein associations of the fixed host transcriptional program in CPTAC-STAD",
        y=1.01,
        fontsize=12,
        fontweight="bold",
    )
    fig.tight_layout(w_pad=2.6, h_pad=2.4)
    save_figure(fig, FIGURES / "Supplementary_Figure_CPTAC_proteomic_validation")


def main() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": 8,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "axes.linewidth": 0.8,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )

    gmt = read_gmt(BASE / "signature/h.all.v2026.1.Hs.symbols.gmt")
    glycolysis = gmt["HALLMARK_GLYCOLYSIS"]
    hypoxia = gmt["HALLMARK_HYPOXIA"]
    signatures = dict(FIXED_SIGNATURES)
    signatures["Hallmark_Glycolysis"] = glycolysis
    signatures["Hallmark_Hypoxia"] = hypoxia

    target_genes = set(MIP1_GENES + MIP2_GENES + INDIVIDUAL_TARGETS)
    for genes in signatures.values():
        target_genes.update(genes)

    annotation_cache = {}
    correlation_rows = []
    coverage_rows = []
    cohort_rows = []

    for cohort, config in COHORTS.items():
        annotation = config["annotation"]
        if annotation not in annotation_cache:
            annotation_cache[annotation] = load_probe_map(annotation, target_genes)
        probe_map, annotation_counts = annotation_cache[annotation]
        expression, sample_qc, observed_probe_counts = read_geo_expression(
            config["matrix"], probe_map, config["filter_excluded_titles"]
        )
        sample_qc.to_csv(TABLES / f"{cohort}_sample_inclusion.csv", index=False)
        z_expression = zscore_rows(expression)

        mip1_score, mip1_available = score_signature(z_expression, MIP1_GENES)
        mip2_score, mip2_available = score_signature(z_expression, MIP2_GENES)
        driver_score = mip1_score - mip2_score
        sample_scores = pd.DataFrame(index=z_expression.columns)
        sample_scores.index.name = "sample_id"
        sample_scores["driver_score"] = driver_score
        sample_scores["MIP1_module_score"] = mip1_score
        sample_scores["MIP2_module_score"] = mip2_score

        for gene in INDIVIDUAL_TARGETS:
            if gene in z_expression.index:
                sample_scores[gene] = z_expression.loc[gene]
            else:
                sample_scores[gene] = np.nan
        for signature, genes in signatures.items():
            score, available = score_signature(z_expression, genes)
            sample_scores[signature] = score
            coverage_rows.append(
                {
                    "dataset": cohort,
                    "data_layer": "transcriptome",
                    "feature": signature,
                    "requested_genes": len(genes),
                    "available_genes": len(available),
                    "coverage_percent": 100.0 * len(available) / len(genes),
                    "genes_available": "|".join(available),
                    "genes_missing": "|".join([g for g in genes if g not in available]),
                }
            )
        coverage_rows.extend(
            [
                {
                    "dataset": cohort,
                    "data_layer": "transcriptome",
                    "feature": "MIP1_driver_module",
                    "requested_genes": len(MIP1_GENES),
                    "available_genes": len(mip1_available),
                    "coverage_percent": 100.0 * len(mip1_available) / len(MIP1_GENES),
                    "genes_available": "|".join(mip1_available),
                    "genes_missing": "|".join(
                        [g for g in MIP1_GENES if g not in mip1_available]
                    ),
                },
                {
                    "dataset": cohort,
                    "data_layer": "transcriptome",
                    "feature": "MIP2_driver_module",
                    "requested_genes": len(MIP2_GENES),
                    "available_genes": len(mip2_available),
                    "coverage_percent": 100.0 * len(mip2_available) / len(MIP2_GENES),
                    "genes_available": "|".join(mip2_available),
                    "genes_missing": "|".join(
                        [g for g in MIP2_GENES if g not in mip2_available]
                    ),
                },
            ]
        )

        sample_scores.reset_index().to_csv(
            SCORES / f"{cohort}_fixed_scores.csv", index=False
        )
        for target in INDIVIDUAL_TARGETS + list(signatures):
            result = correlation(sample_scores["driver_score"], sample_scores[target])
            correlation_rows.append({"cohort": cohort, "target": target, **result})
        cohort_rows.append(
            {
                "cohort": cohort,
                "platform": config["platform"],
                "matrix_samples": len(sample_qc),
                "included_samples": int(sample_qc["included"].sum()),
                "excluded_samples": int((~sample_qc["included"]).sum()),
                "observed_target_genes": len(expression),
                "MIP1_genes_available": len(mip1_available),
                "MIP2_genes_available": len(mip2_available),
                "probe_aggregation": "Median across probes mapping to one prespecified target gene",
                "normalization_for_scoring": "Within-cohort gene-wise z score",
            }
        )

    correlations = pd.DataFrame(correlation_rows)
    correlations["fdr_global"] = bh_adjust(correlations["p_value"])
    correlations["fdr_within_cohort"] = correlations.groupby(
        "cohort", group_keys=False
    )["p_value"].apply(bh_adjust)
    correlations.to_csv(TABLES / "cross_cohort_spearman_correlations.csv", index=False)

    meta_rows = []
    for target, group in correlations.groupby("target", sort=False):
        meta_rows.append({"target": target, **random_effects_meta(group)})
    meta = pd.DataFrame(meta_rows)
    meta["fdr_across_meta_tests"] = bh_adjust(meta["p_value"])
    meta.to_csv(TABLES / "cross_cohort_random_effects_meta_analysis.csv", index=False)
    pd.DataFrame(coverage_rows).to_csv(
        TABLES / "gene_signature_coverage.csv", index=False
    )
    pd.DataFrame(cohort_rows).to_csv(TABLES / "cohort_information.csv", index=False)

    # CPTAC matched RNA-protein associations.
    cptac_targets_rna = set(MIP1_GENES + MIP2_GENES)
    cptac_targets_protein = set(glycolysis + hypoxia + ["PKM"])
    rna = read_cct_targets(
        BASE / "CPTAC_STAD/GC_RNAseq_gene_rsem_rm_circRNA_UQ_log2(x+1)_Tumor.cct",
        cptac_targets_rna,
    )
    protein = read_cct_targets(
        BASE / "CPTAC_STAD/GC_TMT_proteomics_ratio_gene_MD_GencodeV42_Tumor.cct",
        cptac_targets_protein,
    )
    matched_samples = [sample for sample in rna.columns if sample in protein.columns]
    rna = rna.loc[:, matched_samples]
    protein = protein.loc[:, matched_samples]
    rna_z = zscore_rows(rna)
    protein_z = zscore_rows(protein)
    rna_mip1, rna_mip1_available = score_signature(rna_z, MIP1_GENES)
    rna_mip2, rna_mip2_available = score_signature(rna_z, MIP2_GENES)
    rna_driver = rna_mip1 - rna_mip2
    pkm_protein = protein_z.loc["PKM"]
    protein_glycolysis, protein_glycolysis_available = score_signature(
        protein_z, glycolysis
    )
    protein_hypoxia, protein_hypoxia_available = score_signature(protein_z, hypoxia)
    cptac_scores = pd.DataFrame(
        {
            "sample_id": matched_samples,
            "RNA_driver_score": rna_driver.reindex(matched_samples).to_numpy(),
            "PKM_protein_z": pkm_protein.reindex(matched_samples).to_numpy(),
            "Protein_Hallmark_Glycolysis": protein_glycolysis.reindex(
                matched_samples
            ).to_numpy(),
            "Protein_Hallmark_Hypoxia": protein_hypoxia.reindex(
                matched_samples
            ).to_numpy(),
        }
    ).set_index("sample_id")
    driver_median = float(cptac_scores["RNA_driver_score"].median())
    cptac_scores["driver_group"] = np.where(
        cptac_scores["RNA_driver_score"] > driver_median, "High", "Low"
    )
    cptac_scores.reset_index().to_csv(
        SCORES / "CPTAC_STAD_matched_RNA_protein_scores.csv", index=False
    )

    cptac_corr_rows = []
    for readout in [
        "PKM_protein_z",
        "Protein_Hallmark_Glycolysis",
        "Protein_Hallmark_Hypoxia",
    ]:
        spearman = correlation(cptac_scores["RNA_driver_score"], cptac_scores[readout])
        pair = cptac_scores[["RNA_driver_score", readout]].dropna()
        pearson_r, pearson_p = stats.pearsonr(pair.iloc[:, 0], pair.iloc[:, 1])
        cptac_corr_rows.append(
            {
                "protein_readout": readout,
                **spearman,
                "pearson_r": float(pearson_r),
                "pearson_p": float(pearson_p),
            }
        )
    cptac_corr = pd.DataFrame(cptac_corr_rows)
    cptac_corr["spearman_fdr"] = bh_adjust(cptac_corr["p_value"])
    cptac_corr.to_csv(
        TABLES / "CPTAC_RNA_driver_vs_protein_correlations.csv", index=False
    )

    group_rows = []
    for readout in [
        "PKM_protein_z",
        "Protein_Hallmark_Glycolysis",
        "Protein_Hallmark_Hypoxia",
    ]:
        high = (
            cptac_scores.loc[cptac_scores["driver_group"] == "High", readout]
            .dropna()
            .to_numpy()
        )
        low = (
            cptac_scores.loc[cptac_scores["driver_group"] == "Low", readout]
            .dropna()
            .to_numpy()
        )
        u_stat, p_value = stats.mannwhitneyu(high, low, alternative="two-sided")
        g_value = hedges_g(high, low)
        g_low, g_high = bootstrap_hedges_ci(high, low)
        group_rows.append(
            {
                "protein_readout": readout,
                "n_high": len(high),
                "n_low": len(low),
                "median_high": float(np.median(high)),
                "median_low": float(np.median(low)),
                "mann_whitney_u": float(u_stat),
                "mann_whitney_p": float(p_value),
                "hedges_g": g_value,
                "g_ci_low": g_low,
                "g_ci_high": g_high,
                "bootstrap_resamples": 10000,
            }
        )
    group_comparison = pd.DataFrame(group_rows)
    group_comparison["mann_whitney_fdr"] = bh_adjust(group_comparison["mann_whitney_p"])
    group_comparison.to_csv(
        TABLES / "CPTAC_driver_high_low_protein_comparisons.csv", index=False
    )

    cptac_coverage = [
        {
            "dataset": "CPTAC-STAD",
            "data_layer": "RNA",
            "feature": "MIP1_driver_module",
            "requested_genes": len(MIP1_GENES),
            "available_genes": len(rna_mip1_available),
            "coverage_percent": 100 * len(rna_mip1_available) / len(MIP1_GENES),
            "genes_available": "|".join(rna_mip1_available),
            "genes_missing": "|".join(
                [g for g in MIP1_GENES if g not in rna_mip1_available]
            ),
        },
        {
            "dataset": "CPTAC-STAD",
            "data_layer": "RNA",
            "feature": "MIP2_driver_module",
            "requested_genes": len(MIP2_GENES),
            "available_genes": len(rna_mip2_available),
            "coverage_percent": 100 * len(rna_mip2_available) / len(MIP2_GENES),
            "genes_available": "|".join(rna_mip2_available),
            "genes_missing": "|".join(
                [g for g in MIP2_GENES if g not in rna_mip2_available]
            ),
        },
        {
            "dataset": "CPTAC-STAD",
            "data_layer": "protein",
            "feature": "Hallmark_Glycolysis",
            "requested_genes": len(glycolysis),
            "available_genes": len(protein_glycolysis_available),
            "coverage_percent": 100
            * len(protein_glycolysis_available)
            / len(glycolysis),
            "genes_available": "|".join(protein_glycolysis_available),
            "genes_missing": "|".join(
                [g for g in glycolysis if g not in protein_glycolysis_available]
            ),
        },
        {
            "dataset": "CPTAC-STAD",
            "data_layer": "protein",
            "feature": "Hallmark_Hypoxia",
            "requested_genes": len(hypoxia),
            "available_genes": len(protein_hypoxia_available),
            "coverage_percent": 100 * len(protein_hypoxia_available) / len(hypoxia),
            "genes_available": "|".join(protein_hypoxia_available),
            "genes_missing": "|".join(
                [g for g in hypoxia if g not in protein_hypoxia_available]
            ),
        },
    ]
    pd.DataFrame(cptac_coverage).to_csv(TABLES / "CPTAC_gene_coverage.csv", index=False)

    make_cross_cohort_figure(correlations, meta)
    make_cptac_figure(cptac_scores, group_comparison)

    methods = {
        "analysis_name": "MultiSTAD cross-cohort transcriptional and RNA–protein associations",
        "random_seed": RANDOM_SEED,
        "driver_score": "mean(z[MIP1-high available genes]) - mean(z[MIP2-high available genes])",
        "gene_standardization": "Gene-wise z scoring within each cohort; cohorts were never merged",
        "probe_mapping": "Primary symbols and platform aliases were matched to fixed targets; probes mapping to multiple fixed target genes were excluded",
        "probe_aggregation": "Median expression across probes mapping uniquely to the same gene",
        "external_feature_selection": "None",
        "external_direction_reestimation": "None",
        "GSE15459_exclusion": "Samples whose GEO title contains [EXCLUDED] were removed before scoring",
        "association": "Two-sided Spearman correlation",
        "multiple_testing": "Benjamini-Hochberg FDR",
        "meta_analysis": "DerSimonian-Laird random-effects meta-analysis on Fisher-z transformed Spearman correlations",
        "CPTAC_grouping": "Median split of the RNA-derived fixed host-program score among matched RNA-protein tumor samples",
        "CPTAC_group_test": "Two-sided Mann-Whitney U test; Hedges' g with 10,000 bootstrap resamples",
        "MIP1_high_genes": MIP1_GENES,
        "MIP2_high_genes": MIP2_GENES,
        "fixed_signatures": signatures,
    }
    (OUT / "analysis_manifest.json").write_text(
        json.dumps(methods, indent=2), encoding="utf-8"
    )

    summary_lines = [
        "# MultiSTAD cross-cohort association summary",
        "",
        "## Cohorts",
        pd.DataFrame(cohort_rows).to_csv(index=False).strip(),
        "",
        "## Random-effects meta-analysis",
        meta[
            [
                "target",
                "pooled_spearman_rho",
                "ci_low",
                "ci_high",
                "p_value",
                "fdr_across_meta_tests",
                "I2_percent",
            ]
        ]
        .to_csv(index=False)
        .strip(),
        "",
        "## CPTAC RNA-protein correlations",
        cptac_corr.to_csv(index=False).strip(),
        "",
        "## CPTAC high versus low fixed host-program score comparisons",
        group_comparison.to_csv(index=False).strip(),
        "",
        "All analyses used the prespecified gene directions and no external-cohort feature selection.",
    ]
    (OUT / "RESULTS_SUMMARY.md").write_text("\n".join(summary_lines), encoding="utf-8")

    print("Analysis complete")
    print(pd.DataFrame(cohort_rows).to_string(index=False))
    print("\nMeta-analysis:")
    print(
        meta[
            [
                "target",
                "pooled_spearman_rho",
                "ci_low",
                "ci_high",
                "p_value",
                "I2_percent",
            ]
        ].to_string(index=False)
    )
    print("\nCPTAC correlations:")
    print(cptac_corr.to_string(index=False))
    print("\nCPTAC group comparisons:")
    print(group_comparison.to_string(index=False))


if __name__ == "__main__":
    main()
