#!/usr/bin/env python3
"""Specificity and sensitivity analyses of prespecified PKM-associated spatial regions."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import spearmanr, wilcoxon


from .settings import DATA, OUTPUT as OUTPUT_ROOT

ROOT = DATA / "spatial"
INPUT = ROOT / "GSE251950_all_spot_scores.csv.gz"
OUTPUT = OUTPUT_ROOT / "spatial_sensitivity"

OUTCOMES = {
    "Glycolysis_score_no_PKM": "Glycolysis excluding PKM",
    "Hypoxia_score": "Hypoxia",
    "Tcell_inflamed": "T-cell-inflamed",
    "IFNG_score": "IFNG",
    "Immune_fraction": "Tangram immune abundance",
}


def zscore(values: pd.Series) -> pd.Series:
    values = pd.to_numeric(values, errors="coerce")
    sd = values.std(ddof=0)
    if not np.isfinite(sd) or sd <= 1e-12:
        return pd.Series(0.0, index=values.index)
    return (values - values.mean()) / sd


def bh_adjust(values: pd.Series) -> np.ndarray:
    p = pd.to_numeric(values, errors="coerce").to_numpy(float)
    out = np.full(p.shape, np.nan)
    finite = np.isfinite(p)
    if not finite.any():
        return out
    x = p[finite]
    order = np.argsort(x)
    ranked = x[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.clip(adjusted, 0, 1)
    out[finite] = restored
    return out


def signed_rank(values: pd.Series) -> tuple[float, int]:
    x = pd.to_numeric(values, errors="coerce").to_numpy(float)
    x = x[np.isfinite(x)]
    x = x[np.abs(x) > 1e-12]
    if x.size < 3:
        return math.nan, int(x.size)
    try:
        return float(wilcoxon(x, alternative="two-sided", method="auto").pvalue), int(
            x.size
        )
    except ValueError:
        return math.nan, int(x.size)


def neighbor_lists(frame: pd.DataFrame) -> list[np.ndarray]:
    coordinates = frame[["pixel_x_fullres", "pixel_y_fullres"]].to_numpy(float)
    tree = cKDTree(coordinates)
    distances, indices = tree.query(coordinates, k=min(7, len(frame)))
    if indices.ndim == 1:
        indices = indices[:, None]
        distances = distances[:, None]
    spacing = float(np.nanmedian(distances[:, 1])) if distances.shape[1] > 1 else 1.0
    cutoff = spacing * 1.45
    return [
        indices[row][(distances[row] > 0) & (distances[row] <= cutoff)]
        for row in range(len(frame))
    ]


def prepare_sample(sample: pd.DataFrame, qc_mode: str) -> pd.DataFrame:
    sample = sample.copy().reset_index(drop=True)
    if qc_mode == "exclude_low_UMI_10pct":
        cutoff = sample["total_counts"].quantile(0.10)
        sample = sample.loc[sample["total_counts"] >= cutoff].reset_index(drop=True)
    elif qc_mode != "all_tissue_spots":
        raise ValueError(qc_mode)

    for column in [
        "PKM_expression",
        "RPL36AL_expression",
        "TCDP_driver_score",
        "PKM_RPL36AL_axis",
        "Fibroblast_fraction",
        "Fibroblast_marker_score",
        *OUTCOMES,
    ]:
        sample[f"sensitivity_z_{column}"] = zscore(sample[column])

    sample["composite_tangram"] = (
        sample["sensitivity_z_PKM_expression"]
        + sample["sensitivity_z_Fibroblast_fraction"]
    ) / 2
    sample["composite_marker"] = (
        sample["sensitivity_z_PKM_expression"]
        + sample["sensitivity_z_Fibroblast_marker_score"]
    ) / 2
    return sample


def define_groups(
    sample: pd.DataFrame, definition: str
) -> tuple[np.ndarray, np.ndarray]:
    if definition == "composite_tangram_all_spots":
        score = sample["composite_tangram"]
        eligible = np.ones(len(sample), dtype=bool)
    elif definition == "composite_marker_all_spots":
        score = sample["composite_marker"]
        eligible = np.ones(len(sample), dtype=bool)
    elif definition == "PKM_within_tangram_fibroblast_enriched":
        score = sample["PKM_expression"]
        eligible = sample["Fibroblast_fraction"].to_numpy() >= sample[
            "Fibroblast_fraction"
        ].quantile(2 / 3)
    elif definition == "PKM_within_marker_fibroblast_enriched":
        score = sample["PKM_expression"]
        eligible = sample["Fibroblast_marker_score"].to_numpy() >= sample[
            "Fibroblast_marker_score"
        ].quantile(2 / 3)
    elif definition == "TCDP_within_tangram_fibroblast_enriched":
        score = sample["TCDP_driver_score"]
        eligible = sample["Fibroblast_fraction"].to_numpy() >= sample[
            "Fibroblast_fraction"
        ].quantile(2 / 3)
    elif definition == "PKM_RPL36AL_within_tangram_fibroblast_enriched":
        score = sample["PKM_RPL36AL_axis"]
        eligible = sample["Fibroblast_fraction"].to_numpy() >= sample[
            "Fibroblast_fraction"
        ].quantile(2 / 3)
    else:
        raise ValueError(definition)

    cutoff = score.loc[eligible].quantile(0.75)
    high = eligible & (score.to_numpy() >= cutoff)
    reference = eligible & ~high
    return high, reference


def run_effects(data: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    definitions = [
        "composite_tangram_all_spots",
        "composite_marker_all_spots",
        "PKM_within_tangram_fibroblast_enriched",
        "PKM_within_marker_fibroblast_enriched",
        "TCDP_within_tangram_fibroblast_enriched",
        "PKM_RPL36AL_within_tangram_fibroblast_enriched",
    ]
    rows = []
    for qc_mode in ["all_tissue_spots", "exclude_low_UMI_10pct"]:
        for sample_name, raw_sample in data.groupby("sample", sort=True):
            sample = prepare_sample(raw_sample, qc_mode)
            neighbors = neighbor_lists(sample)
            for definition in definitions:
                high, reference = define_groups(sample, definition)
                for outcome, outcome_label in OUTCOMES.items():
                    values = sample[f"sensitivity_z_{outcome}"].to_numpy(float)
                    neighbor_values = np.array(
                        [np.nanmean(values[x]) if len(x) else np.nan for x in neighbors]
                    )
                    for analysis, evaluated in [
                        ("direct", values),
                        ("neighborhood", neighbor_values),
                    ]:
                        rows.append(
                            {
                                "sample": sample_name,
                                "qc_mode": qc_mode,
                                "definition": definition,
                                "analysis": analysis,
                                "outcome": outcome,
                                "outcome_label": outcome_label,
                                "n_spots_after_qc": len(sample),
                                "n_high": int(high.sum()),
                                "n_reference": int(reference.sum()),
                                "high_mean_z": float(np.nanmean(evaluated[high])),
                                "reference_mean_z": float(
                                    np.nanmean(evaluated[reference])
                                ),
                                "delta_mean_z": float(
                                    np.nanmean(evaluated[high])
                                    - np.nanmean(evaluated[reference])
                                ),
                            }
                        )
    effects = pd.DataFrame(rows)
    summary_rows = []
    grouping = ["qc_mode", "definition", "analysis", "outcome", "outcome_label"]
    for keys, group in effects.groupby(grouping, sort=False):
        p_value, n_nonzero = signed_rank(group["delta_mean_z"])
        summary_rows.append(
            dict(
                zip(grouping, keys),
                n_primary_sections=int(group["sample"].nunique()),
                n_nonzero_sections=n_nonzero,
                mean_delta_z=float(group["delta_mean_z"].mean()),
                median_delta_z=float(group["delta_mean_z"].median()),
                min_delta_z=float(group["delta_mean_z"].min()),
                max_delta_z=float(group["delta_mean_z"].max()),
                wilcoxon_p=p_value,
            )
        )
    summary = pd.DataFrame(summary_rows)
    summary["BH_FDR_within_definition"] = summary.groupby(
        ["qc_mode", "definition", "analysis"], sort=False
    )["wilcoxon_p"].transform(lambda x: bh_adjust(x))
    return effects, summary


def run_correlations(data: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    pairs = [
        ("PKM_expression", "Fibroblast_fraction", "PKM vs Tangram fibroblast"),
        ("PKM_expression", "Fibroblast_marker_score", "PKM vs fibroblast markers"),
        (
            "Fibroblast_fraction",
            "Fibroblast_marker_score",
            "Tangram vs fibroblast markers",
        ),
        ("Fibroblast_fraction", "Immune_fraction", "Tangram fibroblast vs immune"),
        ("Fibroblast_fraction", "total_counts", "Tangram fibroblast vs total UMI"),
        ("Fibroblast_marker_score", "total_counts", "Fibroblast markers vs total UMI"),
    ]
    rows = []
    for qc_mode in ["all_tissue_spots", "exclude_low_UMI_10pct"]:
        for sample_name, raw_sample in data.groupby("sample", sort=True):
            sample = prepare_sample(raw_sample, qc_mode)
            for x, y, label in pairs:
                rho, p_value = spearmanr(sample[x], sample[y], nan_policy="omit")
                rows.append(
                    {
                        "sample": sample_name,
                        "qc_mode": qc_mode,
                        "comparison": label,
                        "x": x,
                        "y": y,
                        "rho": float(rho),
                        "spot_level_p": float(p_value),
                    }
                )
    correlations = pd.DataFrame(rows)
    summaries = []
    for keys, group in correlations.groupby(
        ["qc_mode", "comparison", "x", "y"], sort=False
    ):
        p_value, n_nonzero = signed_rank(group["rho"])
        summaries.append(
            {
                "qc_mode": keys[0],
                "comparison": keys[1],
                "x": keys[2],
                "y": keys[3],
                "n_primary_sections": int(group["sample"].nunique()),
                "n_nonzero_sections": n_nonzero,
                "mean_rho": float(group["rho"].mean()),
                "median_rho": float(group["rho"].median()),
                "min_rho": float(group["rho"].min()),
                "max_rho": float(group["rho"].max()),
                "wilcoxon_p": p_value,
            }
        )
    summary = pd.DataFrame(summaries)
    summary["BH_FDR_within_qc_mode"] = summary.groupby("qc_mode", sort=False)[
        "wilcoxon_p"
    ].transform(lambda x: bh_adjust(x))
    return correlations, summary


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(INPUT)
    primary_flag = data["is_primary"]
    if primary_flag.dtype == object:
        primary_flag = primary_flag.astype(str).str.lower().isin(["true", "1"])
    else:
        primary_flag = primary_flag.astype(bool)
    primary = data.loc[primary_flag].copy()

    effects, effect_summary = run_effects(primary)
    correlations, correlation_summary = run_correlations(primary)
    effects.to_csv(OUTPUT / "PKM_spatial_sensitivity_sample_effects.csv", index=False)
    effect_summary.to_csv(OUTPUT / "PKM_spatial_sensitivity_summary.csv", index=False)
    correlations.to_csv(
        OUTPUT / "PKM_spatial_sensitivity_sample_correlations.csv", index=False
    )
    correlation_summary.to_csv(
        OUTPUT / "PKM_spatial_sensitivity_correlation_summary.csv", index=False
    )

    key = effect_summary.loc[
        (effect_summary["qc_mode"] == "all_tissue_spots")
        & (effect_summary["analysis"] == "neighborhood")
    ].copy()
    payload = {
        "input": str(INPUT),
        "n_primary_sections": int(primary["sample"].nunique()),
        "n_primary_spots": int(len(primary)),
        "definitions": sorted(effect_summary["definition"].unique().tolist()),
        "key_neighborhood_results": key.to_dict(orient="records"),
    }
    with (OUTPUT / "PKM_spatial_sensitivity_key_results.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(payload, handle, indent=2)
    print(effect_summary.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
