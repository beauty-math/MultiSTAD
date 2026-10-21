#!/usr/bin/env python3
"""Spatial tissue correlates of the fixed host transcriptional program.

The script uses GSE163558 primary gastric-tumor single-cell data as a Tangram
reference and GSE251950 Visium sections as the spatial target. MultiSTAD gene
programs are fixed before spatial analysis and are never reselected here.
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import os
from pathlib import Path
import shutil
import tarfile
import tempfile
from typing import Iterable

import numpy as np

# Tangram 1.0.4 still references aliases removed from recent NumPy versions.
np.float = float  # type: ignore[attr-defined]
np.int = int  # type: ignore[attr-defined]
np.bool = bool  # type: ignore[attr-defined]

import anndata as ad
import pandas as pd
import scanpy as sc
from scipy.io import mmread
from scipy.spatial import cKDTree
from scipy.stats import spearmanr, wilcoxon
import tangram as tg


from .settings import DATA, OUTPUT

PROJECT_ROOT = DATA
SPATIAL_ROOT = DATA / "spatial"
ARCHIVE = SPATIAL_ROOT / "raw" / "GSE251950_RAW.tar"
WORK_ROOT = OUTPUT / "spatial_mapping/work"
RESULT_ROOT = OUTPUT / "spatial_mapping"
IMAGE_ROOT = RESULT_ROOT / "images"
SAMPLE_ROOT = RESULT_ROOT / "samples"
REFERENCE_ROOT = RESULT_ROOT / "reference"

SCRNA_ROOT = DATA / "single_cell/raw/GSE163558"
SCRNA_META = DATA / "single_cell/GSE163558_coarse_celltype_metadata.csv"

RANDOM_SEED = 2026
TARGET_SUM = 10_000.0
PRIMARY_REFERENCE_SAMPLES = [
    "GSM5004180_PT1",
    "GSM5004181_PT2",
    "GSM5004182_PT3",
]

SAMPLES = {
    "GC1": {
        "accession": "GSM7990473",
        "member": "GSM7990473_20_00331_LI_SING.tar.gz",
        "prefix": "20_00331_LI_SING",
        "primary": True,
    },
    "GC2": {
        "accession": "GSM7990474",
        "member": "GSM7990474_21_00731_LI_SING.tar.gz",
        "prefix": "21_00731_LI_SING",
        "primary": True,
    },
    "GC3": {
        "accession": "GSM7990475",
        "member": "GSM7990475_21_00732_LI_SING.tar.gz",
        "prefix": "21_00732_LI_SING",
        "primary": True,
    },
    "GC4": {
        "accession": "GSM7990476",
        "member": "GSM7990476_21_00733_LI_SING.tar.gz",
        "prefix": "21_00733_LI_SING",
        "primary": True,
    },
    "GC5": {
        "accession": "GSM7990477",
        "member": "GSM7990477_21_00734_LI_SING.tar.gz",
        "prefix": "21_00734_LI_SING",
        "primary": True,
    },
    "GC6": {
        "accession": "GSM7990478",
        "member": "GSM7990478_21_01251_LI_SING.tar.gz",
        "prefix": "21_01251_LI_SING",
        "primary": True,
    },
    "GC6-PM": {
        "accession": "GSM7990479",
        "member": "GSM7990479_21_01252_LI_SING.tar.gz",
        "prefix": "21_01252_LI_SING",
        "primary": False,
    },
    "GC7": {
        "accession": "GSM7990480",
        "member": "GSM7990480_21_01253_LI_SING.tar.gz",
        "prefix": "21_01253_LI_SING",
        "primary": True,
    },
    "GC8": {
        "accession": "GSM7990481",
        "member": "GSM7990481_21_01254_LI_SING.tar.gz",
        "prefix": "21_01254_LI_SING",
        "primary": True,
    },
    "GC9": {
        "accession": "GSM7990482",
        "member": "GSM7990482_21_01675_LI_SING.tar.gz",
        "prefix": "21_01675_LI_SING",
        "primary": True,
    },
}

CELLTYPE_MARKERS = {
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "MSLN", "MUC1"],
    "T_NK": ["CD3D", "CD3E", "TRAC", "NKG7", "GNLY", "KLRD1"],
    "B_Plasma": ["MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN", "CD74"],
    "Myeloid": ["LYZ", "LST1", "FCGR3A", "TYROBP", "S100A8", "S100A9"],
    "Fibroblast": ["COL1A1", "COL1A2", "DCN", "COL3A1", "PDGFRA", "TAGLN"],
    "Endothelial": ["PECAM1", "VWF", "KDR", "EMCN", "ESAM", "RAMP2"],
    "Mast": ["TPSAB1", "TPSB2", "KIT", "CPA3"],
    "Proliferating": ["MKI67", "TOP2A", "UBE2C", "STMN1"],
}

SIGNATURES = {
    "IFNG_score": ["IFNG", "STAT1", "CXCL9", "CXCL10", "IDO1", "HLA-DRA"],
    "Tcell_inflamed": ["CD8A", "CXCL9", "CXCL10", "GZMB", "IFNG", "STAT1", "TBX21"],
    "Glycolysis_score": [
        "HK2",
        "PFKP",
        "ALDOA",
        "GAPDH",
        "ENO1",
        "PKM",
        "LDHA",
        "SLC2A1",
    ],
    "Hypoxia_score": ["HIF1A", "VEGFA", "CA9", "SLC2A1", "ENO1", "LDHA"],
}

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

MAPPING_CELLTYPES = [
    "Epithelial",
    "T_NK",
    "B_Plasma",
    "Myeloid",
    "Fibroblast",
    "Endothelial",
    "Mast",
    "Proliferating",
]


def unique_in_order(values: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    result: list[str] = []
    for value in values:
        value = str(value).upper()
        if value not in seen:
            seen.add(value)
            result.append(value)
    return result


MAPPING_GENES = unique_in_order(
    gene for genes in CELLTYPE_MARKERS.values() for gene in genes
)
ANALYSIS_GENES = unique_in_order(
    MAPPING_GENES
    + [gene for genes in SIGNATURES.values() for gene in genes]
    + MIP1_GENES
    + MIP2_GENES
    + ["PKM", "RPL36AL"]
)


def read_features(path: Path) -> list[str]:
    frame = pd.read_csv(path, sep="\t", header=None, compression="gzip")
    column = 1 if frame.shape[1] > 1 else 0
    return frame.iloc[:, column].astype(str).str.upper().tolist()


def read_barcodes(path: Path) -> list[str]:
    frame = pd.read_csv(path, sep="\t", header=None, compression="gzip")
    return frame.iloc[:, 0].astype(str).tolist()


def read_matrix(path: Path):
    with gzip.open(path, "rb") as handle:
        matrix = mmread(handle).tocsr()
    return matrix


def normalized_selected_expression(
    matrix,
    feature_names: list[str],
    selected_genes: list[str],
) -> tuple[np.ndarray, list[str], np.ndarray]:
    if matrix.shape[0] != len(feature_names):
        raise ValueError(
            f"Feature count mismatch: matrix has {matrix.shape[0]} rows, "
            f"features has {len(feature_names)}"
        )
    totals = np.asarray(matrix.sum(axis=0)).ravel().astype(np.float64)
    totals[totals <= 0] = 1.0
    gene_to_rows: dict[str, list[int]] = {}
    for row_index, gene in enumerate(feature_names):
        if gene in selected_genes:
            gene_to_rows.setdefault(gene, []).append(row_index)

    present = [gene for gene in selected_genes if gene in gene_to_rows]
    expression = np.zeros((matrix.shape[1], len(present)), dtype=np.float32)
    for column_index, gene in enumerate(present):
        rows = gene_to_rows[gene]
        counts = np.asarray(matrix[rows, :].sum(axis=0)).ravel()
        expression[:, column_index] = np.log1p(counts / totals * TARGET_SUM)
    return expression, present, totals


def build_single_cell_reference(force: bool = False) -> ad.AnnData:
    REFERENCE_ROOT.mkdir(parents=True, exist_ok=True)
    reference_path = REFERENCE_ROOT / "GSE163558_primary_reference_44markers.h5ad"
    if reference_path.exists() and not force:
        print(f"Loading cached single-cell reference: {reference_path}", flush=True)
        return sc.read_h5ad(reference_path)

    metadata = pd.read_csv(SCRNA_META, index_col=0)
    metadata["sample_name"] = metadata["sample_name"].astype(str)
    metadata["cell_id_raw"] = metadata["cell_id_raw"].astype(str)
    metadata["predicted_celltype"] = metadata["predicted_celltype"].astype(str)
    metadata = metadata[
        metadata["sample_name"].isin(PRIMARY_REFERENCE_SAMPLES)
        & metadata["predicted_celltype"].isin(MAPPING_CELLTYPES)
    ].copy()

    expression_blocks: list[np.ndarray] = []
    observation_blocks: list[pd.DataFrame] = []
    genes_used: list[str] | None = None

    for sample in PRIMARY_REFERENCE_SAMPLES:
        print(f"Building scRNA reference block: {sample}", flush=True)
        matrix_path = SCRNA_ROOT / f"{sample}_matrix.mtx.gz"
        feature_path = SCRNA_ROOT / f"{sample}_features.tsv.gz"
        barcode_path = SCRNA_ROOT / f"{sample}_barcodes.tsv.gz"
        matrix = read_matrix(matrix_path)
        feature_names = read_features(feature_path)
        barcodes = read_barcodes(barcode_path)
        expression, present, _ = normalized_selected_expression(
            matrix, feature_names, MAPPING_GENES
        )
        if genes_used is None:
            genes_used = present
        elif present != genes_used:
            raise ValueError(f"Marker-gene order differs in {sample}")

        sample_meta = metadata[metadata["sample_name"] == sample].copy()
        sample_meta = sample_meta.set_index("cell_id_raw", drop=False)
        keep_positions = [
            i for i, barcode in enumerate(barcodes) if barcode in sample_meta.index
        ]
        keep_barcodes = [barcodes[i] for i in keep_positions]
        block_meta = sample_meta.loc[keep_barcodes].copy()
        block_meta.index = [f"{sample}__{barcode}" for barcode in keep_barcodes]
        expression_blocks.append(expression[keep_positions, :])
        observation_blocks.append(block_meta)
        del matrix, expression

    if genes_used is None:
        raise RuntimeError("No single-cell reference genes were loaded")

    combined_expression = np.vstack(expression_blocks).astype(np.float32)
    combined_observations = pd.concat(observation_blocks, axis=0)
    reference = ad.AnnData(
        X=combined_expression,
        obs=combined_observations,
        var=pd.DataFrame(index=genes_used),
    )
    reference.obs["coarse_celltype"] = pd.Categorical(
        reference.obs["predicted_celltype"], categories=MAPPING_CELLTYPES
    )
    reference.write_h5ad(reference_path, compression="gzip")

    counts = (
        reference.obs.groupby(["sample_name", "coarse_celltype"], observed=False)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    counts.to_csv(
        REFERENCE_ROOT / "GSE163558_primary_reference_cell_counts.csv", index=False
    )
    profile_rows = []
    for celltype in MAPPING_CELLTYPES:
        mask = np.asarray(reference.obs["coarse_celltype"] == celltype)
        if not mask.any():
            continue
        means = np.asarray(reference.X[mask].mean(axis=0)).ravel()
        for gene, value in zip(reference.var_names, means):
            profile_rows.append(
                {
                    "coarse_celltype": celltype,
                    "gene": gene,
                    "mean_log_expression": value,
                }
            )
    pd.DataFrame(profile_rows).to_csv(
        REFERENCE_ROOT / "GSE163558_primary_reference_profiles.csv", index=False
    )
    print(f"Reference cells: {reference.n_obs}; genes: {reference.n_vars}", flush=True)
    return reference


def extract_required_sample_files(
    outer: tarfile.TarFile, member_name: str, target: Path
) -> None:
    required_suffixes = (
        "_matrix.mtx.gz",
        "_barcodes.tsv.gz",
        "_features.tsv.gz",
        "_tissue_positions_list.csv",
        "_scalefactors_json.json",
        "_tissue_hires_image.png",
        "_tissue_lowres_image.png",
    )
    outer_member = outer.getmember(member_name)
    nested_path = target / Path(member_name).name
    source = outer.extractfile(outer_member)
    if source is None:
        raise RuntimeError(f"Unable to read outer member {member_name}")
    with source, nested_path.open("wb") as destination:
        shutil.copyfileobj(source, destination)

    with tarfile.open(nested_path, "r:gz") as nested:
        found: set[str] = set()
        for member in nested.getmembers():
            if not member.isfile() or not member.name.endswith(required_suffixes):
                continue
            source_file = nested.extractfile(member)
            if source_file is None:
                continue
            output = target / Path(member.name).name
            with source_file, output.open("wb") as destination:
                shutil.copyfileobj(source_file, destination)
            found.add(output.name)
    nested_path.unlink()
    if not any(name.endswith("_matrix.mtx.gz") for name in found):
        raise RuntimeError(f"No Matrix Market file found in {member_name}")


def find_one(directory: Path, suffix: str) -> Path:
    matches = list(directory.glob(f"*{suffix}"))
    if len(matches) != 1:
        raise RuntimeError(
            f"Expected one *{suffix} in {directory}; found {len(matches)}"
        )
    return matches[0]


def mean_available(expression: pd.DataFrame, genes: list[str]) -> pd.Series:
    available = [gene for gene in genes if gene in expression.columns]
    if not available:
        return pd.Series(np.nan, index=expression.index)
    return expression[available].mean(axis=1)


def zscore(values: pd.Series) -> pd.Series:
    values = pd.to_numeric(values, errors="coerce")
    standard_deviation = values.std(ddof=0)
    if not np.isfinite(standard_deviation) or standard_deviation <= 1e-12:
        return pd.Series(0.0, index=values.index)
    return (values - values.mean()) / standard_deviation


def tangram_celltype_fractions(
    reference: ad.AnnData,
    spatial_expression: pd.DataFrame,
    total_counts: np.ndarray,
    epochs: int,
    device: str,
) -> tuple[pd.DataFrame, dict, pd.DataFrame]:
    genes = [
        gene
        for gene in MAPPING_GENES
        if gene in spatial_expression.columns and gene in reference.var_names
    ]
    if len(genes) < 30:
        raise RuntimeError(
            f"Only {len(genes)} of {len(MAPPING_GENES)} mapping genes were found"
        )

    ad_sc = reference[:, genes].copy()
    ad_sp = ad.AnnData(
        X=spatial_expression[genes].to_numpy(dtype=np.float32),
        obs=pd.DataFrame(index=spatial_expression.index.copy()),
        var=pd.DataFrame(index=genes),
    )
    density = total_counts.astype(np.float64)
    density = density / density.sum()
    ad_sp.obs["rna_count_based_density"] = density
    ad_sp.obs["uniform_density"] = 1.0 / ad_sp.n_obs
    tg.pp_adatas(ad_sc, ad_sp, genes=genes)

    mapped = tg.map_cells_to_space(
        ad_sc,
        ad_sp,
        mode="clusters",
        cluster_label="coarse_celltype",
        density_prior="rna_count_based",
        device=device,
        num_epochs=epochs,
        learning_rate=0.1,
        scale=True,
        random_state=RANDOM_SEED,
        verbose=False,
    )

    mapping = np.asarray(mapped.X, dtype=np.float64)
    labels = mapped.obs["coarse_celltype"].astype(str).tolist()

    # Tangram's official annotation projection is mapping.T @ one_hot(labels).
    # In cluster mode each row already represents one coarse cell type, so the
    # projected spot-by-type scores are mapping.T. Normalize within spots to
    # obtain relative abundance scores; do not reapply cluster_density here.
    projected_scores = mapping.T
    denominator = projected_scores.sum(axis=1, keepdims=True)
    denominator[denominator <= 0] = 1.0
    fractions = projected_scores / denominator
    fraction_frame = pd.DataFrame(
        fractions,
        index=spatial_expression.index,
        columns=[f"fraction_{label}" for label in labels],
    )

    projected = tg.project_genes(
        mapped,
        ad_sc,
        cluster_label="coarse_celltype",
        scale=True,
    )
    observed = ad_sp[:, list(projected.var_names)].X
    projected_values = projected.X
    if hasattr(observed, "toarray"):
        observed = observed.toarray()
    if hasattr(projected_values, "toarray"):
        projected_values = projected_values.toarray()
    qc_rows = []
    for index, gene in enumerate(projected.var_names):
        rho, p_value = spearmanr(observed[:, index], projected_values[:, index])
        qc_rows.append({"gene": gene, "spearman_rho": rho, "p_value": p_value})
    qc_frame = pd.DataFrame(qc_rows)

    history = mapped.uns.get("training_history", {})
    main_history = history.get("main_loss", [])
    summary = {
        "mapping_genes": genes,
        "n_mapping_genes": len(genes),
        "epochs": epochs,
        "learning_rate": 0.1,
        "device": device,
        "random_seed": RANDOM_SEED,
        "celltype_projection": "Tangram mapping.T projection, normalized within spot",
        "final_training_objective": (
            float(main_history[-1]) if len(main_history) else None
        ),
        "median_marker_projection_rho": float(qc_frame["spearman_rho"].median()),
    }
    return fraction_frame, summary, qc_frame


def process_spatial_sample(
    sample_name: str,
    reference: ad.AnnData,
    epochs: int,
    device: str,
    force: bool,
) -> Path:
    sample_output = SAMPLE_ROOT / f"{sample_name}_spot_scores.csv.gz"
    summary_output = SAMPLE_ROOT / f"{sample_name}_mapping_summary.json"
    if sample_output.exists() and summary_output.exists() and not force:
        print(f"Skipping completed sample: {sample_name}", flush=True)
        return sample_output

    sample_info = SAMPLES[sample_name]
    print(f"Processing spatial sample: {sample_name}", flush=True)
    with tarfile.open(ARCHIVE, "r") as outer:
        with tempfile.TemporaryDirectory(
            dir=WORK_ROOT, prefix=f"{sample_name}_"
        ) as temp_name:
            temp = Path(temp_name)
            extract_required_sample_files(outer, str(sample_info["member"]), temp)
            matrix_path = find_one(temp, "_matrix.mtx.gz")
            feature_path = find_one(temp, "_features.tsv.gz")
            barcode_path = find_one(temp, "_barcodes.tsv.gz")
            position_path = find_one(temp, "_tissue_positions_list.csv")
            scale_path = find_one(temp, "_scalefactors_json.json")
            hires_path = find_one(temp, "_tissue_hires_image.png")
            lowres_path = find_one(temp, "_tissue_lowres_image.png")

            matrix = read_matrix(matrix_path)
            features = read_features(feature_path)
            barcodes = read_barcodes(barcode_path)
            expression_values, genes_present, total_counts = (
                normalized_selected_expression(matrix, features, ANALYSIS_GENES)
            )
            expression = pd.DataFrame(
                expression_values,
                index=barcodes,
                columns=genes_present,
            )
            del matrix, expression_values

            positions = pd.read_csv(
                position_path,
                header=None,
                names=[
                    "barcode",
                    "in_tissue",
                    "array_row",
                    "array_col",
                    "pixel_y_fullres",
                    "pixel_x_fullres",
                ],
            ).set_index("barcode")
            positions = positions.reindex(barcodes)
            if positions[["pixel_x_fullres", "pixel_y_fullres"]].isna().any().any():
                raise RuntimeError(f"Spatial coordinates are missing for {sample_name}")

            with scale_path.open("r", encoding="utf-8") as handle:
                scale_factors = json.load(handle)
            hires_scale = float(scale_factors.get("tissue_hires_scalef", 1.0))
            lowres_scale = float(scale_factors.get("tissue_lowres_scalef", 1.0))

            IMAGE_ROOT.mkdir(parents=True, exist_ok=True)
            shutil.copy2(hires_path, IMAGE_ROOT / f"{sample_name}_tissue_hires.png")
            shutil.copy2(lowres_path, IMAGE_ROOT / f"{sample_name}_tissue_lowres.png")
            shutil.copy2(scale_path, IMAGE_ROOT / f"{sample_name}_scalefactors.json")

            scores = pd.DataFrame(index=barcodes)
            scores.index.name = "barcode"
            scores["sample"] = sample_name
            scores["accession"] = str(sample_info["accession"])
            scores["is_primary"] = bool(sample_info["primary"])
            scores["total_counts"] = total_counts
            scores["n_genes_detected"] = np.asarray((expression > 0).sum(axis=1))
            for column in positions.columns:
                scores[column] = positions[column].to_numpy()
            scores["pixel_x_hires"] = scores["pixel_x_fullres"] * hires_scale
            scores["pixel_y_hires"] = scores["pixel_y_fullres"] * hires_scale
            scores["pixel_x_lowres"] = scores["pixel_x_fullres"] * lowres_scale
            scores["pixel_y_lowres"] = scores["pixel_y_fullres"] * lowres_scale

            scores["PKM_expression"] = expression.get("PKM", np.nan)
            scores["RPL36AL_expression"] = expression.get("RPL36AL", np.nan)
            for signature_name, genes in SIGNATURES.items():
                scores[signature_name] = mean_available(expression, genes)
            glycolysis_without_pkm = [
                gene for gene in SIGNATURES["Glycolysis_score"] if gene != "PKM"
            ]
            scores["Glycolysis_score_no_PKM"] = mean_available(
                expression, glycolysis_without_pkm
            )
            scores["Fibroblast_marker_score"] = mean_available(
                expression, CELLTYPE_MARKERS["Fibroblast"]
            )
            scores["TCDP_MIP1_score"] = mean_available(expression, MIP1_GENES)
            scores["TCDP_MIP2_score"] = mean_available(expression, MIP2_GENES)
            scores["TCDP_driver_score"] = (
                scores["TCDP_MIP1_score"] - scores["TCDP_MIP2_score"]
            )

            fractions, mapping_summary, marker_qc = tangram_celltype_fractions(
                reference,
                expression,
                total_counts,
                epochs,
                device,
            )
            scores = scores.join(fractions)
            for celltype in MAPPING_CELLTYPES:
                column = f"fraction_{celltype}"
                if column not in scores:
                    scores[column] = 0.0
            scores["Fibroblast_fraction"] = scores["fraction_Fibroblast"]
            scores["Epithelial_fraction"] = scores["fraction_Epithelial"]
            scores["Immune_fraction"] = scores[
                [
                    "fraction_T_NK",
                    "fraction_B_Plasma",
                    "fraction_Myeloid",
                    "fraction_Mast",
                ]
            ].sum(axis=1)
            scores["Other_fraction"] = scores[
                ["fraction_Endothelial", "fraction_Proliferating"]
            ].sum(axis=1)

            broad = scores[
                ["Fibroblast_fraction", "Epithelial_fraction", "Immune_fraction"]
            ]
            broad_names = broad.idxmax(axis=1).str.replace("_fraction", "", regex=False)
            scores["dominant_compartment"] = broad_names
            detailed_columns = [
                f"fraction_{celltype}" for celltype in MAPPING_CELLTYPES
            ]
            scores["dominant_celltype"] = (
                scores[detailed_columns]
                .idxmax(axis=1)
                .str.replace("fraction_", "", regex=False)
            )

            z_columns = [
                "PKM_expression",
                "RPL36AL_expression",
                "Fibroblast_fraction",
                "Epithelial_fraction",
                "Immune_fraction",
                "Fibroblast_marker_score",
                "Glycolysis_score",
                "Glycolysis_score_no_PKM",
                "Hypoxia_score",
                "IFNG_score",
                "Tcell_inflamed",
                "TCDP_driver_score",
            ]
            for column in z_columns:
                scores[f"{column}_z"] = zscore(scores[column])
            scores["PKM_fibroblast_score"] = (
                scores["PKM_expression_z"] + scores["Fibroblast_fraction_z"]
            ) / 2.0
            scores["PKM_RPL36AL_axis"] = (
                scores["PKM_expression_z"] - scores["RPL36AL_expression_z"]
            )

            SAMPLE_ROOT.mkdir(parents=True, exist_ok=True)
            scores.to_csv(sample_output, compression="gzip")
            marker_qc.to_csv(
                SAMPLE_ROOT / f"{sample_name}_tangram_marker_projection_qc.csv",
                index=False,
            )
            mapping_summary.update(
                {
                    "sample": sample_name,
                    "accession": sample_info["accession"],
                    "is_primary": sample_info["primary"],
                    "n_spots": int(scores.shape[0]),
                    "n_analysis_genes_present": len(genes_present),
                    "analysis_genes_present": genes_present,
                    "hires_scale": hires_scale,
                    "lowres_scale": lowres_scale,
                }
            )
            with summary_output.open("w", encoding="utf-8") as handle:
                json.dump(mapping_summary, handle, indent=2)
    return sample_output


def bh_adjust(p_values: Iterable[float]) -> np.ndarray:
    values = np.asarray(list(p_values), dtype=float)
    result = np.full(values.shape, np.nan, dtype=float)
    finite = np.isfinite(values)
    finite_values = values[finite]
    if finite_values.size == 0:
        return result
    order = np.argsort(finite_values)
    ranked = finite_values[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    restored = np.empty_like(adjusted)
    restored[order] = adjusted
    result[finite] = restored
    return result


def json_compatible(value):
    if isinstance(value, dict):
        return {str(key): json_compatible(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_compatible(item) for item in value]
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def signed_rank(values: Iterable[float]) -> tuple[float, int]:
    array = np.asarray(list(values), dtype=float)
    array = array[np.isfinite(array)]
    array = array[np.abs(array) > 1e-12]
    if array.size < 3:
        return math.nan, int(array.size)
    try:
        p_value = float(wilcoxon(array, alternative="two-sided", method="auto").pvalue)
    except ValueError:
        p_value = math.nan
    return p_value, int(array.size)


def add_spatial_neighborhood_columns(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    coordinates = frame[["pixel_x_fullres", "pixel_y_fullres"]].to_numpy(float)
    tree = cKDTree(coordinates)
    distances, indices = tree.query(coordinates, k=min(7, len(frame)))
    if indices.ndim == 1:
        indices = indices[:, None]
        distances = distances[:, None]
    positive_nearest = (
        distances[:, 1] if distances.shape[1] > 1 else np.ones(len(frame))
    )
    spacing = float(np.nanmedian(positive_nearest))
    maximum_distance = spacing * 1.45
    frame["neighbor_spacing_fullres"] = spacing
    frame["neighbor_distance_cutoff_fullres"] = maximum_distance

    outcomes = [
        "Glycolysis_score_no_PKM_z",
        "Hypoxia_score_z",
        "Tcell_inflamed_z",
        "IFNG_score_z",
        "Immune_fraction_z",
    ]
    neighbor_lists: list[np.ndarray] = []
    for row in range(len(frame)):
        mask = (distances[row] > 0) & (distances[row] <= maximum_distance)
        neighbor_lists.append(indices[row][mask])

    for outcome in outcomes:
        values = frame[outcome].to_numpy(float)
        neighbor_means = np.full(len(frame), np.nan)
        for row, neighbors in enumerate(neighbor_lists):
            if len(neighbors):
                neighbor_means[row] = float(np.nanmean(values[neighbors]))
        frame[f"neighbor_{outcome}"] = neighbor_means

    high_cutoff = float(frame["PKM_fibroblast_score"].quantile(0.75))
    frame["PKM_high_CAF"] = frame["PKM_fibroblast_score"] >= high_cutoff
    fibroblast_cutoff = float(frame["Fibroblast_fraction"].quantile(2 / 3))
    epithelial_cutoff = float(frame["Epithelial_fraction"].quantile(2 / 3))
    fibroblast_enriched = frame["Fibroblast_fraction"] >= fibroblast_cutoff
    epithelial_high = frame["Epithelial_fraction"].to_numpy(float) >= epithelial_cutoff
    interface = np.zeros(len(frame), dtype=bool)
    for row, neighbors in enumerate(neighbor_lists):
        interface[row] = bool(
            fibroblast_enriched.iloc[row] and np.any(epithelial_high[neighbors])
        )
    frame["Fibroblast_enriched"] = fibroblast_enriched.to_numpy()
    frame["Tumor_stroma_interface"] = interface
    return frame


def aggregate_results(sample_names: list[str]) -> None:
    frames = []
    for sample_name in sample_names:
        path = SAMPLE_ROOT / f"{sample_name}_spot_scores.csv.gz"
        frame = pd.read_csv(path)
        frame = add_spatial_neighborhood_columns(frame)
        frame.to_csv(path, index=False, compression="gzip")
        frames.append(frame)
    combined = pd.concat(frames, axis=0, ignore_index=True)
    combined.to_csv(
        RESULT_ROOT / "GSE251950_all_spot_scores.csv.gz",
        index=False,
        compression="gzip",
    )

    primary = combined[combined["is_primary"].astype(bool)].copy()
    primary_samples = [
        sample for sample in sample_names if bool(SAMPLES[sample]["primary"])
    ]

    correlation_pairs = [
        (
            "PKM_expression",
            "Fibroblast_fraction",
            "PKM vs fibroblast abundance",
            "positive",
        ),
        (
            "TCDP_driver_score",
            "Fibroblast_fraction",
            "TCDP driver vs fibroblast abundance",
            "positive",
        ),
        (
            "PKM_RPL36AL_axis",
            "Fibroblast_fraction",
            "PKM-RPL36AL axis vs fibroblast abundance",
            "positive",
        ),
        (
            "PKM_fibroblast_score",
            "Glycolysis_score_no_PKM",
            "PKM-fibroblast state vs glycolysis (PKM excluded)",
            "positive",
        ),
        (
            "PKM_fibroblast_score",
            "Hypoxia_score",
            "PKM-fibroblast state vs hypoxia",
            "positive",
        ),
        (
            "PKM_fibroblast_score",
            "Tcell_inflamed",
            "PKM-fibroblast state vs T-cell-inflamed score",
            "negative",
        ),
        (
            "PKM_fibroblast_score",
            "IFNG_score",
            "PKM-fibroblast state vs IFNG score",
            "negative",
        ),
        (
            "PKM_fibroblast_score",
            "Immune_fraction",
            "PKM-fibroblast state vs immune abundance",
            "negative",
        ),
    ]
    correlation_rows = []
    for sample_name in primary_samples:
        sample = primary[primary["sample"] == sample_name]
        for x_name, y_name, label, expected in correlation_pairs:
            rho, p_value = spearmanr(sample[x_name], sample[y_name], nan_policy="omit")
            correlation_rows.append(
                {
                    "sample": sample_name,
                    "comparison": label,
                    "x": x_name,
                    "y": y_name,
                    "expected_direction": expected,
                    "n_spots": int(sample[[x_name, y_name]].dropna().shape[0]),
                    "spearman_rho": float(rho),
                    "spot_level_p_descriptive": float(p_value),
                }
            )
    correlation_frame = pd.DataFrame(correlation_rows)
    correlation_frame.to_csv(
        RESULT_ROOT / "primary_sample_spatial_correlations.csv", index=False
    )
    correlation_summary_rows = []
    for label, group in correlation_frame.groupby("comparison", sort=False):
        p_value, n_samples = signed_rank(group["spearman_rho"])
        correlation_summary_rows.append(
            {
                "comparison": label,
                "expected_direction": group["expected_direction"].iloc[0],
                "n_primary_samples": n_samples,
                "median_spearman_rho": float(group["spearman_rho"].median()),
                "min_spearman_rho": float(group["spearman_rho"].min()),
                "max_spearman_rho": float(group["spearman_rho"].max()),
                "sample_level_wilcoxon_p": p_value,
            }
        )
    correlation_summary = pd.DataFrame(correlation_summary_rows)
    correlation_summary["BH_FDR"] = bh_adjust(
        correlation_summary["sample_level_wilcoxon_p"]
    )
    correlation_summary.to_csv(
        RESULT_ROOT / "primary_sample_spatial_correlation_summary.csv", index=False
    )

    neighborhood_outcomes = [
        ("Glycolysis_score_no_PKM_z", "glycolysis (PKM excluded)", "positive"),
        ("Hypoxia_score_z", "hypoxia", "positive"),
        ("Tcell_inflamed_z", "T-cell-inflamed score", "negative"),
        ("IFNG_score_z", "IFNG score", "negative"),
        ("Immune_fraction_z", "immune abundance", "negative"),
    ]
    neighborhood_rows = []
    direct_rows = []
    for sample_name in primary_samples:
        sample = primary[primary["sample"] == sample_name]
        high = sample["PKM_high_CAF"].astype(bool)
        for outcome, label, expected in neighborhood_outcomes:
            neighbor_column = f"neighbor_{outcome}"
            neighborhood_rows.append(
                {
                    "sample": sample_name,
                    "outcome": label,
                    "expected_direction": expected,
                    "n_high_spots": int(high.sum()),
                    "n_reference_spots": int((~high).sum()),
                    "mean_high_neighborhood": float(
                        sample.loc[high, neighbor_column].mean()
                    ),
                    "mean_reference_neighborhood": float(
                        sample.loc[~high, neighbor_column].mean()
                    ),
                    "delta_z": float(
                        sample.loc[high, neighbor_column].mean()
                        - sample.loc[~high, neighbor_column].mean()
                    ),
                }
            )
            direct_rows.append(
                {
                    "sample": sample_name,
                    "outcome": label,
                    "expected_direction": expected,
                    "mean_high_spots": float(sample.loc[high, outcome].mean()),
                    "mean_reference_spots": float(sample.loc[~high, outcome].mean()),
                    "delta_z": float(
                        sample.loc[high, outcome].mean()
                        - sample.loc[~high, outcome].mean()
                    ),
                }
            )
    neighborhood_frame = pd.DataFrame(neighborhood_rows)
    direct_frame = pd.DataFrame(direct_rows)
    neighborhood_frame.to_csv(
        RESULT_ROOT / "primary_sample_PKM_high_CAF_neighborhood_effects.csv",
        index=False,
    )
    direct_frame.to_csv(
        RESULT_ROOT / "primary_sample_PKM_high_CAF_direct_effects.csv", index=False
    )

    def summarize_effects(frame: pd.DataFrame, output_name: str) -> pd.DataFrame:
        rows = []
        for label, group in frame.groupby("outcome", sort=False):
            p_value, n_samples = signed_rank(group["delta_z"])
            rows.append(
                {
                    "outcome": label,
                    "expected_direction": group["expected_direction"].iloc[0],
                    "n_primary_samples": n_samples,
                    "mean_delta_z": float(group["delta_z"].mean()),
                    "median_delta_z": float(group["delta_z"].median()),
                    "sample_level_wilcoxon_p": p_value,
                }
            )
        summary = pd.DataFrame(rows)
        summary["BH_FDR"] = bh_adjust(summary["sample_level_wilcoxon_p"])
        summary.to_csv(RESULT_ROOT / output_name, index=False)
        return summary

    neighborhood_summary = summarize_effects(
        neighborhood_frame, "primary_sample_PKM_high_CAF_neighborhood_summary.csv"
    )
    direct_summary = summarize_effects(
        direct_frame, "primary_sample_PKM_high_CAF_direct_summary.csv"
    )

    interface_outcomes = [
        ("PKM_expression_z", "PKM expression"),
        ("PKM_fibroblast_score", "PKM-fibroblast state"),
        ("TCDP_driver_score_z", "TCDP driver score"),
    ]
    interface_rows = []
    for sample_name in primary_samples:
        sample = primary[primary["sample"] == sample_name]
        fibroblast_enriched = sample["Fibroblast_enriched"].astype(bool)
        interface = sample["Tumor_stroma_interface"].astype(bool)
        reference_mask = fibroblast_enriched & ~interface
        for outcome, label in interface_outcomes:
            if int(interface.sum()) < 5 or int(reference_mask.sum()) < 5:
                delta = math.nan
            else:
                delta = float(
                    sample.loc[interface, outcome].mean()
                    - sample.loc[reference_mask, outcome].mean()
                )
            interface_rows.append(
                {
                    "sample": sample_name,
                    "outcome": label,
                    "n_interface_spots": int(interface.sum()),
                    "n_noninterface_fibroblast_spots": int(reference_mask.sum()),
                    "delta_z": delta,
                }
            )
    interface_frame = pd.DataFrame(interface_rows)
    interface_frame.to_csv(
        RESULT_ROOT / "primary_sample_tumor_stroma_interface_effects.csv", index=False
    )
    interface_summary_rows = []
    for label, group in interface_frame.groupby("outcome", sort=False):
        p_value, n_samples = signed_rank(group["delta_z"])
        interface_summary_rows.append(
            {
                "outcome": label,
                "n_primary_samples": n_samples,
                "mean_delta_z": float(group["delta_z"].mean()),
                "median_delta_z": float(group["delta_z"].median()),
                "sample_level_wilcoxon_p": p_value,
            }
        )
    interface_summary = pd.DataFrame(interface_summary_rows)
    interface_summary["BH_FDR"] = bh_adjust(
        interface_summary["sample_level_wilcoxon_p"]
    )
    interface_summary.to_csv(
        RESULT_ROOT / "primary_sample_tumor_stroma_interface_summary.csv", index=False
    )

    sample_qc_rows = []
    for sample_name in sample_names:
        sample = combined[combined["sample"] == sample_name]
        with (SAMPLE_ROOT / f"{sample_name}_mapping_summary.json").open(
            "r", encoding="utf-8"
        ) as handle:
            mapping_summary = json.load(handle)
        sample_qc_rows.append(
            {
                "sample": sample_name,
                "accession": SAMPLES[sample_name]["accession"],
                "is_primary": SAMPLES[sample_name]["primary"],
                "n_spots": len(sample),
                "median_total_counts": float(sample["total_counts"].median()),
                "median_selected_genes_detected": float(
                    sample["n_genes_detected"].median()
                ),
                "median_marker_projection_rho": mapping_summary[
                    "median_marker_projection_rho"
                ],
                "final_training_objective": mapping_summary["final_training_objective"],
            }
        )
    pd.DataFrame(sample_qc_rows).to_csv(
        RESULT_ROOT / "GSE251950_sample_QC_summary.csv", index=False
    )

    key_summary = {
        "dataset": "GSE251950",
        "primary_samples": primary_samples,
        "sensitivity_sample": "GC6-PM" if "GC6-PM" in sample_names else None,
        "n_primary_samples": len(primary_samples),
        "n_primary_spots": int(len(primary)),
        "correlation_summary": correlation_summary.to_dict(orient="records"),
        "direct_summary": direct_summary.to_dict(orient="records"),
        "neighborhood_summary": neighborhood_summary.to_dict(orient="records"),
        "interface_summary": interface_summary.to_dict(orient="records"),
    }
    with (RESULT_ROOT / "spatial_validation_key_results.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(json_compatible(key_summary), handle, indent=2, allow_nan=False)


def write_parameters(epochs: int, device: str) -> None:
    package_versions = {
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scanpy": sc.__version__,
        "anndata": ad.__version__,
        "tangram": getattr(tg, "__version__", "1.0.4"),
    }
    parameters = {
        "dataset": "GSE251950",
        "primary_analysis_samples": [
            sample for sample, info in SAMPLES.items() if bool(info["primary"])
        ],
        "sensitivity_sample": "GC6-PM",
        "single_cell_reference": "GSE163558 primary gastric tumors PT1/PT2/PT3",
        "single_cell_reference_samples": PRIMARY_REFERENCE_SAMPLES,
        "normalization": "library-size normalization to 10000 counts followed by log1p",
        "mapping_method": "Tangram cluster mode using fixed MultiSTAD coarse-cell marker genes",
        "mapping_genes": MAPPING_GENES,
        "mapping_celltypes": MAPPING_CELLTYPES,
        "tangram_epochs": epochs,
        "tangram_learning_rate": 0.1,
        "tangram_density_prior": "rna_count_based",
        "device": device,
        "random_seed": RANDOM_SEED,
        "celltype_markers": CELLTYPE_MARKERS,
        "signatures": SIGNATURES,
        "MIP1_driver_genes": MIP1_GENES,
        "MIP2_driver_genes": MIP2_GENES,
        "TCDP_driver_score": "mean(MIP1-high genes) - mean(MIP2-high genes)",
        "PKM_fibroblast_score": "mean(within-section z(PKM), within-section z(Tangram fibroblast fraction))",
        "PKM_high_CAF_definition": "top within-section quartile of PKM_fibroblast_score",
        "neighborhood_definition": "up to six nearest Visium spots within 1.45x median nearest-neighbor distance; index spot excluded",
        "interface_definition": "upper-tertile fibroblast-fraction spot with at least one upper-tertile epithelial-fraction neighbor",
        "sample_level_inference": "two-sided Wilcoxon signed-rank test across nine primary sections, BH-adjusted within analysis family",
        "package_versions": package_versions,
    }
    RESULT_ROOT.mkdir(parents=True, exist_ok=True)
    with (RESULT_ROOT / "analysis_parameters.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(parameters, handle, indent=2)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--samples",
        nargs="+",
        default=list(SAMPLES),
        choices=list(SAMPLES),
        help="Spatial samples to process; defaults to all 10 sections.",
    )
    parser.add_argument("--epochs", type=int, default=1000)
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--rebuild-reference", action="store_true")
    parser.add_argument("--skip-aggregate", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    np.random.seed(RANDOM_SEED)
    for directory in [WORK_ROOT, RESULT_ROOT, IMAGE_ROOT, SAMPLE_ROOT, REFERENCE_ROOT]:
        directory.mkdir(parents=True, exist_ok=True)
    if not ARCHIVE.exists():
        raise FileNotFoundError(f"Missing spatial archive: {ARCHIVE}")

    reference = build_single_cell_reference(force=args.rebuild_reference)
    write_parameters(args.epochs, args.device)
    for sample_name in args.samples:
        process_spatial_sample(
            sample_name,
            reference,
            epochs=args.epochs,
            device=args.device,
            force=args.force,
        )
    if not args.skip_aggregate:
        missing = [
            sample
            for sample in args.samples
            if not (SAMPLE_ROOT / f"{sample}_spot_scores.csv.gz").exists()
        ]
        if missing:
            raise RuntimeError(f"Missing processed samples: {missing}")
        aggregate_results(args.samples)
    print("Spatial association analysis completed.", flush=True)


if __name__ == "__main__":
    main()
