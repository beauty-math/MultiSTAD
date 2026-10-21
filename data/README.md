# Data preparation

Only `definitions/program.csv` and `definitions/signatures.json` are distributed here. They contain the frozen signed modules, common core and target signatures. Supply original or audited intermediate inputs externally; do not substitute example data for manuscript observations.

```text
DATA_ROOT/
├── bulk/
│   ├── GSE62254_expression.csv.gz       # genes × expression records
│   ├── GSE15459_expression.csv.gz
│   ├── GSE84433_expression.csv.gz
│   ├── GSE26253_expression.csv.gz
│   ├── CPTAC_RNA.cct                   # genes × sample IDs; tab separated
│   ├── CPTAC_meta.tsv                  # original metatable, including CAT row
│   └── CPTAC_matching.csv              # sample_id and case_id
├── discovery/
│   ├── expression.csv                 # cases × genes, original unscaled matrix
│   └── split.csv                      # case_id, split; train/validation frozen
├── external/                          # original GEO/CPTAC download layout
│   ├── bulk/<GEO>/<GEO>_series_matrix.txt.gz
│   ├── annotation/                    # original GPL570/GPL6947/GPL8432 annotation
│   ├── signature/h.all.v2026.1.Hs.symbols.gmt
│   └── CPTAC_STAD/                     # original RNA and protein .cct files
├── metadata/
│   ├── GSE183904_official.soft.txt
│   ├── GSE163558_official.soft.txt
│   └── GSE251950_official.soft.txt
├── single_cell/
│   ├── <GEO>_panel_expression.csv
│   ├── <GEO>_coarse_celltype_metadata.csv
│   ├── <GEO>_cell_scores.csv
│   ├── <GEO>_full_QC.csv.gz
│   └── raw/GSE163558/                  # optional original mapping reference
├── spatial/
│   ├── gate.json                      # saved decision: GO or CONDITIONAL GO
│   ├── GSE251950_sample_QC_summary.csv
│   ├── GSE251950_all_spot_scores.csv.gz # optional sensitivity input
│   ├── samples/GC*_spot_scores.csv.gz
│   ├── sensitivity/PKM_spatial_sensitivity_sample_effects.csv
│   ├── sensitivity/PKM_spatial_sensitivity_summary.csv
│   └── raw/GSE251950_RAW.tar           # optional original mapping input
├── histology/
│   ├── oof_predictions.csv            # 2640 rows; case_id, repeat,
│   │                                  # measured_ecosystem_score,
│   │                                  # locked_outer_fold_prediction
│   ├── features/                      # optional original train/val_dist_features.csv
│   └── targets/                       # optional original step25 target matrices
└── immunotherapy/
    ├── tumor_scores.csv               # 45 unique patient_id; response_binary;
    │                                  # ecosystem_score; PKM, IFNG,
    │                                  # Tcell_inflamed, CD8, Merck18
    └── comparators.csv                # optional archived delta-AUC CIs and tests
```

The filenames and column contracts follow the archived analysis. Histology and immunotherapy analyses require archived out-of-fold predictions and tumor-level scores, respectively. These intermediate inputs are not distributed with the code repository. The optional comparator table uses `program`, `estimate`, `p_value`, `ci_low` and `ci_high`; a mismatch with recomputed point estimates raises an error. Without this table, comparator CIs remain NA.

Bulk arrays require the original target-aware probe mapping, sample exclusions and transforms; `cohorts.py` preserves the original loader. The specificity inputs additionally contain the all-gene pool used for random controls, so targeted scores alone cannot replace expression matrices. Single-cell matrices and labels must be the original frozen objects; no new cell annotation is inferred by the main entrypoint.

The tree above is an **input contract**, not a list of files shipped in this repository. Only the two files under `definitions/` are supplied. Public downloads and historical intermediates are distinct: GEO series/platform records, GDC materials and the original CPTAC release are external sources; the exact frozen split, processed cell objects, Tangram outputs, OOF records and tumor-level score tables require archived project inputs. Public raw-data access alone does not provide these intermediates.

The exact CPTAC filenames expected by `cohorts` are `GC_RNAseq_gene_rsem_rm_circRNA_UQ_log2(x+1)_Tumor.cct` and `GC_TMT_proteomics_ratio_gene_MD_GencodeV42_Tumor.cct`. The annotation files are `GPL570.annot.gz`, `GPL6947.annot.gz` and `GPL8432_HUMANREF-8_V3_0_R1_11282963_A_WGDASL.txt.gz`. `cohorts` produces targeted gene scores; it does not automatically construct the all-gene `bulk` specificity inputs. Stage inputs are prepared separately, not inferred from preceding output folders.

Sources are GEO GSE62254, GSE15459, GSE84433, GSE26253, GSE183904, GSE163558 and GSE251950; the original CPTAC gastric multi-omics release; TCGA-STAD/GDC; and PRJEB25780. Obtain data under their respective access conditions. For optional R analyses, retain the historical `scRNA/` and `processed_v2/step13_scrna_celltype_anchor/` input layout under `MULTISTAD_DATA_ROOT`. The manuscript is accompanied by one Supplementary Materials document; no separate Supplementary Data CSV archive is supplied. Historical patient-level intermediate files are not distributed here.
