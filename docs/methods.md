# Analysis settings

## Historical discovery and fixed program

Original MIP states were discovered jointly from host and microbial features. They are frozen inputs to the current study, not newly inferred independent microbial states. Microbial features were filtered on Discovery zero fraction (retain ≤0.80), transformed as `log(x+1) − mean_taxa(log(x+1))` within each sample, then z-scored per feature using Discovery means and population standard deviations. The internal partition reused those parameters. Thus the clustering inputs had **CLR followed by training-derived z-scoring**, whereas descriptive Bray–Curtis analysis used relative abundances from nonnegative counts.

The frozen historical data retained 191 of 250 microbial features. Re-evaluation using the archived 211/53 split reproduced the saved 211 × 191 and 53 × 191 transformed matrices to floating-point precision. The number 191 is an observed property of those inputs, not a hard-coded output size for arbitrary data. `preprocessing.py` preserves the transform; the complete original MIP clustering and plotting workflow is not included in this curated package.

The historical gene-selection criteria were parallel: 19531 tested genes; 1188 FDR-qualified; 2614 effect-size-qualified; 974 satisfying both. Measurement availability in the external PRJEB data and ranking reduced candidates to 300. One hundred 80% draws without replacement, L1 logistic regression (liblinear, C=0.2, balanced classes, seed 2026), selection frequency and coefficient magnitude yielded the fixed 20 genes. `selection.py` preserves the historical algorithm, including its original fallback, but is not called by the current entrypoint. No new selection was performed for this repository.

Bulk scores are mean available fixed positive-gene z minus mean available fixed negative-gene z, with cohort-specific gene standardization (ddof=0). Selection coefficients are not score weights. The historical 11-gene core and directions remain fixed. `program.py` retains historical helpers; `bulk.py` enforces the actual variant/coverage rules used for the specificity analysis.

## Independent cohorts and specificity

Four GEO cohorts contribute 300, 192, 357 and 432 expression records under the historical one-record-per-patient assumption; cross-study patient deduplication was not established. CPTAC has 151 unique cases. Available-gene coverage remains cohort-specific.

`bulk.py` preserves PKM omission, target-overlap removal, limited covariate controls and 200 matched random programs per cohort (seed 20261011). Matching uses expression/detection/variance ranks and Discovery detection, preserves measured module sizes and does not use target association results. The empirical comparison is `(1 + count(|null rho| ≥ |observed rho|))/201`; null percentiles are not patient-effect CIs.

Patient/record-paired bootstrap uses 5000 draws. Association BH families contain 25 program–target tests per cohort; omission deltas use the 20 nonreference contrasts per cohort. Overlap-control q refers to the post-removal correlation, not its paired delta; no delta-q was archived. The 25 random-program tests form one BH family. Seven-gene comparisons use HK2/PFKP/ALDOA/GAPDH/ENO1/LDHA/SLC2A1 against four targets, with four-test BH within cohort. Delta is the difference in **signed** correlations, not absolute association strength.

## Single-cell and spatial context

The paired primary single-cell comparisons use 25 patients in GSE183904 and 2 in GSE163558, with at least 20 cells per sample/type before aggregation. Labels are frozen marker-panel inferences, not recovered author gold standards. The second-cohort original-program fibroblast–epithelial contrast is negative. Counts, QC residuals, balanced-cell sensitivities and all comparison directions remain in the analysis.

Spatial inference uses GC1–GC9 primary patients, excluding GC6-PM; 29808 spots are measurements. The main region is the upper quartile of mean section z(PKM) and z(mapped fibroblast fraction). The historical sensitivity takes PKM's upper quartile within the mapped-fibroblast upper third. These definitions must not be exchanged. There are 60 comparisons per QC stratum (6 definitions × 2 analysis types × 5 outcomes), 120 saved rows overall. `spatial.py` consumes existing region/section outputs and the original gate; it does not remap cells. The optional mapping and sensitivity modules preserve the original upstream implementations.

Recorded mapping versions: NumPy 1.26.4, pandas 2.3.3, scanpy 1.11.5, anndata 0.11.4, Tangram 1.0.4. Use a separate compatible environment when invoking those optional modules. `MULTISTAD_DATA_ROOT` and `MULTISTAD_OUTPUT_ROOT` configure their paths. Both optional R scripts directly load Matrix, data.table, matrixStats, ggplot2, ggrepel, dplyr, jsonlite and scTenifoldKnk; historical exact R/package versions were not recovered. The selected computational KO runs were pooledAuto for PKM and pooledAuto2 for RPL36AL; other historical script outputs are not additional current claims.

The optional Python mapping imports `anndata`, `scanpy` and `tangram` (distribution name `tangram-sc`), in addition to NumPy, pandas and SciPy. Tangram's PyTorch dependency must be installed in that separate environment; no UNI weights or feature encoder are bundled. These optional dependencies are not required for the six-stage core interface. R base/recommended packages and transitive dependencies are provided by the R installation/package manager; the eight packages above are the complete explicit library/namespace dependencies found in the two scripts. Neither R network analysis nor Tangram mapping was rerun for this release check.

## Histology and immunotherapy

UNI embeddings remain frozen. `features.py` computes mean, standard deviation and five quantiles per 1024 dimensions plus two patch-count descriptors (7170 features). `histology.py` retains the actual 5 × 10 outer/5-fold inner Ridge procedure for method transparency; the main entrypoint only evaluates saved OOF predictions. Each of 264 patients contributes the mean of 10 predictions. Prediction R² is 1 − SSE/SST. Historical 5000-patient bootstrap uses seed 20261006; default evaluation reports point estimates without substituting r² for prediction R².

ICB analysis uses 45 unique tumor patients (12 responders / 33 nonresponders), excluding 33 normal records from the original 78. Use tumor-only-normalized `ecosystem_score`, not legacy `driver_score`. The univariable OR is per population SD. Five paired DeLong comparisons use Holm correction; archived paired-bootstrap delta-AUC CIs are carried only when their point estimates/tests match. No threshold/model optimization is performed. Exact upstream TIGER definitions and the individual molecular-subtype mapping remain incompletely recovered.

## Reproduction scope

This is a curated source release, not an archive of every intermediate output. Core algorithms were retained or extracted from actual analysis scripts; paths, module imports and entrypoints were adapted. Current-code syntax and selected frozen-data regression checks were performed locally. Large raw-data stages, training, Tangram mapping and R network inference were not rerun for packaging. Historical input/environment gaps prevent a claim of complete end-to-end reproduction. The overview PDF is provided for orientation; final manually edited figure production is not an additional analysis module.

The bundled `requirements.txt` defines the core Python environment. Historical dependency versions are not universally recoverable; a successful interface check under current compatible versions is not proof of bitwise reproduction of every historical analysis. Frozen score definitions and numerical algorithms remain unchanged. Legacy variable names, output filenames and data columns such as `driver_score` and `proteomic_validation` are retained for compatibility; descriptive titles refer to associations, not mechanistic or clinical validation.
