# MultiSTAD

MultiSTAD evaluates a fixed host transcriptional program across independent gastric cancer cohorts and cellular, spatial and histological contexts. The signed 20-gene program was historically selected from joint host–microbial MIP states. The current workflow focuses on reproducible metabolic and partial immune RNA associations, specificity controls and complementary multiscale evidence.

Associated manuscript: *Cross-cohort metabolic RNA associations and multiscale tissue correlates of a fixed host transcriptional program in gastric cancer*.

[Study overview (PDF)](assets/Figure1.pdf)

## Highlights

- Fixed signed program with no outcome-dependent gene reselection.
- Four GEO cohorts and CPTAC, with PKM omission, shared-gene removal and matched random-program controls.
- Patient-level single-cell and spatial comparisons.
- Frozen histology predictions and tumor-only exploratory immunotherapy associations.

## Repository structure

```text
MultiSTAD/
├── main.py                     # Main analysis entrypoint
├── src/
│   ├── preprocessing.py        # Historical CLR and training-derived scaling
│   ├── program.py             # Signed-score and differential-expression helpers
│   ├── selection.py           # Historical selection; not run by main.py
│   ├── cohorts.py             # GEO/CPTAC loading and cross-cohort associations
│   ├── bulk.py                # Fixed-program specificity analyses
│   ├── specificity.py         # Paired seven-gene glycolysis comparison
│   ├── single_cell.py         # Patient-level cell-type contrasts and QC
│   ├── spatial.py             # Patient-level spatial effects
│   ├── spatial_mapping.py     # Optional original Tangram workflow
│   ├── spatial_sensitivity.py # Original fixed spatial-region sensitivities
│   ├── features.py            # Distribution summaries of frozen UNI features
│   ├── histology.py           # Mean-OOF metrics and original nested-CV function
│   ├── immunotherapy.py       # Tumor-only association and paired comparators
│   ├── virtual_knockout_*.R   # Supplementary computational perturbation
│   └── statistics.py          # Shared resampling and statistical utilities
├── data/                      # Input specification and fixed gene definitions
├── assets/                    # Study overview
├── docs/methods.md             # Analysis settings and reproduction notes
└── requirements.txt
```

## Installation

```bash
conda create -n multistad python=3.12 -y
conda activate multistad
pip install -r requirements.txt
```

Optional Tangram and R dependencies are listed in [methods](docs/methods.md). The repository excludes raw datasets, processed patient records, UNI embeddings, WSI files, trained weights, logs and result tables. Prepare the inputs described in [data/README.md](data/README.md); data can remain outside the repository.

## Usage

```bash
python main.py --list
python main.py --stage cohorts --data-root /path/to/data --output /path/to/results/cohort_run
python main.py --stage bulk --data-root /path/to/data --output /path/to/results/bulk_run
python main.py --stage single_cell --data-root /path/to/data --output /path/to/results/cell_run
python main.py --stage spatial --data-root /path/to/data --output /path/to/results/spatial_run
python main.py --stage histology --data-root /path/to/data --output /path/to/results/histology_run
python main.py --stage immunotherapy --data-root /path/to/data --output /path/to/results/icb_run
```

Each output directory must be new and separate from the data directory. `bulk` also runs the prespecified simple-glycolysis comparison; `single_cell` runs both the primary and QC sensitivity passes. `histology` evaluates existing patient-mean OOF predictions and does not retrain a model. `cohorts` provides the original external-cohort loading and association workflow.

## Inputs and execution scope

The repository contains analysis source, frozen gene/signature definitions and the overview PDF. Public raw datasets must be obtained from their original repositories. Historical processed matrices, cell annotations, spatial mappings, OOF predictions and tumor-level immunotherapy scores are **not distributed**. The paper has one Supplementary Materials document and does not supply a separate Supplementary Data CSV archive.

`python main.py --list` works without scientific datasets. None of the six analysis stages is a bundled-data demonstration or a one-command raw-data reproduction of the manuscript.

| Stage | Required inputs | Outputs under the new `--output` directory | Scope |
|---|---|---|---|
| `cohorts` | GEO matrices/platform annotations, Hallmark GMT, CPTAC RNA/protein matrices | `cohorts/` tables, scores and analysis figures | External inputs required; historical loading and association implementation |
| `bulk` | All-gene cohort expression, frozen Discovery split/expression, CPTAC matching/metadata, GMT | `bulk/` specificity and sensitivity tables; `logs/` | External inputs required; no new gene selection |
| `single_cell` | Frozen cell expression panels, labels, scores, QC and GEO metadata | `single_cell/` patient contrasts and QC tables | External intermediate inputs required; does not reconstruct cell labels |
| `spatial` | Existing spot/section effects, QC, GEO metadata and saved gate decision | `spatial/` nine-patient effect tables | Archived output evaluation; does not perform Tangram mapping |
| `histology` | 2,640 archived OOF predictions for 264 patients | `histology/` patient-mean predictions and metrics | Archived output evaluation; no training or WSI processing |
| `immunotherapy` | Frozen scores for 45 tumors; optional archived comparator intervals | `immunotherapy/` association and paired-comparison tables | Archived score evaluation; does not reconstruct upstream normalization |

Use `main.py` for these protected entrypoints. Direct execution/import of historical analysis modules can read inputs or create output directories and does not provide the entrypoint's overwrite protection. Historical training, mapping and R scripts require separately prepared inputs/environments and are not run by the six-stage interface. See [data preparation](data/README.md) and [methods](docs/methods.md).

## Interpretation and reproducibility

The code is organized from the actual analysis implementations; path handling and entrypoints were adapted for this repository. The fixed program is observational and has not established consistent superiority over general metabolic controls or independent CAF specificity. Histology performance is internal OOF evaluation; ICB results are exploratory in 45 tumor patients. Historical raw inputs and some exact software versions are not included or fully recovered, so the package does not claim complete end-to-end reproduction. See [methods](docs/methods.md) for the applicable units, settings and remaining limitations.

## Citation

Please cite the associated manuscript when it becomes available. This source package does not imply redistribution rights for third-party data or model weights.
