"""Current manuscript analysis entrypoint. Data are external; existing outputs are protected."""

import argparse, os, runpy, sys
from pathlib import Path

STAGES = {
    "cohorts": "GEO/CPTAC loading, fixed-program associations and meta-analysis",
    "bulk": "Fixed-program ablations, overlap removal, covariates and matched controls",
    "single_cell": "Patient-level paired contrasts and prespecified QC sensitivities",
    "spatial": "Nine-patient effects from existing fixed spatial region outputs",
    "histology": "Patient-mean OOF performance; no model training",
    "immunotherapy": "45-tumor exploratory association and five paired comparisons",
}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--list", action="store_true")
    p.add_argument("--stage", choices=STAGES)
    p.add_argument(
        "--data-root",
        type=Path,
        default=Path(os.environ.get("MULTISTAD_DATA_ROOT", "data")),
    )
    p.add_argument(
        "--output",
        type=Path,
        help="New output directory; existing directories are rejected",
    )
    a = p.parse_args()
    if a.list or not a.stage:
        for name, description in STAGES.items():
            print(f"{name:16} {description}")
        return
    if a.output is None:
        p.error("--output is required")
    data = a.data_root.resolve()
    out = a.output.resolve()
    if not data.is_dir():
        p.error("Data root does not exist; see data/README.md")
    if out == data or out in data.parents or data in out.parents:
        p.error("Output must be separate from the input data tree")
    if out.exists():
        p.error("Output already exists; choose a new directory")
    os.environ["MULTISTAD_DATA_ROOT"] = str(data)
    os.environ["MULTISTAD_OUTPUT_ROOT"] = str(out)
    out.mkdir(parents=True)
    try:
        if a.stage == "histology":
            from multistad.histology import evaluate_oof

            print(
                evaluate_oof(data / "histology/oof_predictions.csv", out / "histology")
            )
        elif a.stage == "immunotherapy":
            from multistad.immunotherapy import evaluate_tumors

            archived = data / "immunotherapy/comparators.csv"
            result, _ = evaluate_tumors(
                data / "immunotherapy/tumor_scores.csv",
                archived if archived.is_file() else None,
                out / "immunotherapy",
            )
            print(result)
        elif a.stage == "single_cell":
            for mode in ["main", "sensitivity"]:
                sys.argv = ["single_cell", mode]
                runpy.run_module("multistad.single_cell", run_name="__main__")
        elif a.stage == "bulk":
            runpy.run_module("multistad.bulk", run_name="__main__")
            runpy.run_module("multistad.specificity", run_name="__main__")
        else:
            runpy.run_module("multistad." + a.stage, run_name="__main__")
    except Exception:
        (out / "INCOMPLETE.txt").write_text(
            "Run failed. Outputs are incomplete; inspect the raised error.\n"
        )
        raise
    (out / "COMPLETED.txt").write_text(a.stage + "\n")


if __name__ == "__main__":
    main()
