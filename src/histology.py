"""Internal out-of-fold histological evaluation and original nested-CV function; no training in default workflow."""

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error
from sklearn.model_selection import RepeatedKFold, KFold, GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.feature_selection import VarianceThreshold
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
from .settings import DATA, OUTPUT

OUT = OUTPUT / "histology"
SEED = 20261006
N_BOOT = 5000


def metric_values(y, pred):
    y = np.asarray(y, float)
    pred = np.asarray(pred, float)
    return {
        "pearson_r": float(stats.pearsonr(y, pred).statistic),
        "spearman_rho": float(stats.spearmanr(y, pred).statistic),
        "r2": float(r2_score(y, pred)),
        "mae": float(mean_absolute_error(y, pred)),
        "rmse": float(mean_squared_error(y, pred) ** 0.5),
    }


def bootstrap_metrics(y, pred, n_boot=N_BOOT, seed=SEED):
    rng = np.random.default_rng(seed)
    y = np.asarray(y, float)
    pred = np.asarray(pred, float)
    names = list(metric_values(y, pred))
    values = {name: [] for name in names}
    for _ in range(n_boot):
        idx = rng.integers(0, len(y), len(y))
        if np.std(y[idx]) == 0 or np.std(pred[idx]) == 0:
            continue
        row = metric_values(y[idx], pred[idx])
        for name in names:
            values[name].append(row[name])
    rows = []
    point = metric_values(y, pred)
    for name in names:
        arr = np.asarray(values[name])
        rows.append(
            {
                "metric": name,
                "estimate": point[name],
                "ci_low": float(np.quantile(arr, 0.025)),
                "ci_high": float(np.quantile(arr, 0.975)),
                "bootstrap_valid": int(len(arr)),
            }
        )
    return pd.DataFrame(rows)


def pathology_nested_cv():
    OUT.mkdir(parents=True, exist_ok=True)
    feature_dir = DATA / "histology/features"
    target_dir = DATA / "histology/targets"
    x_train = pd.read_csv(feature_dir / "train_dist_features.csv", index_col=0)
    x_val = pd.read_csv(feature_dir / "val_dist_features.csv", index_col=0)
    y_train = pd.read_csv(target_dir / "step25_target_matrix_train.csv", index_col=0)[
        "driver_score"
    ]
    y_val = pd.read_csv(
        target_dir / "step25_target_matrix_validation.csv", index_col=0
    )["driver_score"]
    x = pd.concat([x_train, x_val], axis=0)
    y = pd.concat([y_train, y_val], axis=0).reindex(x.index)
    if len(x) != 264 or y.isna().any() or not x.index.is_unique:
        raise RuntimeError("Unexpected pathology matrix alignment")

    alphas = np.logspace(-2, 6, 13)
    outer = RepeatedKFold(n_splits=5, n_repeats=10, random_state=SEED)
    prediction_rows = []
    selection_rows = []
    x_values = x.to_numpy(dtype=np.float32)
    y_values = y.to_numpy(float)
    for split_index, (train_idx, test_idx) in enumerate(outer.split(x_values)):
        repeat = split_index // 5 + 1
        fold = split_index % 5 + 1
        inner = KFold(n_splits=5, shuffle=True, random_state=SEED + split_index)
        pipe = Pipeline(
            [
                ("impute", SimpleImputer(strategy="median")),
                ("variance", VarianceThreshold()),
                ("scale", StandardScaler()),
                ("ridge", Ridge(solver="lsqr", max_iter=20000)),
            ]
        )
        grid = GridSearchCV(
            pipe,
            {"ridge__alpha": alphas},
            scoring="neg_mean_squared_error",
            cv=inner,
            n_jobs=4,
            refit=True,
            error_score="raise",
        )
        grid.fit(x_values[train_idx], y_values[train_idx])
        pred = grid.predict(x_values[test_idx])
        for idx, value in zip(test_idx, pred):
            prediction_rows.append(
                {
                    "case_id": x.index[idx],
                    "repeat": repeat,
                    "outer_fold": fold,
                    "measured_ecosystem_score": y_values[idx],
                    "locked_outer_fold_prediction": float(value),
                }
            )
        selection_rows.append(
            {
                "repeat": repeat,
                "outer_fold": fold,
                "n_outer_train": len(train_idx),
                "n_outer_test": len(test_idx),
                "selected_alpha": float(grid.best_params_["ridge__alpha"]),
                "inner_cv_neg_mse": float(grid.best_score_),
            }
        )

    long = pd.DataFrame(prediction_rows)
    long.to_csv(OUT / "pathology_nested_cv_predictions_long.csv", index=False)
    pd.DataFrame(selection_rows).to_csv(
        OUT / "pathology_nested_cv_inner_selection.csv", index=False
    )
    aggregated = long.groupby("case_id", as_index=False).agg(
        measured_ecosystem_score=("measured_ecosystem_score", "first"),
        morphology_ecosystem_score=("locked_outer_fold_prediction", "mean"),
        prediction_sd_across_repeats=("locked_outer_fold_prediction", "std"),
        n_locked_predictions=("locked_outer_fold_prediction", "size"),
    )
    aggregated.to_csv(
        OUT / "pathology_nested_cv_patient_oof_predictions.csv", index=False
    )
    summary = bootstrap_metrics(
        aggregated["measured_ecosystem_score"], aggregated["morphology_ecosystem_score"]
    )
    summary.insert(0, "n_patients", len(aggregated))
    summary.to_csv(OUT / "pathology_nested_cv_metrics_bootstrap.csv", index=False)

    per_repeat = []
    for repeat, frame in long.groupby("repeat"):
        m = metric_values(
            frame["measured_ecosystem_score"], frame["locked_outer_fold_prediction"]
        )
        per_repeat.append({"repeat": repeat, **m})
    pd.DataFrame(per_repeat).to_csv(
        OUT / "pathology_nested_cv_metrics_by_repeat.csv", index=False
    )


def evaluate_oof(path, output):
    x = pd.read_csv(path)
    required = [
        "case_id",
        "repeat",
        "measured_ecosystem_score",
        "locked_outer_fold_prediction",
    ]
    if (
        x[required].isna().any().any()
        or not np.isfinite(x[required[2:]].to_numpy()).all()
    ):
        raise ValueError("Missing/nonfinite patient identity, repeat or frozen score")
    if (
        len(x) != 2640
        or x.case_id.nunique() != 264
        or x.duplicated(["case_id", "repeat"]).any()
    ):
        raise ValueError(
            "Expected 264 patients and 2640 unique patient-repeat OOF records"
        )
    if (
        not (x.groupby("case_id").size() == 10).all()
        or not (x.groupby("case_id").measured_ecosystem_score.nunique() == 1).all()
    ):
        raise ValueError("Inconsistent patient scores or repeat counts")
    p = x.groupby("case_id").agg(
        observed=("measured_ecosystem_score", "first"),
        predicted=("locked_outer_fold_prediction", "mean"),
    )
    result = metric_values(p.observed, p.predicted)
    output.mkdir(parents=True, exist_ok=True)
    p.to_csv(output / "patient_mean_oof.csv")
    pd.DataFrame([dict(n=264, **result)]).to_csv(output / "metrics.csv", index=False)
    return result
