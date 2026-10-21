"""Historical selection only; never invoked by the current main workflow.
Inputs are frozen Discovery expression/labels and historical PRJEB measurability.
Old gene-microbe networks and 78-record ICB evaluation are excluded.
"""

import os
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from .program import diff_table_numeric


def historical_selection(bulk_train, labels, prjeb_expr, output_dir):
    OUT_DIR = os.fspath(output_dir)
    os.makedirs(OUT_DIR, exist_ok=False)
    gene_diff = diff_table_numeric(
        bulk_train, labels, "MIP-1", "MIP-2", min_abs_effect=0.30
    )
    gene_sig = gene_diff[gene_diff["fdr"] < 0.05].copy()
    shared_genes = sorted(
        set(gene_sig["feature"]).intersection(set(prjeb_expr.columns))
    )
    candidate_genes = gene_sig[gene_sig["feature"].isin(shared_genes)].copy()

    candidate_genes = (
        candidate_genes.sort_values(["fdr", "abs_effect"], ascending=[True, False])
        .head(300)
        .copy()
    )
    candidate_gene_list = candidate_genes["feature"].tolist()

    pd.DataFrame({"candidate_gene": candidate_gene_list}).to_csv(
        os.path.join(OUT_DIR, "candidate_driver_genes_shared_with_prjeb.csv"),
        index=False,
    )

    print("Candidate driver genes shared with PRJEB:", len(candidate_gene_list))

    # =========================================================
    # 7. 稳定特征选择
    # =========================================================
    y_train = (labels == "MIP-1").astype(int).copy()

    X_train = bulk_train[candidate_gene_list].copy()
    for c in X_train.columns:
        X_train[c] = pd.to_numeric(X_train[c], errors="coerce")
    X_train = X_train.fillna(X_train.median(axis=0))

    N_BOOT = 100
    SUBSAMPLE_FRAC = 0.80
    C_VALUE = 0.2

    selection_count = pd.Series(0, index=X_train.columns, dtype=float)
    coef_sum = pd.Series(0.0, index=X_train.columns, dtype=float)

    rng = np.random.RandomState(2026)

    for b in range(N_BOOT):
        idx = rng.choice(
            X_train.index, size=int(len(X_train) * SUBSAMPLE_FRAC), replace=False
        )
        X_sub = X_train.loc[idx].copy()
        y_sub = y_train.loc[idx].copy()

        scaler = StandardScaler()
        X_sub_z = scaler.fit_transform(X_sub)

        clf = LogisticRegression(
            penalty="l1",
            solver="liblinear",
            C=C_VALUE,
            class_weight="balanced",
            max_iter=2000,
            random_state=2026 + b,
        )
        clf.fit(X_sub_z, y_sub)
        coef = pd.Series(clf.coef_.flatten(), index=X_train.columns)

        selected = coef[coef != 0].index.tolist()
        selection_count.loc[selected] += 1
        coef_sum += coef.abs()

    stability_df = (
        pd.DataFrame(
            {
                "gene": X_train.columns,
                "selection_freq": selection_count.values / N_BOOT,
                "mean_abs_coef_proxy": coef_sum.values / N_BOOT,
            }
        )
        .sort_values(
            ["selection_freq", "mean_abs_coef_proxy"], ascending=[False, False]
        )
        .reset_index(drop=True)
    )

    stability_df.to_csv(
        os.path.join(OUT_DIR, "driver_gene_stability_selection.csv"), index=False
    )

    selected_driver = stability_df[stability_df["selection_freq"] >= 0.60].copy()

    if len(selected_driver) < 6:
        selected_driver = stability_df.head(10).copy()

    if len(selected_driver) > 20:
        selected_driver = selected_driver.head(20).copy()

    selected_driver.to_csv(
        os.path.join(OUT_DIR, "selected_driver_genes.csv"), index=False
    )

    selected_gene_list = selected_driver["gene"].tolist()
    print("Selected driver genes:", len(selected_gene_list))

    # =========================================================
    # 8. 方向定义
    # =========================================================
    gene_dir = candidate_genes[
        candidate_genes["feature"].isin(selected_gene_list)
    ].copy()
    gene_dir["higher_in"] = np.where(gene_dir["effect_mean_diff"] > 0, "MIP-1", "MIP-2")

    up_genes = gene_dir.loc[gene_dir["higher_in"] == "MIP-1", "feature"].tolist()
    down_genes = gene_dir.loc[gene_dir["higher_in"] == "MIP-2", "feature"].tolist()

    with open(os.path.join(OUT_DIR, "driver_up_genes_MIP1.txt"), "w") as f:
        f.write("\n".join(up_genes))

    with open(os.path.join(OUT_DIR, "driver_up_genes_MIP2.txt"), "w") as f:
        f.write("\n".join(down_genes))

    return selected_driver, up_genes, down_genes
