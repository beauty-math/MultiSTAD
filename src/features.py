"""Distributional UNI feature calculation from the original feature-generation code."""

import numpy as np

QUANTILES = [10, 25, 50, 75, 90]


def distribution_features_from_bag(bag, mode="dist"):
    bag = np.asarray(bag, dtype=np.float32)
    n_patches, dim = bag.shape
    mean_feat = np.mean(bag, axis=0)
    if mode == "mean":
        vals = [mean_feat]
        names = [f"mean_{i}" for i in range(dim)]
    elif mode == "dist":
        std_feat = np.std(bag, axis=0)
        q_feats = np.percentile(bag, QUANTILES, axis=0)
        vals = [mean_feat, std_feat] + [q_feats[i, :] for i in range(len(QUANTILES))]
        names = [f"mean_{i}" for i in range(dim)] + [f"std_{i}" for i in range(dim)]
        for q in QUANTILES:
            names += [f"p{q}_{i}" for i in range(dim)]
    else:
        raise ValueError(mode)

    feat = np.concatenate(vals).astype(np.float32)
    qc = np.array([n_patches, np.log1p(n_patches)], dtype=np.float32)
    feat = np.concatenate([feat, qc])
    names = names + ["qc_n_patches", "qc_log_n_patches"]
    return feat, names
