"""Historical microbiome transform; split and taxa are fixed from Discovery."""

import numpy as np
import pandas as pd


def fit_zscore(df):
    mean_ = df.mean(axis=0)
    std_ = df.std(axis=0, ddof=0)
    std_[std_ == 0] = 1.0
    std_ = std_.fillna(1.0)
    return mean_, std_


def apply_zscore(df, mean_, std_):
    return (df - mean_) / std_


def clr_transform(df, pseudocount=1.0):
    arr = df.astype(float).values + pseudocount
    log_arr = np.log(arr)
    gm = log_arr.mean(axis=1, keepdims=True)
    clr = log_arr - gm
    return pd.DataFrame(clr, index=df.index, columns=df.columns)


def prepare_microbiome(train, validation):
    zero_ratio = (train == 0).mean(axis=0)
    taxa = zero_ratio[zero_ratio <= 0.80].index.tolist()
    train_clr = clr_transform(train[taxa], pseudocount=1.0)
    validation_clr = clr_transform(validation[taxa], pseudocount=1.0)
    mean, std = fit_zscore(train_clr)
    return (
        apply_zscore(train_clr, mean, std),
        apply_zscore(validation_clr, mean, std),
        taxa,
    )
