#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np

import xgboost as xgb
from sklearn.metrics import mean_squared_error, r2_score
from scipy.stats import pearsonr

def main():
    parser = argparse.ArgumentParser(
        description="Train & evaluate a GBM using an external test set CSV."
    )
    parser.add_argument(
        "--features_csv", "-f", required=True,
        help="CSV with columns: PDBID, pK, plus feature columns (training+full set)."
    )
    parser.add_argument(
        "--test_csv", "-e", required=True,
        help="CSV with the same schema whose rows define the test set."
    )
    parser.add_argument(
        "--random_state", "-s", type=int, default=42,
        help="Random seed for reproducibility."
    )
    args = parser.parse_args()

    # --- load full & test
    df_full = pd.read_csv(args.features_csv).dropna(subset=["pK"])
    df_test = pd.read_csv(args.test_csv).dropna(subset=["pK"])
    for col in ("PDBID", "pK"):
        if col not in df_full.columns or col not in df_test.columns:
            raise ValueError(f"Both CSVs must include '{col}' column")

    # identify feature columns
    feature_cols = [c for c in df_full.columns if c not in ("PDBID", "pK")]

    # prepare test set
    df_test = df_test.set_index("PDBID")
    X_test = df_test[feature_cols]
    y_test = df_test["pK"].values

    # prepare training set = rows in full not in test
    train_mask = ~df_full["PDBID"].isin(df_test.index)
    df_train = df_full[train_mask].set_index("PDBID")
    X_train = df_train[feature_cols]
    y_train = df_train["pK"].values

    print(f"Training on {len(X_train)} samples; testing on {len(X_test)} samples")

    # --- train GBM
    model = xgb.XGBRegressor(
        n_estimators=10000,
        max_depth=7,
        learning_rate=0.01,
        eval_metric='rmse',
        subsample = .7,
        random_state=args.random_state,
        colsample_bytree=.6,
        n_jobs=-1,               # all your cores
        tree_method="approx",
        early_stopping_rounds=50
    )

    model.fit(
        X_train, y_train,
        eval_set=[(X_test, y_test)],
        verbose=True
    )

    # --- evaluate
    print("\n=== Test Set Evaluation ===")
    y_pred = model.predict(X_test)
    print(f" R²        = {r2_score(y_test, y_pred):.4f}")
    print(f" RMSE      = {mean_squared_error(y_test, y_pred):.4f}")
    r, p = pearsonr(y_test, y_pred)
    print(f" Pearson r = {r:.4f} (p={p:.2e})")
    

if __name__ == "__main__":
    main()
