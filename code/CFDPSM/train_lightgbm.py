"""
@description: Perform 5fold-cross validation on the training set to obtain the optimal parameter
"""

import argparse
import os.path
import joblib
import numpy as np
import pandas as pd
import time
import re
import itertools as it
from sklearn.model_selection import KFold
from lightgbm import early_stopping
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score, precision_recall_curve, auc
import lightgbm as lgb
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
callbacks = [early_stopping(stopping_rounds=10)]
callbacks1 = [early_stopping(stopping_rounds=5)]
callbacks2 = [early_stopping(stopping_rounds=50)]


def load_data(df_train, y_train):
    train_x, train_y = df_train, y_train
    X, val_X, y, val_y = train_test_split(
        train_x,
        train_y,
        test_size=0.2,
        random_state=2024,
        stratify=train_y
    )
    lgb_train = lgb.Dataset(X, y)
    lgb_eval = lgb.Dataset(val_X, val_y, reference=lgb_train)
    return lgb_train, lgb_eval


def perform_cross_validation(X, y, params, n_folds=5):
    kf = KFold(n_splits=n_folds, shuffle=True, random_state=2024)
    auc_scores = []
    aupr_scores = []

    for train_index, test_index in kf.split(X):
        X_train, X_test = X.iloc[train_index], X.iloc[test_index]
        y_train, y_test = y.iloc[train_index], y.iloc[test_index]

        lgb_train = lgb.Dataset(X_train, y_train)
        lgb_eval = lgb.Dataset(X_test, y_test, reference=lgb_train)

        gbm = lgb.train(params,
                        lgb_train,
                        num_boost_round=1000,
                        valid_sets=lgb_eval,
                        callbacks=callbacks2)

        y_pred = gbm.predict(X_test)

        # Calculate AUC
        fold_auc = roc_auc_score(y_test, y_pred)
        auc_scores.append(fold_auc)

        # Calculate AUPR
        precision, recall, _ = precision_recall_curve(y_test, y_pred)
        fold_aupr = auc(recall, precision)
        aupr_scores.append(fold_aupr)

    return np.mean(auc_scores), np.mean(aupr_scores)


def step_training(X_train, Y_train, model_file):
    X_train = X_train.astype("float64")

    # initial params
    params = {
        'boosting_type': 'gbdt',
        'objective': 'binary',
        'metric': 'auc',
        'nthread': 4,
        'learning_rate': 0.1
    }
    max_auc = 0.0
    best_params = {}

    # ======= Param search 0: num_leaves & max_depth =======
    for num_leaves in range(5, 100, 5):
        for max_depth in range(3, 8, 1):
            params['num_leaves'] = num_leaves
            params['max_depth'] = max_depth

            cv_results = lgb.cv(
                params,
                lgb.Dataset(X_train, Y_train),  # 每次创建新Dataset
                seed=2024,
                nfold=5,
                metrics=['auc'],
                callbacks=callbacks,
                eval_train_metric=True
            )

            mean_auc = pd.Series(cv_results['valid auc-mean']).max()
            if mean_auc >= max_auc:
                max_auc = mean_auc
                best_params['num_leaves'] = num_leaves
                best_params['max_depth'] = max_depth

    # ======= Param search 1: max_bin & min_data_in_leaf =======
    for max_bin in range(5, 256, 10):
        for min_data_in_leaf in range(1, 102, 10):
            params['max_bin'] = max_bin
            params['min_data_in_leaf'] = min_data_in_leaf

            cv_results = lgb.cv(
                params,
                lgb.Dataset(X_train, Y_train),  # 重新创建Dataset
                seed=2024,
                nfold=5,
                metrics=['auc'],
                callbacks=callbacks,
                eval_train_metric=True
            )

            mean_auc = pd.Series(cv_results['valid auc-mean']).max()
            if mean_auc >= max_auc:
                max_auc = mean_auc
                best_params['max_bin'] = max_bin
                best_params['min_data_in_leaf'] = min_data_in_leaf

    # ======= Param search 2: feature_fraction & bagging =======
    for feature_fraction in [0.6, 0.7, 0.8, 0.9, 1.0]:
        for bagging_fraction in [0.6, 0.7, 0.8, 0.9, 1.0]:
            for bagging_freq in range(0, 50, 5):
                params['feature_fraction'] = feature_fraction
                params['bagging_fraction'] = bagging_fraction
                params['bagging_freq'] = bagging_freq

                cv_results = lgb.cv(
                    params,
                    lgb.Dataset(X_train, Y_train),
                    seed=2024,
                    nfold=5,
                    metrics=['auc'],
                    callbacks=callbacks,
                    eval_train_metric=True
                )

                mean_auc = pd.Series(cv_results['valid auc-mean']).max()
                if mean_auc >= max_auc:
                    max_auc = mean_auc
                    best_params['feature_fraction'] = feature_fraction
                    best_params['bagging_fraction'] = bagging_fraction
                    best_params['bagging_freq'] = bagging_freq

    # ======= Param search 3: lambda_l1 & lambda_l2 =======
    for lambda_l1 in [1e-5, 1e-3, 1e-1, 0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0]:
        for lambda_l2 in [1e-5, 1e-3, 1e-1, 0.0, 0.1, 0.4, 0.6, 0.7, 0.9, 1.0]:
            params['lambda_l1'] = lambda_l1
            params['lambda_l2'] = lambda_l2

            cv_results = lgb.cv(
                params,
                lgb.Dataset(X_train, Y_train),
                seed=2024,
                nfold=5,
                metrics=['auc'],
                callbacks=callbacks,
                eval_train_metric=True
            )

            mean_auc = pd.Series(cv_results['valid auc-mean']).max()
            if mean_auc >= max_auc:
                max_auc = mean_auc
                best_params['lambda_l1'] = lambda_l1
                best_params['lambda_l2'] = lambda_l2

    # ======= Param search 4: min_split_gain =======
    for min_split_gain in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
        params['min_split_gain'] = min_split_gain

        cv_results = lgb.cv(
            params,
            lgb.Dataset(X_train, Y_train),
            seed=2024,
            nfold=5,
            metrics=['auc'],
            callbacks=callbacks,
            eval_train_metric=True
        )

        mean_auc = pd.Series(cv_results['valid auc-mean']).max()
        if mean_auc >= max_auc:
            max_auc = mean_auc
            best_params['min_split_gain'] = min_split_gain

    # ======= Final params =======
    final_params = {
        'boosting_type': 'gbdt',
        'objective': 'binary',
        'metric': {'binary_logloss', 'auc'},
        'nthread': 4,
        'learning_rate': 0.1,
        'verbose': 5,
        'is_unbalance': False
    }
    final_params.update(best_params)

    # set defaults if missing
    defaults = {
        'num_leaves': 31,
        'max_depth': -1,
        'max_bin': 255,
        'min_data_in_leaf': 20,
        'feature_fraction': 1.0,
        'bagging_fraction': 1.0,
        'bagging_freq': 0,
        'lambda_l1': 0.0,
        'lambda_l2': 0.0,
        'min_split_gain': 0.0
    }
    for k, v in defaults.items():
        final_params.setdefault(k, v)

    # ======= Train final model =======
    lgb_train, lgb_eval = load_data(X_train, Y_train)
    
    # 创建目录
    os.makedirs(os.path.dirname(model_file), exist_ok=True)
    
    gbm = lgb.train(final_params,
                    lgb_train,
                    num_boost_round=1000,
                    valid_sets=lgb_eval,
                    callbacks=callbacks2)
    
    # 保存模型
    joblib.dump(gbm, model_file)

    # 5-fold CV evaluation
    mean_auc, mean_aupr = perform_cross_validation(X_train, Y_train, final_params)
    print(f"\nCross-validation results with best parameters:")
    print(f"Mean AUC: {mean_auc:.4f}")
    print(f"Mean AUPR: {mean_aupr:.4f}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--train_file", type=str)
    parser.add_argument("--model_file", type=str)
    args = parser.parse_args()

    df_train = pd.read_csv(args.train_file) 
    train_data = df_train.iloc[:, 5:] # feature columns
    Y_train = df_train["label"] # label column
    step_training(train_data, Y_train, args.model_file)