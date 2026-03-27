import argparse
import joblib
import pandas as pd
import numpy as np

from sklearn.neighbors import KNeighborsClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.metrics import roc_auc_score, precision_recall_curve, auc

import warnings
warnings.filterwarnings("ignore")


def train_knn(X, y, model_file):
    pipe = Pipeline([
        ("scaler", StandardScaler()),
        ("knn", KNeighborsClassifier())
    ])

    param_grid = {
        "knn__n_neighbors": [3, 5, 7, 9, 11, 15],
        "knn__weights": ["uniform", "distance"],
        "knn__metric": ["euclidean", "manhattan"]
    }

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=2024)

    grid = GridSearchCV(
        pipe,
        param_grid=param_grid,
        scoring="roc_auc",
        cv=cv,
        n_jobs=-1,
        verbose=1
    )

    grid.fit(X, y)

    best_model = grid.best_estimator_
    print("\n********** Best Parameters (KNN) *********")
    print(grid.best_params_)
    print("Best CV AUC:", grid.best_score_)
    print("*****************************************\n")

    joblib.dump(best_model, model_file)
    print(f"[SUCCESS] KNN model saved to: {model_file}")

    return best_model


def cross_validation_eval(model, X, y):
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=2024)
    auc_list, aupr_list = [], []

    for train_idx, test_idx in cv.split(X, y):
        X_tr, X_te = X.iloc[train_idx], X.iloc[test_idx]
        y_tr, y_te = y.iloc[train_idx], y.iloc[test_idx]

        model.fit(X_tr, y_tr)
        y_prob = model.predict_proba(X_te)[:, 1]

        auc_list.append(roc_auc_score(y_te, y_prob))
        precision, recall, _ = precision_recall_curve(y_te, y_prob)
        aupr_list.append(auc(recall, precision))

    print("5-fold CV Results:")
    print(f"Mean AUC : {np.mean(auc_list):.4f}")
    print(f"Mean AUPR: {np.mean(aupr_list):.4f}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--train_file", type=str, required=True)
    parser.add_argument("--model_file", type=str, required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.train_file)
    X = df.iloc[:, 5:]
    y = df["label"]

    model = train_knn(X, y, args.model_file)
    cross_validation_eval(model, X, y)
