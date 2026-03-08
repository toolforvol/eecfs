"""
@description: evaluate the final model.
"""


import joblib
import pandas as pd
from sklearn.metrics import roc_auc_score, accuracy_score, precision_recall_curve, auc


def evaluate_model(model_path, test_file, prediction_output_file=None):
    """
    Desc:
        load model, predict, and evaluate.
    Args:
        - model_path: model path.
        - test_file: test data file path.
        - error_output_file: wrong prediction of results
    Returns:
        - predictions: the predicted logits.
        - metrics: metrics include auc, aupr, etc.
    """
    # load test data
    test_data = pd.read_csv(test_file)
    X_test = test_data.iloc[:, 5:]  # feature columns
    y_test = test_data["label"]     # label column

    # load model
    model = joblib.load(model_path)

    # predict probes
    y_pred_prob = model.predict(X_test)

    # the predicted label
    y_pred = (y_pred_prob >= 0.5).astype(int)

    # evaluation metrics
    auc_score = roc_auc_score(y_test, y_pred_prob)
    accuracy = accuracy_score(y_test, y_pred)
    precision, recall, _ = precision_recall_curve(y_test, y_pred_prob)
    aupr = auc(recall, precision)

    metrics = {
        "AUC": auc_score,
        "AUPR": aupr,
        "Accuracy": accuracy
    }

    print("Evaluation Metrics:")
    print(f"AUC: {auc_score:.4f}")
    print(f"AUPR: {aupr:.4f}")
    print(f"Accuracy: {accuracy:.4f}")

    # save the results
    if prediction_output_file:
        full_result_df = test_data.copy()
        full_result_df["y_pred_prob"] = y_pred_prob
        full_result_df.to_csv(prediction_output_file, index=False)
        print(f"Save result to: {prediction_output_file}")

    return y_pred_prob, metrics

model_path = "../model/train_2362_CFS_MMPC_EECFS_EECFS.model"
test_file = "../data/MB_feature/"
predictions, metrics = evaluate_model(model_path, test_file)