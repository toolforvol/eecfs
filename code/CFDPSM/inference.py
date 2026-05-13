import joblib
import pandas as pd
from sklearn.metrics import roc_auc_score, accuracy_score, precision_recall_curve, auc

def evaluate_model(model_path, test_file):
	# 1. Load test data
	test_data = pd.read_csv(test_file)
	X_test = test_data.iloc[:, 5:]
	y_test = test_data["label"]

	# 2. Load model
	model = joblib.load(model_path)

	# 3. Predict probabilities
	y_pred_prob = model.predict(X_test)
	y_pred = (y_pred_prob >= 0.5).astype(int)

	# 4. Evaluation metrics
	auc_score = roc_auc_score(y_test, y_pred_prob)
	acc_score = accuracy_score(y_test, y_pred)
	precision, recall, _ = precision_recall_curve(y_test, y_pred_prob)
	aupr_score = auc(recall, precision)

	print("===== Test Results =====")
	print(f"AUC      : {auc_score:.4f}")
	print(f"AUPR     : {aupr_score:.4f}")
	print(f"Accuracy : {acc_score:.4f}")

	return {
		"AUC": auc_score,
		"AUPR": aupr_score,
		"Accuracy": acc_score
	}

def inference_model(model_path, test_file):
	# 1. Load test data
	test_data = pd.read_csv(test_file, sep='\t')
	X_test = test_data.iloc[:, 1:]
	# 2. Load model
	model = joblib.load(model_path)
	# 3. Predict probabilities
	y_pred_prob = model.predict(X_test)
	# 4. return result
	df_result = pd.DataFrame({'Variant38': test_data['Variant38'], 'CFDPSM': y_pred_prob})
	return df_result


model_path = "../model/train_2362_CFS_MMPC_EECFS_EECFS.model"
test_file = f"{YOUR_TASK_DIR}/All.features.txt"
df_result = inference_model(model_path, test_file)
df_result.to_csv(f"{YOUR_TASK_DIR}/CFDPSM.pred.txt", sep='\t', index=False)