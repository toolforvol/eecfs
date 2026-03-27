import subprocess
import pandas as pd
import re
import sys
import os

alg_name = "EECFS"
molecular_type = "DNA"

# 1 configuration
train_file = f"../../result{molecular_type}{alg_name}_train_2362_selectedMB.csv"
test_file = f"../../result{molecular_type}{alg_name}_test_238_selectedMB.csv"
models = [
    {
        "name": "KNN",
        "train_script": "knn_train.py",
        "model_file": f"../../result{molecular_type}{alg_name}_KNN.model"
    },
    {
        "name": "SVM",
        "train_script": "svm_train.py",
        "model_file": f"../../result{molecular_type}{alg_name}_SVM.model"
    },
    {
        "name": "LightGBM",
        "train_script": "lightgbm_train.py",
        "model_file": f"../../result{molecular_type}{alg_name}_LightGBM.model"
    }
]

# 2 parse function
def parse_cv_results(stdout):
    auc_match = re.search(
        r"Mean AUC\s*:\s*([0-9.]+)",
        stdout
    )
    aupr_match = re.search(
        r"Mean AUPR\s*:\s*([0-9.]+)",
        stdout
    )
    if auc_match and aupr_match:
        mean_auc = float(auc_match.group(1))
        mean_aupr = float(aupr_match.group(1))
    else:
        mean_auc = None
        mean_aupr = None
    return mean_auc, mean_aupr

# 3 run every clf
results = []
for model in models:
    print(f"\n========== Running {model['name']} ==========")
    command = [
        sys.executable,
        model["train_script"],
        "--train_file", train_file,
        "--model_file", model["model_file"]
    ]
    try:
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            check=True
        )
        stdout = result.stdout
        print(stdout)
        mean_auc, mean_aupr = parse_cv_results(stdout)
        results.append({
            "Model": model["name"],
            "Mean_AUC": mean_auc,
            "Mean_AUPR": mean_aupr
        })
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] {model['name']} failed")
        print(e.stderr)
        results.append({
            "Model": model["name"],
            "Mean_AUC": None,
            "Mean_AUPR": None
        })


# 4 save result
df_results = pd.DataFrame(results)
output_file = f"../../result{molecular_type}{alg_name}_CV_results.csv"
df_results.to_csv(
    output_file,
    index=False
)
print("\nSaved CV results to:")
print(output_file)
print("\nFinal Results:")
print(df_results)