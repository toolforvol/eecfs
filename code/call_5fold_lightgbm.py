"""
@description: Call 5fold_lightgbm.py (batch process for DNA, RNA, protein levels and different algos.)
"""

import subprocess
import os

# 1 5-fold best params search
train_file = "../data/MB_feature/train_2362_DNA_RNA_protein_merged_CFS_MMPC.csv"  # train file path
out_importance = "../importance/features_importance.csv"  # model importance path
model_file = "../model./train_2362_CFS_MMPC_EEDMB_EEDMB.model" # model save path

# 2 run cmd
command = [
    "python", 
    "5fold_lightgbm.py",  
    "--train_file", train_file,  
    "--out_importance", out_importance, 
    "--model_file", model_file  
]

try:
    result = subprocess.run(command, capture_output=True, text=True, check=True)
    print("Output: ")
    print(result.stdout) 
except subprocess.CalledProcessError as e:
    print("Error: ")
    print(e.stderr) 