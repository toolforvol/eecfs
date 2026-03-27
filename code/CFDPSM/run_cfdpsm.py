import subprocess
import os


# 1. Define file paths
train_file = f"../../data/MB_feature/train_2362_DNA_RNA_protein_merged_CFS_MMPC.csv"
test_file = f"../../data/MB_feature/test_238_DNA_RNA_protein_merged_CFS_MMPC.csv"
model_file = f"../../result/model/CFDPSM.model"

# =========Train=========
# 2. Construct command line call - Train
command = [
    "python",                      
    "train_lightgbm.py",                
    "--train_file", train_file,    
    "--model_file", model_file     
]

# 3. Execute command - Train
try:
    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=True
    )
    print("[Train] Script output:")
    print(result.stdout)

except subprocess.CalledProcessError as e:
    print("Script execution failed!")
    print(e.stderr)

# =========Test=========
# Construct command line call - Test
command = [
    "python",                      
    "test_lightgbm.py",                
    "--test_file", test_file,    
    "--model_path", model_file     
]

# Execute command - Test
try:
    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=True
    )
    print("[Test] Script output:")
    print(result.stdout)

except subprocess.CalledProcessError as e:
    print("Script execution failed!")
    print(e.stderr)