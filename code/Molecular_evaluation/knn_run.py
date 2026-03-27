import subprocess
import os


alg_name = "EECFS"

# 1. Define file paths
train_file = f"../../result/DNA/{alg_name}_train_2362_selectedMB.csv"
test_file = f"../../result/DNA/{alg_name}_test_238_selectedMB.csv"
model_file = f"../../result/DNA/{alg_name}_KNN.model"

# =========Train=========
# 2. Construct command line call - Train
command = [
    "python",                      
    "knn_train.py",                
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
        errors="replace"  
    )
    print("[Train] KNN script output:")
    print(result.stdout)

except subprocess.CalledProcessError as e:
    print("KNN script execution failed!")
    print(e.stderr)

# =========Test=========
# Construct command line call - Test
command = [
    "python",                      
    "common_test.py",                
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
        errors="replace"  
    )
    print("[Test] KNN script output:")
    print(result.stdout)

except subprocess.CalledProcessError as e:
    print("KNN script execution failed!")
    print(e.stderr)