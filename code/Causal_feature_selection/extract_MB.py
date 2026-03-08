"""
@description: This script help to extract the MB feature data.
"""


import pandas as pd

# input file path
input_file = '../../data/Protein/test_238_protein_17aa_final_feature_with_label_renamed.csv'

# output file path
output_file = '../../data/Protein/test_238_protein_17aa_EECFS_features_2362.csv'

# MATLAB indice
# matlab_indices = [11, 40, 83, 84, 2111, 2524, 5371, 5373, 5388,
#                  5393, 6332, 6519, 9283, 10580, 10582, 10959, 11538]   # EECFS DNA

# matlab_indices = [11, 40, 83, 84, 2100, 2111, 2514, 2524, 4220, 5371, 5373, 5388,
#                   6332, 6519, 9498, 10580, 10582, 11538]   # EECFS DNA 2362

# matlab_indices = [18, 21, 22, 3512, 3550, 4400, 5133, 5220, 6257]  # EECFS RNA
# matlab_indices = [2, 3, 4, 5, 18, 21, 22, 1299, 1564, 5710]  # EECFS new RNA
# matlab_indices = [2, 3, 4, 5, 18, 21, 22, 1299, 5710, 5968]  # EECFS new RNA 2362
# matlab_indices = [3789, 3811, 4032, 4055, 4058, 4079, 4242]  # EECFS PROTEIN 2362 


# matlab_indices = [11, 40, 84]  #EECFS DNA no seq
# matlab_indices = [1519, 2007, 2015, 2424, 5029, 5275, 5278, 5292, 5297, 5332, 6349, 6423, 6434,
#                   8688, 8780, 9000, 9193, 9297, 9405, 9508, 9832, 10866, 11387, 11442]  #EECFS DNA seq

# matlab_indices = [18, 21, 22]  # RNA no seq
# matlab_indices = [328, 332, 1245, 1248, 1256, 1262, 1265, 1273, 2074, 4274, 4278, 4522, 5937, 5973,
#                   5977, 6127, 6194, 6204, 6236, 6463, 6466, 6933, 7644]  # EECFS RNA seq

# matlab_indices = [26, 27, 28]  #EECFS PROTEIN no seq
# matlab_indices = [3789, 3811, 4032, 4055, 4058, 4079, 4242]  #EECFS PROTEIN seq

# matlab_indices = [8, 11, 84, 2103, 2111, 2524, 4218, 5369, 5371, 5393, 6344, 6453,
#                   6522, 6567, 7636, 7655, 8318, 8748, 8882, 9283, 10580, 10582, 10583, 11537, 11538]  # MMMB DNA

# matlab_indices = [8, 11, 84, 2100, 2111, 2524, 4218, 5369, 5371, 5393, 6344, 6453, 6522,
#                   6567, 7636, 7655, 8318, 8748, 8882, 9498, 9500, 10580, 10582, 10583, 11537, 11538]  # MMMB 2362 DNA

# matlab_indices = [10, 14, 17, 18, 21, 22, 1284, 1294, 1559, 4122, 4293, 4331,
#                   4378, 5129, 5133, 5148, 5157, 5220, 5222, 5974, 6257]  # MMMB RNA

# matlab_indices = [2, 3, 4, 5, 10, 14, 18, 21, 22, 26, 1289, 1299, 1564, 1706, 4127, 4323,
#                  5137, 6338, 7670]  # MMMB new RNA

# matlab_indices = [2, 3, 4, 5, 10, 14, 18, 21, 22, 26, 814, 1289, 1299, 1564, 1706, 4127, 4323,
#                  5137, 5968, 6338, 7670]  # MMMB new RNA 2362

# matlab_indices = [87, 107, 127, 187, 200, 207, 233, 3107, 3780, 3877, 4028, 4032, 4038, 4039, 4055, 4079, 4249] # MMMB PROTEIN 2362 

# matlab_indices = [11, 84, 2103, 2111, 2520, 2524, 5369, 5371, 5388, 6260, 6523,
#                   7602, 8934, 8942, 9074, 9106, 9283, 9393, 9394, 9500, 10253, 11538]  # HITONMB DNA 2362 

# matlab_indices = [10, 14, 17, 18, 21, 22, 1294, 1559, 3512, 4122, 4177, 4225, 4331, 4332,
#                   5129, 5140, 5141, 5220, 5222, 5974, 6123, 6230, 6257, 6484]  # HITONMB RNA

# matlab_indices = [2, 3, 4, 5, 10, 14, 18, 21, 22, 26, 33, 1002, 1299, 1564, 1616, 3448, 4303, 6338, 6489]  # HITONMB new RNA

# matlab_indices = [2, 3, 4, 5, 10, 14, 18, 21, 22, 26, 33, 1299, 1564, 3448, 4303, 5968, 6338, 6489]  # HITONMB new RNA 2362

# matlab_indices = [187, 193, 207, 213, 3789, 4026, 4029, 4032, 4055, 4079, 4242, 4243, 4245, 4247, 4358]  # HITONMB PROTEIN

# matlab_indices = [187, 193, 207, 213, 3789, 4026, 4029, 4032, 4055, 4079, 4242, 4243, 4244, 4245, 4247, 4358]  # HITONMB PROTEIN 2362

# matlab_indices = [8, 11, 84]  # PCMB DNA 2362 
# matlab_indices = [10, 14, 21, 22]  # PCMB RNA
# matlab_indices = [2, 3, 4, 5, 10, 18, 21, 26]  # PCMB new RNA
# matlab_indices = [2, 3, 4, 5, 10, 21, 26]  # PCMB new RNA 2362
# matlab_indices = [4032]  # PCMB PROTEIN 2362 

#matlab_indices = [11, 40, 83, 84, 1615, 1647, 2111, 2198, 2520, 2524, 5371, 6395, 7541
#                 , 7619, 7645, 8942, 9074, 9106, 9174, 9283, 9393, 9498, 11538]   # BAMB DNA

# matlab_indices = [11, 40, 83, 84, 1615, 2111, 2127, 2198, 2520, 2524, 4701, 5371, 6395, 7541
#                 , 8942, 9074, 9106, 9174, 9283, 9289, 9393, 9498, 11538]   # BAMB DNA 2362

# matlab_indices = [10, 17, 21, 1294, 1559, 3512, 3525, 4249, 4281, 4543, 5142, 5220, 5705, 5963, 6135, 6484]  # BAMB RNA
# matlab_indices = [2, 3, 4, 5, 10, 21, 1289, 1299, 1564, 4310, 5132, 5710, 6489]  # BAMB new RNA
# matlab_indices = [2, 3, 4, 5, 10, 21, 1299, 1564, 5225, 5710, 6489]  # BAMB new RNA 2362
# matlab_indices = [190, 207, 2958, 3780, 3789, 3796, 4026, 4029, 4032, 4055, 4079, 4242, 4247]  # BAMB PROTEIN 2362 

# matlab_indices = [11, 83, 84, 1615, 2112, 2520, 5371, 5388, 6338, 6519, 7871, 8942, 9106, 9283, 9289, 11538]   # EEMB DNA
# matlab_indices = [11, 83, 84, 1615, 2112, 2520, 2524, 5371, 5388, 6338, 6519, 7871, 8942, 9106, 9283, 9289, 11538]   # EEMB DNA 2362

# matlab_indices = [18, 21, 22, 1294, 1559, 3525, 4262, 4332, 5220, 6135, 6257]  # EEMB RNA
# matlab_indices = [2, 3, 4, 5, 10, 17, 21, 1299, 1564, 5710, 6489]  # EEMB new RNA
# matlab_indices = [2, 3, 4, 5, 10, 17, 21, 1299, 1564, 5968, 6489]  # EEMB new RNA 2362
# matlab_indices = [207, 213, 3789, 3794, 4026, 4029, 4032, 4055, 4058, 4079, 4242]  # EEMB PROTEIN 2362 

# matlab_indices = [8, 11, 84, 2111, 5371, 11538]   # EDMB DNA
# matlab_indices = [8, 11, 84, 2100, 4218, 5371, 10582, 11538]   # EDMB DNA 2362
# matlab_indices = [10, 14, 18, 21, 22, 6257]  # EDMB RNA
# matlab_indices = [2, 3, 4, 5, 10, 18, 21, 22, 25, 26, 1299, 1564]  # EDMB new RNA
# matlab_indices = [2, 3, 4, 5, 10, 21, 22, 25, 26, 1299, 5968]  # EDMB new RNA 2362
# matlab_indices = [4032, 4055, 4079]  # EDMB PROTEIN 2362 

# matlab_indices = [11, 40, 83, 84, 1615, 2524, 5371, 5373, 5388, 6519, 9042, 9283, 9497, 10580, 10582, 11538]   # CFS_MMPC DNA
# matlab_indices = [11, 40, 83, 84, 2111, 2524, 5371, 5373, 5388, 6519, 10580, 10582, 11538]   # CFS_MMPC DNA 2362
# matlab_indices = [18, 21, 22, 1284, 1294, 1559, 4543, 5133, 5220, 6257]  # CFS_MMPC RNA
# matlab_indices = [2, 3, 4, 5, 10, 18, 21, 1289, 1299, 1564, 5710, 7670]  # CFS_MMPC new RNA
# matlab_indices = [2, 3, 4, 5, 10, 18, 21, 1289, 1299, 1564, 5710, 5968, 7670]  # CFS_MMPC new RNA 2362
# matlab_indices = [187, 207, 233, 3789, 4032, 4055, 4058, 4079]  # CFS_MMPC PROTEIN 2362 

# label column index
label_col_index = 4
real_indices = [i + label_col_index - 1 for i in matlab_indices]  # NOTE: python 0-based index

# read data
df = pd.read_csv(input_file)

print("The selected columns are:")
for idx in real_indices:
    print(f"{idx}: {df.columns[idx]}")

# extract the columns with 5 meta-infos
info_cols = df.iloc[:, :5]
selected_features = df.iloc[:, real_indices]

# concat and save
final_df = pd.concat([info_cols, selected_features], axis=1)
final_df.to_csv(output_file, index=False)