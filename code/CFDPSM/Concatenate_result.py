import pandas as pd

root_path = YOUR_TASK_DIR # 🔻 Modify this
protein_fea = pd.read_csv(f"{root_path}/protein_seq/ifo.fea.protein.csv")
print(f"[INFO] protein ifo: {len(protein_fea)}")

dna_bio_fea = pd.read_csv(f"{root_path}/SynMall.output.txt", sep='\t')
result_df = protein_fea.merge(dna_bio_fea, how='left', left_on='Variant38', right_on='Variant38', suffixes=('', '_m'))
print(f"[INFO] merge dna bio: {len(result_df)}")

dna_ifo_fea = pd.read_csv(f"{root_path}/DNA_seq/ifo.fea.DNA.csv")
result_df = result_df.merge(dna_ifo_fea, how='left', left_on='Variant38', right_on='Variant38', suffixes=('', '_m'))
print(f"[INFO] merge dna ifo: {len(result_df)}")

dna_mf_fea = pd.read_csv(f"{root_path}/DNA_seq/mathfeature.csv")
result_df = result_df.merge(dna_mf_fea, how='left', left_on='Variant38', right_on='Variant38', suffixes=('', '_m'))
print(f"[INFO] merge dna mf: {len(result_df)}")

rna_mms_fea = pd.read_csv(f"{root_path}/MMSplice/MMSplice.result.txt", sep='\t')
rna_mms_fea['Variant38'] = rna_mms_fea['Variant38'].apply(lambda x: '/'.join(x.rsplit('_', 1)))
result_df = result_df.merge(rna_mms_fea, how='left', left_on='Variant38', right_on='Variant38', suffixes=('', '_m'))
print(f"[INFO] merge rna mms: {len(result_df)}")

rna_ifo_fea = pd.read_csv(f"{root_path}/RNA_seq/ifo.fea.RNA.csv")
result_df = result_df.merge(rna_ifo_fea, how='left', left_on='Variant38', right_on='Variant38', suffixes=('', '_m'))
print(f"[INFO] merge rna ifo: {len(result_df)}")

rna_ftr_fea = pd.read_csv(f"{root_path}/RNA_seq/ftrCOOL.csv")
result_df = result_df.merge(rna_ftr_fea, how='left', left_on='Variant38', right_on='Variant38', suffixes=('', '_m'))
print(f"[INFO] merge rna ftr: {len(result_df)}")

header_info = pd.read_csv("./config/combined_header.txt", header=None)[0].tolist()
result_df = result_df[header_info]
result_df.to_csv(f"{root_path}/All.features.txt", sep='\t', index=False)