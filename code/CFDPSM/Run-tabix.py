"""
@description: synmall local annotation
@author: yechen
@date: 2025-04-28
"""

import argparse
import os
import subprocess
import pandas as pd
from datetime import datetime
import uuid
import warnings
warnings.filterwarnings('ignore', category=pd.errors.DtypeWarning)


anno_root_path = "" # where you download the annotation, download from https://bioinfo.ahu.edu.cn/synMall/#/download
anno_root_header = "./config/Annotation.header"


def retrieve_variant(df, chrom_col, pos_col, ref_col, alt_col, prefix=''):
    variant = [f'{prefix}{chro}_{pos}_{ref}/{alt}' for chro, pos, ref, alt in zip(df[chrom_col], df[pos_col], df[ref_col], df[alt_col])]
    return variant


def merge_output(original_vid_list, output_dir, sto_path):
    files = [os.path.join(output_dir, f) for f in os.listdir(output_dir) if f.endswith('_out.txt')]
    df_list = [pd.read_csv(f, sep='\t', header=None) for f in files]
    df_concat = pd.concat(df_list, ignore_index=True)

    # dedup
    df_concat['Variant38'] = retrieve_variant(df_concat, 0, 1, 2, 3)
    df_concat = df_concat[df_concat['Variant38'].isin(original_vid_list)]

    # write final result
    with open(sto_path, 'w') as f_out:
        with open(anno_root_header, 'r') as header_in:
            f_out.write(header_in.read())  # write header
        df_concat.drop(columns=['Variant38']).to_csv(f_out, sep='\t', index=False, header=False)


def annotate_vcf(input_vcf, output_file):
    try:
        folder_name = str(uuid.uuid4())
        os.makedirs(folder_name, exist_ok=True)

        input_vcf['Variant38'] = retrieve_variant(input_vcf, '#CHROM', 'POS', 'REF', 'ALT')
        input_vcf = input_vcf.drop_duplicates(subset='Variant38')
        original_vid_list = input_vcf['Variant38'].tolist()

        dfs = input_vcf.groupby('#CHROM')

        print(f"[INFO] Annotating...")
        temp_files = []
        for chr_num, tdf in dfs:
            chr_file = f"{folder_name}/chr{chr_num}.txt"
            tdf.to_csv(chr_file, sep='\t', index=False)
            temp_files.append(chr_file)

            out_file = f"{folder_name}/chr{chr_num}_out.txt"
            cmd = ["tabix", "-R", chr_file, f"{anno_root_path}chr{chr_num}.tsv.gz"]
            with open(out_file, "w") as fout:
                subprocess.run(cmd, stdout=fout, check=True)
            temp_files.append(out_file)

        print(f"[INFO] Merging output...")
        merge_output(original_vid_list, folder_name, output_file)

        # clean cache
        for temp_file in temp_files:
            os.remove(temp_file)
        os.rmdir(folder_name)
        print(f"[INFO] Finished successfully.")

    except Exception as e:
        print(f"[ERROR] {str(e)}")
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Annotate VCF with bgzip indexed files per chromosome.")
    parser.add_argument("-i", "--input", required=True, help="Input VCF file path")
    parser.add_argument("-o", "--output", required=True, help="Output annotated VCF file path")
    args = parser.parse_args()

    try:
        input_vcf = pd.read_csv(args.input, sep='\t')
        if '#CHROM' not in input_vcf.columns:
            raise ValueError("Input file missing #CHROM header.")
    except Exception as e:
        print(f"[ERROR] Failed to read input VCF: {str(e)}")
        exit(1)

    annotate_vcf(input_vcf, args.output)
