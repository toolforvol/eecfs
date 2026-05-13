import pandas as pd
import argparse
from Bio import SeqIO


def get_vid_from_fasta(fasta_file):
    vid_list = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        vid = str(record.id).split('|')[0]
        vid_list.append(vid)
    return vid_list


def select_features(fasta_file, ref_fea, alt_fea, output_file):
	vid_list = get_vid_from_fasta(fasta_file)
	ref_fea_df = pd.read_csv(ref_fea, sep=' ', usecols=['CCUU', 'AAUA'])
	ref_fea_df.columns = [item+'_ref' for item in ref_fea_df.columns]
	alt_fea_df = pd.read_csv(alt_fea, sep=' ', usecols=['AAUA'])
	alt_fea_df.columns = [item+'_alt' for item in alt_fea_df.columns]

	result_df = pd.concat([ref_fea_df, alt_fea_df], axis=1)
	result_df['AAUA_diff'] = result_df['AAUA_alt'] - result_df['AAUA_ref']
	result_df.insert(0, 'Variant38', vid_list)
	selected_features = ['Variant38', 'CCUU_ref', 'AAUA_diff']
	result_df = result_df[selected_features]
	result_df.to_csv(output_file, index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Merge feautrues for ftrCOOL")
    parser.add_argument("--fasta_file", type=str, required=True, help="Reference fasta file")
    parser.add_argument("--ref_fea", type=str, help="reference feauture file")
    parser.add_argument("--alt_fea", type=str, help="alternate feauture file")
    parser.add_argument("--output_file", type=str, required=True, help="Output feature csv path")

    args = parser.parse_args()
    fasta_file = args.fasta_file
    ref_fea = args.ref_fea
    alt_fea = args.alt_fea
    output_file = args.output_file
    select_features(fasta_file, ref_fea, alt_fea, output_file)


"""
python Filter_features.py \
--fasta_file /data4/jinfangfang/Project/ftrCOOL/EECFS/RNA/ref_sequences_rna.fasta \
--ref_fea /data4/jinfangfang/Project/ftrCOOL/EECFS/RNA/ExpectedValKmerNUC_RNA_ref.xlsx \
--alt_fea /data4/jinfangfang/Project/ftrCOOL/EECFS/RNA/ExpectedValKmerNUC_RNA_alt.xlsx \
--output_file /data4/jinfangfang/Project/ftrCOOL/EECFS/RNA/ftrCOOL.csv
"""