import argparse
import pandas as pd


def retrieve_variant(df, chrom_col, pos_col, ref_col, alt_col, prefix=''):
    variant = [f'{prefix}{chro}_{pos}_{ref}_{alt}' for chro, pos, ref, alt in zip(df[chrom_col], df[pos_col], df[ref_col], df[alt_col])]
    return variant


def process(inputfile, source_file, outfile):
    df = pd.read_csv(inputfile, sep="\t")
    df['key'] = df['key'].str.split('@').str[1]
    df_source = pd.read_csv(source_file, sep="\t", header=None)
    df_source['key'] = retrieve_variant(df_source, 0, 1, 3, 4)
    df_source = df_source.merge(df, how='left', left_on='key', right_on='key', suffixes=('', '_m'))
    df_source = df_source[['key', 'MMS_alt_donor', 'MMS_pathogenicity', 'MMS_efficiency']]
    df_source.rename(columns={'key': 'Variant38'}, inplace=True)
    df_source.to_csv(outfile, sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputfile", type=str, help="输入文件")
    parser.add_argument("--source", type=str, help="测试集文件")
    parser.add_argument("--outfile", type=str, help="输出文件")
    args = parser.parse_args()
    process(args.inputfile, args.source, args.outfile)