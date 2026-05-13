import pandas as pd
import os
from Bio import SeqIO
import subprocess

# ============================================================================
# CONFIGURATION
# ============================================================================
task_dir = ""      # Your task directory
genome = "hg38"
input_path = ""    # Your input VCF file
context = 55
parent_path = './bedtools/'  # Bedtools directory

print(f'[INFO] New task starts at {task_dir}')
os.makedirs(task_dir, exist_ok=True)


# ============================================================================
# FUNCTIONS
# ============================================================================
def vcf_to_bed(input_path, context, task_dir):
    """
    Convert VCF to BED format with flanking regions
    """
    df = pd.read_csv(input_path, sep='\t')
    bed = pd.DataFrame()
    bed['#CHR'] = "chr" + df['#CHROM'].astype(str)
    bed['START'] = df['POS'] - context - 1
    bed['END'] = df['POS'] + context
    
    df['region'] = bed['#CHR'] + ':' + bed['START'].astype(str) + '-' + bed['END'].astype(str)
    
    bed.to_csv(f"{task_dir}/input.bed", sep='\t', index=False)
    df.to_csv(f"{input_path}.withRegion", sep='\t', index=False)


def bed_to_fasta(bed_file, task_dir):
    """
    Extract DNA sequences from BED regions using bedtools
    """
    bed_tool = "/data5/yechen1/9_Tools/3-bedtools/bedtools"
    reference_genome = f"/data5/yechen1/9_Tools/3-bedtools/reference_genome/{genome}.fa"
    fasta_out = f"{task_dir}/bedtool_output.fasta"
    
    cmd = [bed_tool, "getfasta", "-fi", reference_genome, "-bed", bed_file, "-fo", fasta_out]
    print(" ".join(cmd))
    subprocess.run(cmd, check=True, text=True, capture_output=True)


def retrieve_sequence(input_path, task_dir, context):
    """
    Match sequences to DataFrame and generate ref/alt sequences
    """
    # Load VCF data
    df = pd.read_csv(f"{input_path}.withRegion", sep='\t')
    
    # Load FASTA sequences
    fasta_dict = {}
    for record in SeqIO.parse(f"{task_dir}/bedtool_output.fasta", "fasta"):
        fasta_dict[record.id] = str(record.seq).upper()
    
    # Map sequences to regions
    df['ref_fasta'] = df['region'].map(fasta_dict)
    
    # Generate ALT sequences by replacing center base
    center_idx = context
    alt_fasta_list = []
    for _, row in df.iterrows():
        ref_seq = row['ref_fasta']
        if pd.isna(ref_seq):
            alt_fasta_list.append(None)
            continue
        
        ref_list = list(ref_seq)
        ref_list[center_idx] = row['ALT']
        alt_fasta_list.append("".join(ref_list))
    
    df['alt_fasta'] = alt_fasta_list
    
    # Save FASTA files
    with open(f"{task_dir}/ref_sequences.fasta", "w") as ref_out, \
         open(f"{task_dir}/alt_sequences.fasta", "w") as alt_out:
        for _, row in df.iterrows():
            if pd.isna(row['ref_fasta']):
                continue
            ref_out.write(f">{row['Variant38']}\n{row['ref_fasta']}\n")
            alt_out.write(f">{row['Variant38']}\n{row['alt_fasta']}\n")
    
    print("[INFO] FASTA files saved")


# ============================================================================
# MAIN EXECUTION
# ============================================================================
vcf_to_bed(input_path, context, task_dir)
bed_to_fasta(f"{task_dir}/input.bed", task_dir)
retrieve_sequence(input_path, task_dir, context)