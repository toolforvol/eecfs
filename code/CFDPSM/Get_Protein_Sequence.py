import pandas as pd
import warnings
import os

# ============================================================================
# CONFIGURATION
# ============================================================================
input_file = ""   # Generated automatically after running AnnoSynMall.py
your_dir = ""     # Directory where resources are downloaded
task_dir = ""     # Your task directory

print(f"[INFO] Loading data...")
transcript_df = pd.read_csv(f"{your_dir}/transcript_id_info_summary.txt", sep='\t')
df = pd.read_csv(input_file, sep='\t')

print(f"[INFO] Processing data...")
os.makedirs(task_dir, exist_ok=True)


# ============================================================================
# FUNCTION
# ============================================================================
def retrieve_context_sequence_flankX(df, transcript_df, flanking, transcript_db, save_path=None):
    """
    Extract protein sequence with flanking regions, padding with X for out-of-bound indices
    
    Parameters:
    -----------
    df : DataFrame
        Contains refseq_id and aa_pos (or cds_pos)
    transcript_df : DataFrame
        Contains refseq_id and aa_sequence
    flanking : int
        Flanking length, final sequence length = 2 * flanking + 1
    transcript_db : str
        'refseq_id' (RefSeq) or 'transcript_id' (ENSEMBL)
    save_path : str, optional
        Output path for FASTA file
    
    Returns:
    --------
    df_result : DataFrame
        Processed results with flanking sequences
    warning_samples : list
        Samples with invalid positions
    """
    # Calculate aa_pos from cds_pos if not present
    if 'aa_pos' not in df.columns:
        if 'cds_pos' not in df.columns:
            raise ValueError("DataFrame missing both 'aa_pos' and 'cds_pos' columns")
        
        print("[INFO] Computing aa_pos from cds_pos")
        df = df.copy()
        df['aa_pos'] = ((df['cds_pos'] - 1) // 3) + 1
        
        invalid_mask = df['cds_pos'].isna()
        if invalid_mask.any():
            warnings.warn(f"Removing {invalid_mask.sum()} samples with NaN cds_pos")
            df = df[~invalid_mask]
    
    # Process sequences
    target_length = 2 * flanking + 1
    transcript_info = transcript_df[[transcript_db, 'aa_sequence']].drop_duplicates()
    df_merged = pd.merge(df, transcript_info, on=transcript_db, how='inner')
    
    print(f"Original samples: {len(df)}")
    print(f"Samples removed (no {transcript_db}): {len(df) - len(df_merged)}")
    
    df_merged['aa_seq_flanking'] = None
    warning_samples = []
    total_left_pad = 0
    total_right_pad = 0
    
    for index, row in df_merged.iterrows():
        sequence = row['aa_sequence']
        if pd.isna(sequence) or sequence == "":
            continue
        
        aa_pos = row['aa_pos']
        seq_len = len(sequence)
        center_index = aa_pos - 1  # Convert to 0-based
        
        # Validate position
        if center_index < 0 or center_index >= seq_len:
            warning_samples.append({
                transcript_db: row[transcript_db],
                'aa_pos': aa_pos,
                'reason': f"Position out of range (pos={aa_pos}, len={seq_len})"
            })
            continue
        
        # Calculate padding and extract
        start_index = center_index - flanking
        end_index = center_index + flanking + 1
        
        left_pad = max(0, -start_index)
        right_pad = max(0, end_index - seq_len)
        total_left_pad += left_pad
        total_right_pad += right_pad
        
        real_start = max(0, start_index)
        real_end = min(seq_len, end_index)
        extracted_seq = sequence[real_start:real_end]
        extracted_seq = 'X' * left_pad + extracted_seq + 'X' * right_pad
        
        # Validate length
        if len(extracted_seq) != target_length:
            warning_samples.append({
                transcript_db: row[transcript_db],
                'aa_pos': aa_pos,
                'reason': f"Length mismatch (got {len(extracted_seq)}, expected {target_length})"
            })
            continue
        
        df_merged.at[index, 'aa_seq_flanking'] = extracted_seq
    
    # Clean up results
    df_result = df_merged.dropna(subset=['aa_seq_flanking'])
    df_result = df_result.drop(columns=['aa_sequence'])
    
    # Statistics
    print(f"Final samples retained: {len(df_result)}")
    if len(df_result) > 0:
        print(f"Average left padding (X): {total_left_pad / len(df_result):.2f}")
        print(f"Average right padding (X): {total_right_pad / len(df_result):.2f}")
    
    if len(warning_samples) > 0:
        warnings.warn(f"Found {len(warning_samples)} samples with invalid aa_pos")
    
    # Save FASTA
    if save_path:
        with open(f"{save_path}/ref_sequences.fasta", "w") as f:
            for _, row in df_result.iterrows():
                record_id = str(row["Variant38"])
                seq = str(row["aa_seq_flanking"])
                if pd.isna(seq) or seq == "nan":
                    continue
                f.write(f">{record_id}\n{seq}\n")
    
    return df_result, warning_samples


# ============================================================================
# MAIN EXECUTION
# ============================================================================
retrieve_context_sequence_flankX(
    df=df,
    transcript_df=transcript_df,
    flanking=8,
    transcript_db='transcript_id',
    save_path=task_dir
)