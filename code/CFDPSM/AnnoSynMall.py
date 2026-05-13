# ============================================================================
# CONFIGURATION
# ============================================================================
YOUR_ROOT_PATH = ""  # where you downloaded the resources
TASK_DIR = ""        # where you run the task
INPUT_FILE = ""      # your input vcf file

# ============================================================================
# IMPORTS
# ============================================================================
import pandas as pd
import subprocess
import os
import warnings
warnings.filterwarnings("ignore")

# ============================================================================
# CONSTANTS
# ============================================================================
SYNMALL_TABIX_SCRIPT = "./Run-tabix.py"

# Transcript meta info columns
TRANSCRIPT_META_LIST = [
    'transcript_id', 'cds_pos', 'chr', 'pos', 'ref', 'alt', 
    'mutation_source', 'codon_alternate', 'delta_ess', 'delta_ese', 
    'delta_esr', 'splice_site_acceptor', 'splice_site_donor'
]

# ============================================================================
# INITIALIZATION
# ============================================================================
print("[INFO] Loading basic information...")
original_header = pd.read_csv(
    "./config/synmall_original_header.info", sep='\t', header=None
)[0].tolist()

selected_header = pd.read_csv(
    "./config/synmall_reformat_header.info", sep='\t', header=None
)[0].tolist()

canonical_transcript_set = set(
    pd.read_csv(
        f"{YOUR_ROOT_PATH}/Canonical_transcript.txt", sep='\t', 
        usecols=['Transcript_stable_ID']
    )['Transcript_stable_ID'].tolist()
)

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================
def retrieve_raw_synmall(input_file, output_file):
    """
    Retrieve raw tabix synmall file
    
    Parameters:
    -----------
    input_file : str
        Path to input VCF file
    output_file : str
        Path where raw output will be saved
    """
    cmd = f"python {SYNMALL_TABIX_SCRIPT} -i {input_file} -o {output_file}"
    print(f'[INFO] Exec synmall cmd: {cmd}')
    subprocess.run(cmd, shell=True)


def retrieve_variant(df, chrom_col, pos_col, ref_col, alt_col, prefix=''):
    """
    Generate a variant ID column format as chrom_pos_ref/alt
    
    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame containing variant data
    chrom_col, pos_col, ref_col, alt_col : str
        Column names for chromosome, position, reference, and alternate alleles
    prefix : str
        Optional prefix for variant ID
    
    Returns:
    --------
    list
        List of variant ID strings
    """
    variant = [
        f'{prefix}{chro}_{pos}_{ref}/{alt}' 
        for chro, pos, ref, alt in zip(
            df[chrom_col], df[pos_col], df[ref_col], df[alt_col]
        )
    ]
    return variant


def select_transcript_info(info_str, canonical_set, meta_list):
    """
    Select the best transcript record from annotation string
    
    Parameters:
    -----------
    info_str : str
        Raw transcript info string
    canonical_set : set
        Set of canonical transcript IDs
    meta_list : list
        List of metadata column names
    
    Returns:
    --------
    dict
        Dictionary with selected transcript information
    """
    import numpy as np
    
    if pd.isna(info_str):
        return {col: np.nan for col in meta_list}
    
    # Split all records
    records = [
        dict(zip(meta_list, record.split('|'))) 
        for record in info_str.split(';')
    ]
    
    # Filter canonical transcripts
    canonical_records = [
        r for r in records if r['transcript_id'] in canonical_set
    ]
    
    selected_record = None
    
    if len(canonical_records) == 1:
        selected_record = canonical_records[0]
    elif len(canonical_records) > 1:
        # Multiple records - select longest mutation_source (most confident)
        selected_record = max(
            canonical_records, 
            key=lambda x: len(x.get('mutation_source', '')), 
            default=canonical_records[0]
        )
    else:  # len(canonical_records) == 0
        if records:
            selected_record = max(
                records, 
                key=lambda x: len(x.get('mutation_source', '')), 
                default=records[0]
            )
        else:
            return {col: np.nan for col in meta_list}
    
    # Retrieve necessary fields
    return {
        col: selected_record.get(col, np.nan) 
        for col in ['transcript_id', 'cds_pos', 'codon_alternate', 
                    'd_ess', 'd_ese', 'd_esr', 'splice_site_acceptor', 
                    'splice_site_donor']
    }


def process_synmall_output(df_path, original_header, selected_header, 
                           canonical_transcript_set, save_dir):
    """
    Process raw SynMall output into formatted results
    
    Parameters:
    -----------
    df_path : str
        Path to raw SynMall output file
    original_header : list
        Original column headers
    selected_header : list
        Selected column headers to keep
    canonical_transcript_set : set
        Set of canonical transcript IDs
    save_dir : str
        Directory to save output files
    
    Returns:
    --------
    pd.DataFrame
        Processed DataFrame with SynMall annotations
    """
    import numpy as np
    
    # Read data
    df = pd.read_csv(df_path, sep='\t', header=4)
    df.drop_duplicates(inplace=True)
    df.columns = original_header
    
    # Expand transcript-wise info
    extracted_cols = df['transcript_wise_info'].apply(
        lambda x: pd.Series(
            select_transcript_info(x, canonical_transcript_set, TRANSCRIPT_META_LIST)
        )
    )
    
    # Concatenate extracted columns
    df = pd.concat([df, extracted_cols], axis=1)
    
    # Generate variant ID
    df['Variant38'] = retrieve_variant(df, '#CHROM', 'POS', 'REF', 'ALT')
    df.drop_duplicates(inplace=True)
    
    # Create metadata DataFrame
    df_meta = df[['Variant38', 'transcript_id']]
    df_meta.drop_duplicates(inplace=True)
    
    # Create CDS info DataFrame
    df_cds = df[[
        '#CHROM', 'POS', 'REF', 'ALT', 'Variant38', 
        'transcript_id', 'cds_pos', 'codon_alternate'
    ]]
    df_cds.drop_duplicates(inplace=True)
    
    # Drop temporary columns
    df.drop(['transcript_id', 'cds_pos'], axis=1, inplace=True)
    
    # Select and rename columns
    available_header = [col for col in selected_header if col in df.columns]
    df = df[available_header]
    df.columns = [
        'Variant38', 'verPhyloP', 'GerpS', 'ZooPriPhyloP', 
        'ZooVerPhyloP', 'DS_AG', 'DS_AL', 'DS_DG', 'DS_DL'
    ]
    
    # Convert numeric columns
    numeric_cols = [col for col in df.columns if col != 'Variant38']
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    
    # Save results
    if save_dir:
        df.to_csv(f"{save_dir}/SynMall.output.txt", sep='\t', index=False)
        df_meta.to_csv(f"{save_dir}/input.meta.txt", sep='\t', index=False)
        df_cds.to_csv(f"{save_dir}/input.cds.txt", sep='\t', index=False)
    
    return df


# ============================================================================
# MAIN EXECUTION
# ============================================================================
print(f"[INFO] Task starts at {TASK_DIR}...")
os.makedirs(TASK_DIR, exist_ok=True)

# Annotate SynMall info
retrieve_raw_synmall(
    input_file=INPUT_FILE, 
    output_file=f"{TASK_DIR}/SynMall.raw.txt"
)

process_synmall_output(
    df_path=f"{TASK_DIR}/SynMall.raw.txt",
    original_header=original_header,
    selected_header=selected_header,
    canonical_transcript_set=canonical_transcript_set,
    save_dir=TASK_DIR
)