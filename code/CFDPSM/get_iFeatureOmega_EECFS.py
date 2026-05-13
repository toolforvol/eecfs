import pandas as pd
import iFeatureOmegaCLI
import os
import argparse
import sys
from Bio import SeqIO


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================
def get_vid_from_fasta(fasta_file):
    """Extract variant IDs from FASTA file"""
    vid_list = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        vid = str(record.id).split('|')[0]
        vid_list.append(vid)
    return vid_list


def extract_features(inputfile, feature_list, molecular):
    """
    Extract features using iFeatureOmegaCLI
    
    Parameters:
    -----------
    inputfile : str
        Path to input FASTA file
    feature_list : list
        List of feature names to extract
    molecular : str
        Type of molecule ('DNA', 'RNA', or 'protein')
    
    Returns:
    --------
    pd.DataFrame
        Extracted features
    """
    # Validate input file
    if not os.path.exists(inputfile):
        print(f"[Error] Input file {inputfile} does not exist", file=sys.stderr)
        return None
    
    # Initialize feature extractor based on molecule type
    param_file = ""
    if molecular == 'DNA':
        fea_extractor = iFeatureOmegaCLI.iDNA(inputfile)
        param_file = "./parameters/DNA_parameters_setting.json"
        print("[Success] Built iDNA object")
    elif molecular == 'RNA':
        fea_extractor = iFeatureOmegaCLI.iRNA(inputfile)
        param_file = "./parameters/RNA_parameters_setting.json"
        print("[Success] Built iRNA object")
    else:  # protein
        fea_extractor = iFeatureOmegaCLI.iProtein(inputfile)
        param_file = "./parameters/Protein_parameters_setting.json"
        print("[Success] Built iProtein object")
    
    if fea_extractor is None:
        print(f"[Error] Cannot create {molecular} object for {inputfile}", file=sys.stderr)
        return None
    
    # Load parameters
    if not os.path.exists(param_file):
        print(f"[Error] Parameter file {param_file} does not exist", file=sys.stderr)
        return None
    
    fea_extractor.import_parameters(param_file)
    
    # Extract features
    df = pd.DataFrame()
    for feature in feature_list:
        try:
            print(f"\n[INFO] Processing: {feature}")
            fea_extractor.get_descriptor(feature)
            
            if fea_extractor.encodings is None:
                print(f"[Warning] Feature {feature} returned None", file=sys.stderr)
                continue
            
            print(f"[INFO] Shape of {feature}: {fea_extractor.encodings.shape}")
            print(fea_extractor.encodings.head())
            df = pd.concat([df, fea_extractor.encodings], axis=1)
            
        except Exception as e:
            print(f"[Error] Processing {feature}: {str(e)}", file=sys.stderr)
            continue
    
    if df.empty:
        print("[Error] No features were extracted", file=sys.stderr)
        return None
    
    return df


def merge_ref_alt_features(df_alt, df_ref):
    """
    Merge reference and alternate features with differences
    
    Parameters:
    -----------
    df_alt : pd.DataFrame
        Alternate allele features
    df_ref : pd.DataFrame
        Reference allele features
    
    Returns:
    --------
    pd.DataFrame
        Merged DataFrame with ref, alt, and diff columns
    """
    # Convert to numeric
    df_alt = df_alt.apply(pd.to_numeric, errors='coerce')
    df_ref = df_ref.apply(pd.to_numeric, errors='coerce')
    
    # Align columns
    columns = df_alt.columns.tolist()
    df_ref = df_ref.reindex(columns=columns)
    
    # Calculate differences
    df_diff = df_alt - df_ref
    
    # Rename columns
    df_ref.columns = [f"{col}_ref" for col in columns]
    df_alt.columns = [f"{col}_alt" for col in columns]
    df_diff.columns = [f"{col}_diff" for col in columns]
    
    # Concatenate
    df_result = pd.concat([df_ref, df_alt, df_diff], axis=1)
    return df_result


def get_ifo_selected(ref_fasta_file, alt_fasta_file, output_file, molecular, target_features, selected_features):
    """
    Main function to extract and select features
    
    Parameters:
    -----------
    ref_fasta_file : str
        Reference FASTA file path
    alt_fasta_file : str
        Alternate FASTA file path
    output_file : str
        Output CSV file path
    molecular : str
        Molecule type (DNA/RNA/protein)
    target_features : list
        List of features to extract
    selected_features : list
        List of features to keep in output
    """
    # Extract variant IDs
    vid_list = get_vid_from_fasta(ref_fasta_file)
    
    # Extract reference features
    ref_df = extract_features(ref_fasta_file, target_features, molecular)
    if ref_df is None:
        return
    
    # Extract alternate features (skip for protein)
    if molecular.lower() != 'protein':
        alt_df = extract_features(alt_fasta_file, target_features, molecular)
        if alt_df is None:
            return
        merge_df = merge_ref_alt_features(alt_df, ref_df)
    else:
        merge_df = ref_df.copy()
    
    # Select and reorder columns
    merge_df = merge_df[selected_features]
    merge_df.insert(0, 'Variant38', vid_list)
    merge_df.to_csv(output_file, index=False)
    print(f"[INFO] Results saved to {output_file}")


# ============================================================================
# MAIN
# ============================================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Feature Extraction for iFeatureOmega")
    parser.add_argument("--ref_fasta_file", type=str, required=True, help="Reference FASTA file")
    parser.add_argument("--alt_fasta_file", type=str, help="Alternate FASTA file")
    parser.add_argument("--output_file", type=str, required=True, help="Output CSV file path")
    parser.add_argument("--molecular", type=str, required=True, help="Molecule type: DNA | RNA | protein")
    parser.add_argument("--target_features", nargs="+", required=True, help="Standard descriptor names to extract")
    parser.add_argument("--selected_features", nargs="+", required=True, help="Specific features to keep in output")
    
    args = parser.parse_args()
    
    get_ifo_selected(
        ref_fasta_file=args.ref_fasta_file,
        alt_fasta_file=args.alt_fasta_file,
        output_file=args.output_file,
        molecular=args.molecular,
        target_features=args.target_features,
        selected_features=args.selected_features
    )