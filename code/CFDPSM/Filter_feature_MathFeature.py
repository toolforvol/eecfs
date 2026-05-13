import pandas as pd
import argparse


# ============================================================================
# FUNCTION
# ============================================================================
def filter_output(i1, i2, i3, output):
    """
    Merge SE and CCGR features and calculate SE3 difference
    
    Parameters:
    -----------
    i1 : str
        SE reference feature file path
    i2 : str
        SE alternate feature file path
    i3 : str
        CCGR feature file path
    output : str
        Output CSV file path
    """
    # Load SE reference features
    df1 = pd.read_csv(i1, skiprows=1)
    df1 = df1.iloc[:, :-1]  # Remove label column
    df1.columns = ['Variant38'] + [f'SE{i}_ref' for i in range(1, 6)]
    df1.drop_duplicates(inplace=True)
    
    # Load SE alternate features
    df2 = pd.read_csv(i2, skiprows=1)
    df2 = df2.iloc[:, :-1]  # Remove label column
    df2.columns = ['Variant38'] + [f'SE{i}_alt' for i in range(1, 6)]
    df2.drop_duplicates(inplace=True)
    
    # Load CCGR features
    df3 = pd.read_csv(i3, header=None)
    df3 = df3.iloc[:, :-1]  # Remove label column
    df3.columns = ['Variant38'] + [f'CCGR{i}_alt' for i in range(1, df3.shape[1])]
    df3.drop_duplicates(inplace=True)
    
    # Merge dataframes
    print(f"[INFO] Before merge: {len(df1)}")
    df = df1.merge(df2, how='left', on='Variant38')
    print(f"[INFO] After df2 merge: {len(df)}")
    df = df.merge(df3, how='left', on='Variant38')
    print(f"[INFO] After df3 merge: {len(df)}")
    
    # Calculate SE3 difference (alt - ref)
    df['SE3_diff'] = df['SE3_alt'] - df['SE3_ref']
    
    # Select output columns
    selected_columns = ['Variant38', 'SE3_diff', 'CCGR168_alt', 'CCGR166_alt']
    df = df[selected_columns]
    df.to_csv(output, index=False)
    print(f"[INFO] Results saved to {output}")


# ============================================================================
# MAIN
# ============================================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge features from MathFeature")
    parser.add_argument("--i1", type=str, required=True, help="SE reference feature file")
    parser.add_argument("--i2", type=str, required=True, help="SE alternate feature file")
    parser.add_argument("--i3", type=str, help="CCGR feature file")
    parser.add_argument("--output", type=str, required=True, help="Output CSV file path")
    
    args = parser.parse_args()
    
    filter_output(
        i1=args.i1,
        i2=args.i2,
        i3=args.i3,
        output=args.output
    )