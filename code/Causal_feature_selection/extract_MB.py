import pandas as pd


def retrieve_markov_blanket(input_file, mb_indice_path, output_file):
    """
    Retrieve MB indices from file and extract corresponding columns.
    """

    # Read MB indices
    with open(mb_indice_path) as f:
        mb_indice = list(map(int, f.read().split()))

    # First 5 columns are meta-info
    start_col_index = 4

    # MB index -> Python index
    real_indices = [i + start_col_index - 1 for i in mb_indice]

    # Read data
    df = pd.read_csv(input_file)

    # Extract meta + selected features
    info_cols = df.iloc[:, :5]
    selected_features = df.iloc[:, real_indices]

    # Save result
    final_df = pd.concat([info_cols, selected_features], axis=1)
    final_df.to_csv(output_file, index=False)


def main():

    molecular_type = "DNA"
    alg_name = "EECFS"

    # input_file = f'../../data/{molecular_type}/train_2362_filled_data_{molecular_type}.csv'
    input_file_train = f'../../data/{molecular_type}/train_2362_filled_data_{molecular_type}_toy.csv'
    input_file_test = f'../../data/{molecular_type}/test_238_filled_data_{molecular_type}_toy.csv'
    output_file_train = f'../../result/{molecular_type}/{alg_name}_train_2362_selectedMB.csv'
    output_file_test = f'../../result/{molecular_type}/{alg_name}_test_238_selectedMB.csv'
    mb_indices = f'../../result/{molecular_type}/{alg_name}_MB_indice.txt'

    # train
    retrieve_markov_blanket(
        input_file=input_file_train,
        mb_indice_path=mb_indices,
        output_file=output_file_train
    )
    # test
    retrieve_markov_blanket(
        input_file=input_file_test,
        mb_indice_path=mb_indices,
        output_file=output_file_test
    )

if __name__ == "__main__":
    main()