import pandas as pd
import os


def search_dat_file(filename, search_value):

    # Construct the absolute path to the file in the "data/Polars" folder
    filepath = os.path.join(r"C:\Users\tomva\pythonProject\DUUC\data\Polars",
                            filename)

    # Normalize the filepath to ensure it resolves correctly
    filepath = os.path.abspath(filepath)

    if not os.path.exists(filepath):
        print(f"File {filename} not found in {filepath}.")
        return None

    # Read the .dat file assuming it has multiple columns
    df = pd.read_csv(filepath, sep=r'\s+', header=None)

    # Search for the value in the first column
    result = df[df[0] == search_value]

    if not result.empty:
        # search alpha return cl, cdp and cm
        return result.iloc[0, [1, 3, 4]].values
    else:
        print("Value not found.")
        return None
