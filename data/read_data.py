import os
import numpy as np


def airfoil_polar(filename, search_value):

    # Construct the absolute path to the file in the "data/Polars" folder
    filepath = os.path.join(r"C:\Users\tomva\pythonProject\DUUC\data\Polars",
                            filename)

    # Normalize the filepath to ensure it resolves correctly
    filepath = os.path.abspath(filepath)

    start_line = 13
    data_polar = np.loadtxt(filepath, skiprows=start_line - 1)
    alpha = data_polar[:, 0]
    cl_polar = data_polar[:, 1]
    cd_polar = data_polar[:, 2]
    cdp_polar = data_polar[:, 3]
    cm_polar = data_polar[:, 4]

    cl_val = np.interp(search_value, alpha, cl_polar)
    cd_val = np.interp(search_value, alpha, cd_polar)
    cdp_val = np.interp(search_value, alpha, cdp_polar)
    cm_val = np.interp(search_value, alpha, cm_polar)

    cd = cdp_val + cd_val

    # index = np.where(alpha == search_value)[0]
    """
    if index.size > 0:
        cl_val = cl[index][0]
        cd_val = cd[index][0]
        cdp_val = cdp[index][0]
        cm_val = cm[index][0]
        return cl_val, cd_val, cdp_val, cm_val
    else:
        return None, print("Value not found in the file") """

    return cl_val, cd, cm_val


""" Test case for this function """
# a = airfoil_polar("pylon0012.txt", 0)
# print(a)
