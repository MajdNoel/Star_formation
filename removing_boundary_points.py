## This code is set to remove the boundary points from the cube, so they do not affect the averaging

import pandas as pd
import numpy as np
import os
Resolution = 12


## I  this function, we will remove all the boundaries of the cube, so we can do a better finite difference method.
def removing_boundaries(df, Resolution_K):
    counterr = 0
    max = np.max(np.array(df['X']))
    min = np.min(np.array(df['X']))
    for i in range(len(df['X'])):
        if df['X'][i] == max or df['X'][i] == min or df['Y'][i] == max or df['Y'][i] == min or df['Z'][i] == max or df['Z'][i] == min:
            df = df.drop(index=i)
    return df

#directory = "/Users/majd_noel/Documents/Western_University_Overgraduate_studies/Academics/Research/Python_project/Driving_forces" # directory of the file
file_name = "ALL_and_net_forces_cgs_L_" + str(Resolution)+ ".xlsx" # the name of the file. It has to be in the extension of .xlsx
# Loading the file using pandas. The file has
df_origional = pd.read_excel(file_name)
#vv = load_file(file_directory, f_Density_name)
v = removing_boundaries(df_origional, Resolution)
v.to_excel('boundaries_Removed_L_' + str(Resolution) + '.xlsx', index=False, engine='openpyxl')