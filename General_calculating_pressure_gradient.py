## Importing the modules and libraries:

import pandas as pd
import os
import numpy as np
'''
Resolution_L = 12 # The resolution of the file. This the zoomed in file
file_Tempreature = 'Temp_s'+str(Resolution_L)+'_excel.xltx' # the name of the file
file_density = 'Density_s'+str(Resolution_L)+'__excel_.xltx' # This is the file that we will use to get the number density.

# opening the file and assign it to be a pandas data frame
df = pd.read_excel(file_Tempreature) # the data frame 'df' of the data
df_density = pd.read_excel(file_density) # the data frame 'df' of the data


Density = np.array(df_density['Density'])
n = Density/(2.3*(1.67*(10**-24)))
df_density['n'] = n
df['n'] = n # adding it to the temperature file

df_density.to_csv("number_density_file_L_" + str(Resolution_L) + ".txt", index=False, sep='\t')

df_pressure = df

# What we need for these functions is 'df_pressure' which has the temperature and density.
# We also need to identify all the constants, like background temperature and so on.
'''

def compute_pressure_gradient(df_pressure, component, B_T, size, number_density, Boltzman_constant): # The component should be a string
    #global df_pressure # Entering the pressure Value from above
    choices = ['X', 'Z', 'Y']  # The choices for the coordinates. It is going to be used later in the code
    choices.remove(component) # removing the component from the choice list
    component_2 = choices[0]  # choosing which coordinate will change the second fastest, and thr slowest
    component_3 = choices[1]
    df_pressure = df_pressure.sort_values(by=[component_3, component_2, component]) # sorting the values from
    presssure_col = np.array(df_pressure['Pressure'])
    #n_col = n
    # the fastest to the slowest in the change
    P_component_list_ = [(presssure_col[1] - (B_T * number_density[0] * Boltzman_constant)) / (-2 * size)] # boundary of the cube
    counter = 1 # this is a counter for the loop
    for j in presssure_col[1:]:
        if counter == len(presssure_col) - 1:
            P_component_list_.append(((B_T * number_density[-1] * Boltzman_constant) - presssure_col[counter - 1]) / (-2 * size))
            break # if it reaches the boundary of the cube, break out of the loop
        P_component_list_.append((presssure_col[counter + 1] - presssure_col[counter - 1]) / (-2 * size))
        counter = counter + 1
    df_pressure['Grad_P_' + str(component)] = P_component_list_ # Adding the gradient list to the data frame
    # the name of the coloumn is going to be depending on the component name
    return df_pressure # The returned data from the function

