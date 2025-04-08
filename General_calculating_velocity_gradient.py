## This is a code to calculate the velocity gradient.
# We have the velocity vectors in the Excel files. The process is going to be similar to calculating
# the lorentz force, by doing finite differences.

# importing thr libraries and modules
import pandas as pd
import os
import numpy as np

import os
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
#from General__calculating_gravitational_Force import compute_gravitational_force
#from File_load import load_file
from Cylindrical_cartesian import cartesian_to_cylindrical
#from General_calculating_pressure_gradient import compute_pressure_gradient
#from General_calculating_Lorentz_force import compute_Lorentz_force
#from Summation_net_forces import summation_and_net_force
from Plot_cylindrical_file import plotting_vector_field
from removing_boundary_points import removing_boundaries, Resolution

# Defining all the constants and parameters in cgs units
G = 6.67 * 10 **(-8) # The gravitational constant in (dyne((cm^2)/(g^2)
M_ps = 0.32542335 * 1.989 * 10**33 # mass of the protostar in g at this stage.
k = 1.3807 * (10**(-16)) # Boltzmann constant in cm^2, g s^-1
Background_Temperature = 10 # The temperature in Kelvin
choices = ['X', 'Y', 'Z'] # The choices for the coordinates. It is going to be used later in the code
AU_TO_CM = 1.496e13 # conversion from cm to AU
mean_molecular_weight = 2.3 # the mean molecular weight taken
mass_of_hydrogen = 1.67*(10**-24) # the mass of the hydrogen atom in g
Resolution = 12 # This is the parameter that will determine which files are we dealing with





Resolution = 12
# Creating a function to calculate the Lorentz force:
def compute_velocity_gradient(df_velocity, size, Resolution_L, df_den_col):
    # Putting all the magnetic field values as arrays
    v_x_col = np.array(df_velocity['V_x'])
    v_y_col = np.array(df_velocity['V_y'])
    v_z_col = np.array(df_velocity['V_z'])

    # For the Lorentz force, we will not calculate the boundary as of now. It is very hard to do that. but also, we will
    # not include that in the plotting
    # Therefore, we will make a list of items that starts with 0

    terms = {} # Creating an empty dictionary
    for i in range(1,10):
        terms["Term_" + str(i) + "_list"] = [0]  # Create keys that has lists to append it

    counter_group_1 = 1 # counter for the loop

    '''This is a loop to calculate the component of the Lorentz force It is very difficult to calculate the Lorentz 
    force as a general formula that applies to all the components, since each component of the vector has 
    a specific different equations from the others. We have the data as x is changing the fastest, then y and then z.
    So, whenever we have change in magnetic field over change in x, we have to use this form of data. Whenever we have 
    change in magnetic field over change in y, we have to use the form where y is changing the fastest, and similar for 
    z. In the equation for calculating the Lorentz force, we end up with 12 terms. The code groups the term that has 
    similar change in x together, the terms that have similar change in y together, and terms that have similar change
    in z together. This is to reduce the amount of steps that we change the fastest changing component. So, when x is
    Changing the fastest, we can calculate term 2, 3, 5 and 12 from the long equation that we have, and so on for the 
    other terms. It doesnt matter which two variables will change the second fastest of the slowest, as long as we are 
    passing by all of the variables (which means going acros the entire cube of data). Therefore, the code will have 3 
    groups. Each one them will calculate 4 terms, so we end up calculating all the terms in the equation'''

# Starting the loop. Since the data already starts with x as changing the fastest. This group is going to be for the
# components that has partial-partial-x.


    for j in v_x_col[1:]:
        if (counter_group_1 == len(v_x_col) - 1):
            terms["Term_1_list"].append(0)
            terms["Term_4_list"].append(0)
            terms["Term_7_list"].append(0)
            break  # if it reaches the boundary of the cube, break out of the loop
        terms["Term_1_list"].append(v_x_col[counter_group_1] * ((v_x_col[counter_group_1 + 1] - v_x_col[counter_group_1 - 1]) / (2 * size)))
        terms["Term_4_list"].append(v_x_col[counter_group_1] * ((v_y_col[counter_group_1 + 1] - v_y_col[counter_group_1 - 1]) / (2 * size)))
        terms["Term_7_list"].append(v_x_col[counter_group_1] * ((v_z_col[counter_group_1 + 1] - v_z_col[counter_group_1 - 1]) / (2 * size)))
        counter_group_1 = counter_group_1 + 1

        # appending all the results to the data frame
    df_velocity['Vx_dVx_dx'] = terms["Term_1_list"]  # Adding the gradient list to the data frame
    df_velocity['Vx_dVy_dx'] = terms["Term_4_list"]  # Adding the gradient list to the data frame
    df_velocity['Vx_dVz_dx'] = terms["Term_7_list"]  # Adding the gradient list to the data frame



    ## We sort the values of the Excel file (that we just appended the new terms into, so when the Excel files is sorted
    # out to be changing the fastest in Y, all the new terms stick with grid cells they were assigned to) to be changing
    # in the Y direction as the fastest
    df_velocity = df_velocity.sort_values(by=['X', 'Z', 'Y'])  # sorting the values from

    # We define all the data frame again, to use them in the loop. Because we sorted the original data, so now we need to
    # define our new columns
    v_x_col = np.array(df_velocity['V_x'])
    v_y_col = np.array(df_velocity['V_y'])
    v_z_col = np.array(df_velocity['V_z'])


    counter_group_1 = 1
    for jj in v_y_col[1:]:
        if (counter_group_1 == len(v_y_col) - 1):
            terms["Term_2_list"].append(0)
            terms["Term_5_list"].append(0)
            terms["Term_8_list"].append(0)
            break  # if it reaches the boundary of the cube, break out of the loop
        terms["Term_2_list"].append(v_y_col[counter_group_1] * ((v_x_col[counter_group_1 + 1] - v_x_col[counter_group_1 - 1]) / (2 * size)))
        terms["Term_5_list"].append(v_y_col[counter_group_1] * ((v_y_col[counter_group_1 + 1] - v_y_col[counter_group_1 - 1]) / (2 * size)))
        terms["Term_8_list"].append(v_y_col[counter_group_1] * ((v_z_col[counter_group_1 + 1] - v_z_col[counter_group_1 - 1]) / (2 * size)))
        counter_group_1 = counter_group_1 + 1

    ## Appending everything to the list to the data frame

    df_velocity['Vy_dVx_dy'] = terms["Term_2_list"]
    df_velocity['Vy_dVy_dy'] = terms["Term_5_list"]  # Adding the gradient list to the data frame
    df_velocity['Vy_dVz_dy'] = terms["Term_8_list"]  # Adding the gradient list to the data frame


    ## We sort the values of the Excel file (that we just appended the new terms into, so when the Excel files is sorted
    # out to be changing the fastest in Z, all the new terms stick with grid cells they were assigned to) to be changing
    # in the Z direction as the fastest
    df_velocity = df_velocity.sort_values(by=['X', 'Y', 'Z'])  # sorting the values from

    # We define all the data frame again, to use them in the loop. Because we sorted the original data, so now we need to
    # define our new columns
    v_x_col = np.array(df_velocity['V_x'])
    v_y_col = np.array(df_velocity['V_y'])
    v_z_col = np.array(df_velocity['V_z'])

    counter_group_1 = 1
    for jjj in v_z_col[1:]:
        if (counter_group_1 == len(v_z_col) - 1):
            terms["Term_3_list"].append(0)
            terms["Term_6_list"].append(0)
            terms["Term_9_list"].append(0)
            break  # if it reaches the boundary of the cube, break out of the loop
        terms["Term_3_list"].append(v_z_col[counter_group_1] * ((v_x_col[counter_group_1 + 1] - v_x_col[counter_group_1 - 1]) / (2 * size)))
        terms["Term_6_list"].append(v_z_col[counter_group_1] * ((v_y_col[counter_group_1 + 1] - v_y_col[counter_group_1 - 1]) / (2 * size)))
        terms["Term_9_list"].append(v_z_col[counter_group_1] * ((v_z_col[counter_group_1 + 1] - v_z_col[counter_group_1 - 1]) / (2 * size)))
        counter_group_1 = counter_group_1 + 1


    # Appending to the data frame list
    df_velocity['Vz_dVx_dz'] = terms["Term_3_list"]  # Adding the gradient list to the data frame
    df_velocity['Vz_dVy_dz'] = terms["Term_6_list"]  # Adding the gradient list to the data frame
    df_velocity['Vz_dVz_dz'] = terms["Term_9_list"]  # Adding the gradient list to the data frame

    ## Sort the values again as X changes the fastest Y the second fastest, and Z the slowest (Just like the original
    # files) We do this because we want to compare the data from the other files as well.

    df_velocity = df_velocity.sort_values(by=['Z', 'Y', 'X'])  # sorting the values as X changes the fastest like the original

    ## Switching all the components to numpy variables
    Term_1_np = np.array(df_velocity['Vx_dVx_dx'])
    Term_2_np = np.array(df_velocity['Vy_dVx_dy'])
    Term_3_np = np.array(df_velocity['Vz_dVx_dz'])
    Term_4_np = np.array(df_velocity['Vx_dVy_dx'])
    Term_5_np = np.array(df_velocity['Vy_dVy_dy'])
    Term_6_np = np.array(df_velocity['Vz_dVy_dz'])
    Term_7_np = np.array(df_velocity['Vx_dVz_dx'])
    Term_8_np = np.array(df_velocity['Vy_dVz_dy'])
    Term_9_np = np.array(df_velocity['Vz_dVz_dz'])

    ## Adding it to the data frame
    x_component =  df_den_col*(Term_1_np + Term_2_np + Term_3_np)
    y_component =  df_den_col*(Term_4_np + Term_5_np + Term_6_np)
    z_component =  df_den_col*(Term_7_np + Term_8_np + Term_9_np)
    V_mag = np.sqrt((x_component**2)+(y_component**2)+(z_component**2))


    ## Adding it to the data frame
    df_velocity['Grad_vx_component'] = x_component  # Adding the gradient list to the data frame
    df_velocity['Grad_vy_component'] = y_component  # Adding the gradient list to the data frame
    df_velocity['Grad_vz_component'] = z_component  # Adding the gradient list to the data frame
    df_velocity['Grad_V_magnitude'] = V_mag  # Adding the gradient list to the data frame

    ## Saving the data frame with positions, resolution, magnetic fields, gradient in every component, and the lorentz force components
    df_velocity.to_excel('Velocity_gradient_vector_field_cgs_L_' + str(Resolution_L) + '.xlsx', index=False, engine='openpyxl')
    return df_velocity # returning the data frame



file_density = 'Density_s'+str(Resolution)+'_excel.xlsx'
df_den = pd.read_excel(file_density)
df_den_col = np.array(df_den['Density'])
file_name = "V_field_s" + str(Resolution)+ "_excel.xlsx" # the name of the file. It has to be in the extension of .xlsx
df_or = pd.read_excel(file_name)
dff = pd.DataFrame(df_or)
## Calling the function
df_X_SI = dff['X'] # the X position in cm
for i in range(len(df_X_SI)): # This is a loop to calculate the grid size.
    cell_size = abs(df_X_SI.iloc[i+1] - df_X_SI.iloc[i]) # the size of the grid cell in m
    if cell_size != 0.0:
        break

vel = compute_velocity_gradient(dff,cell_size,Resolution, df_den_col)



# The dictionary for plotting

All_force_data_dictionary = {
    'X': vel['X'],
    'Y': vel['Y'],
    'Z': vel['Z'],
    'L': vel['L'],
    'F_x': vel['Grad_vx_component'],
    'F_y': vel['Grad_vy_component'],
    'F_z': vel['Grad_vz_component'],
    'F_magnitude': vel['Grad_V_magnitude']
}

All_force_data = pd.DataFrame(All_force_data_dictionary) # Transferring the dictionary into a data frame
# First, in order to plot the data, we need to transfer the valuables in cylindrical coordinates
#All_force_data.to_excel('All_force_boundaries_Removed_cartesian_cgs_L_' + str(Resolution) + '.xlsx', index=False, engine='openpyxl')
## Removing the boundaries:
the_velocity_grad = removing_boundaries(All_force_data, Resolution)



## Calling the cartesian to cylindrical function to switch the cartesian data into cylindrical
Net_force_Cylindrical_df = cartesian_to_cylindrical(the_velocity_grad, 'Velocity_gradient_fields', Resolution)
# Saving the data frame
#Net_force_Cylindrical_df.to_excel('Net_forces_cylindrical_coordinates_cgs_L_' + str(Resolution) + '.xlsx', index=False, engine='openpyxl')

# Taking care of the boundaries. We will create a data frame to exclude all the boundary points, because it is hard
# to calculate

##In order for use to plot the function, we need take the data and average it on the angle
Averaged_data_frame = Net_force_Cylindrical_df.groupby(['r', 'z']).mean().reset_index() # Averaging the quantities over the azimuthal direction

# Saving the data frame of averaged data
#Averaged_data_frame.to_excel('Net_force_averaged_cylindrical_coordinates_L_' + str(Resolution) + '.xlsx', index=False, engine='openpyxl')

# Recalling the data frame into another one
Averaged_df = Averaged_data_frame

#Plotting the function.
# To have good visualization, we need to transfer the r and z to AU, and take scale the forces by their magnitude

Plotting_variables = plotting_vector_field(Averaged_df, AU_TO_CM, Resolution)
fig1, fig2 = Plotting_variables  # Now both figures are stored












