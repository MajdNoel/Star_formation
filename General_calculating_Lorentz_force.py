## This is a code to calculate the lorentz force.


# importing thr libraries and modules
import pandas as pd
import os
import numpy as np



# Creating a function to calculate the Lorentz force:
def compute_Lorentz_force(df_magnetic, size, Resolution_L):
    # Putting all the magnetic field values as arrays
    B_x_col = np.array(df_magnetic['B_x'])
    B_y_col = np.array(df_magnetic['B_y'])
    B_z_col = np.array(df_magnetic['B_z'])

    # For the Lorentz force, we will not calculate the boundary as of now. It is very hard to do that. but also, we will
    # not include that in the plotting
    # Therefore, we will make a list of items that starts with 0

    terms = {} # Creating an empty dictionary
    for i in range(1,13):
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


    for j in B_x_col[1:]:
        if (counter_group_1 == len(B_x_col) - 1):
            terms["Term_2_list"].append(0)
            terms["Term_3_list"].append(0)
            terms["Term_5_list"].append(0)
            terms["Term_12_list"].append(0)
            break  # if it reaches the boundary of the cube, break out of the loop
        terms["Term_2_list"].append(B_z_col[counter_group_1]* ((B_z_col[counter_group_1 + 1] - B_z_col[counter_group_1 - 1]) / (2 * size)))
        terms["Term_3_list"].append(B_y_col[counter_group_1] * ((B_y_col[counter_group_1 + 1] - B_y_col[counter_group_1 - 1]) / (2 * size)))
        terms["Term_5_list"].append(B_x_col[counter_group_1] * ((B_y_col[counter_group_1 + 1] - B_y_col[counter_group_1 - 1]) / (2 * size)))
        terms["Term_12_list"].append(B_x_col[counter_group_1] * ((B_z_col[counter_group_1 + 1] - B_z_col[counter_group_1 - 1]) / (2 * size)))
        counter_group_1 = counter_group_1 + 1

        # appending all the results to the data frame
    df_magnetic['Bz_dBz_dx'] = terms["Term_2_list"]  # Adding the gradient list to the data frame
    df_magnetic['By_dBy_dx'] = terms["Term_3_list"]  # Adding the gradient list to the data frame
    df_magnetic['Bx_dBy_dx'] = terms["Term_5_list"]  # Adding the gradient list to the data frame
    df_magnetic['Bx_dBz_dx'] = terms["Term_12_list"]  # Adding the gradient list to the data frame



    ## We sort the values of the Excel file (that we just appended the new terms into, so when the Excel files is sorted
    # out to be changing the fastest in Y, all the new terms stick with grid cells they were assigned to) to be changing
    # in the Y direction as the fastest
    df_magnetic = df_magnetic.sort_values(by=['X', 'Z', 'Y'])  # sorting the values from

    # We define all the data frame again, to use them in the loop. Because we sorted the original data, so now we need to
    # define our new columns
    B_x_col = np.array(df_magnetic['B_x'])
    B_y_col = np.array(df_magnetic['B_y'])
    B_z_col = np.array(df_magnetic['B_z'])


    counter_group_1 = 1
    for jj in B_y_col[1:]:
        if (counter_group_1 == len(B_y_col) - 1):
            terms["Term_4_list"].append(0)
            terms["Term_6_list"].append(0)
            terms["Term_7_list"].append(0)
            terms["Term_9_list"].append(0)
            break  # if it reaches the boundary of the cube, break out of the loop
        terms["Term_4_list"].append(B_y_col[counter_group_1] * ((B_x_col[counter_group_1 + 1] - B_x_col[counter_group_1 - 1]) / (2 * size)))
        terms["Term_6_list"].append(B_x_col[counter_group_1] * ((B_x_col[counter_group_1 + 1] - B_x_col[counter_group_1 - 1]) / (2 * size)))
        terms["Term_7_list"].append(B_z_col[counter_group_1] * ((B_z_col[counter_group_1 + 1] - B_z_col[counter_group_1 - 1]) / (2 * size)))
        terms["Term_9_list"].append(B_y_col[counter_group_1] * ((B_z_col[counter_group_1 + 1] - B_z_col[counter_group_1 - 1]) / (2 * size)))
        counter_group_1 = counter_group_1 + 1

    ## Appending everything to the list to the data frame

    df_magnetic['By_dBx_dy'] = terms["Term_4_list"]
    df_magnetic['Bx_dBx_dy'] = terms["Term_6_list"]  # Adding the gradient list to the data frame
    df_magnetic['Bz_dBz_dy'] = terms["Term_7_list"]  # Adding the gradient list to the data frame
    df_magnetic['By_dBz_dy'] = terms["Term_9_list"]  # Adding the gradient list to the data frame


    ## We sort the values of the Excel file (that we just appended the new terms into, so when the Excel files is sorted
    # out to be changing the fastest in Z, all the new terms stick with grid cells they were assigned to) to be changing
    # in the Z direction as the fastest
    df_magnetic = df_magnetic.sort_values(by=['X', 'Y', 'Z'])  # sorting the values from

    # We define all the data frame again, to use them in the loop. Because we sorted the original data, so now we need to
    # define our new columns
    B_x_col = np.array(df_magnetic['B_x'])
    B_y_col = np.array(df_magnetic['B_y'])
    B_z_col = np.array(df_magnetic['B_z'])

    counter_group_1 = 1
    for jjj in B_z_col[1:]:
        if (counter_group_1 == len(B_z_col) - 1):
            terms["Term_1_list"].append(0)
            terms["Term_8_list"].append(0)
            terms["Term_10_list"].append(0)
            terms["Term_11_list"].append(0)
            break  # if it reaches the boundary of the cube, break out of the loop
        terms["Term_1_list"].append(B_z_col[counter_group_1] * ((B_x_col[counter_group_1 + 1] - B_x_col[counter_group_1 - 1]) / (2 * size)))
        terms["Term_8_list"].append(B_z_col[counter_group_1] * ( (B_y_col[counter_group_1 + 1] - B_y_col[counter_group_1 - 1]) / (2 * size)))
        terms["Term_10_list"].append(B_y_col[counter_group_1] * ((B_y_col[counter_group_1 + 1] - B_y_col[counter_group_1 - 1]) / (2 * size)))
        terms["Term_11_list"].append(B_x_col[counter_group_1] * ((B_x_col[counter_group_1 + 1] - B_x_col[counter_group_1 - 1]) / (2 * size)))
        counter_group_1 = counter_group_1 + 1


    # Appending to the data frame list
    df_magnetic['Bz_dBx_dz'] = terms["Term_1_list"]  # Adding the gradient list to the data frame
    df_magnetic['Bz_dBy_dz'] = terms["Term_8_list"]  # Adding the gradient list to the data frame
    df_magnetic['By_dBy_dz'] = terms["Term_10_list"]  # Adding the gradient list to the data frame
    df_magnetic['Bx_dBx_dz'] = terms["Term_11_list"]  # Adding the gradient list to the data frame

    ## Sort the values again as X changes the fastest Y the second fastest, and Z the slowest (Just like the original
    # files) We do this because we want to compare the data from the other files as well.

    df_magnetic = df_magnetic.sort_values(by=['Z', 'Y', 'X'])  # sorting the values as X changes the fastest like the original

    ## Switching all the components to numpy variables
    Term_1_np = np.array(df_magnetic['Bz_dBx_dz'])
    Term_2_np = np.array(df_magnetic['Bz_dBz_dx'])
    Term_3_np = np.array(df_magnetic['By_dBy_dx'])
    Term_4_np = np.array(df_magnetic['By_dBx_dy'])
    Term_5_np = np.array(df_magnetic['Bx_dBy_dx'])
    Term_6_np = np.array(df_magnetic['Bx_dBx_dy'])
    Term_7_np = np.array(df_magnetic['Bz_dBz_dy'])
    Term_8_np = np.array(df_magnetic['Bz_dBy_dz'])
    Term_9_np = np.array(df_magnetic['By_dBz_dy'])
    Term_10_np = np.array(df_magnetic['By_dBy_dz'])
    Term_11_np = np.array(df_magnetic['Bx_dBx_dz'])
    Term_12_np = np.array(df_magnetic['Bx_dBz_dx'])

    ## Adding it to the data frame
    x_component = (1 / (4 * np.pi)) * (Term_1_np - Term_2_np - Term_3_np + Term_4_np)
    y_component = (1 / (4 * np.pi)) * (Term_5_np - Term_6_np - Term_7_np + Term_8_np)
    z_component = (1 / (4 * np.pi)) * (Term_9_np - Term_10_np - Term_11_np + Term_12_np)
    F_L_mag = np.sqrt((x_component**2)+(y_component**2)+(z_component**2))


    ## Adding it to the data frame
    df_magnetic['L_F_x_component'] = x_component  # Adding the gradient list to the data frame
    df_magnetic['L_F_y_component'] = y_component  # Adding the gradient list to the data frame
    df_magnetic['L_F_z_component'] = z_component  # Adding the gradient list to the data frame
    df_magnetic['L_F_magnitude'] = F_L_mag  # Adding the gradient list to the data frame

    ## Saving the data frame with positions, resolution, magnetic fields, gradient in every component, and the lorentz force components
    df_magnetic.to_excel('cartesian_magnetic_Lorentz_cgs_L_' + str(Resolution_L) + '.xlsx', index=False, engine='openpyxl')




    return df_magnetic # returning the data frame











