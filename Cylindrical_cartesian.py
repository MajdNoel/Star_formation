## This code calculates switches the cartesian coordinates to cylindrical coordinates

## Importing the modules and libraries
import numpy as np
import pandas as pd

# pass dictionary as data for the argument
## Transferring the cartesian coordinates to cylindrical coordinates

def cartesian_to_cylindrical(data, force_type, Resolution_L): # This function takes data argument and force type argument.
    # 'data' should be a dictionary, with X, Y and Z and 'force_type' should be a string
    'The data should have the positions as X,Y and Z along with L. The forces, f_x, f_y and f_z'

    # other quantities, the vector fields

    # calculating the azimuthal angle
    phi_arctan2 = np.arctan2(data['Y'], data['X'])  # We use arctan2 to handle quadrants correctly. The azimuthal coordinate
    Phi = np.mod(phi_arctan2, 2 * np.pi) # ensure that the angle is within the range of 0 to 2pi

    # Calculating the radial distance
    r = np.sqrt(data['X']**2 + data['Y']**2)

    # The height stays the same, so we will just keep it as it is

    lists_of_variables = []
    # This loop is to convert the strings to the names f_r, f_phi and f_z

    # Loop through each column in the DataFrame
    for col in data.columns:
        if col not in ['X', 'Y', 'Z', 'L']:  # Exclude these columns
            lists_of_variables.append(np.array(data[col]))  # Append the remaining columns
    '''Make sure that in the excel, file, the force in the cartesian coordinate are on the order of f_X -> f_Y -> f_Z'''
    f_x = np.array(lists_of_variables[0])
    f_y = np.array(lists_of_variables[1])
    f_z = np.array(lists_of_variables[2])

    # Calculating the function in the cylindrical coordinates
    f_r = (f_x * np.cos(Phi)) + (f_y * np.sin(Phi)) # the function in the radial coordinate
    f_phi =  (f_y * np.cos(Phi)) - (f_x * np.sin(Phi)) # The function the azimuthal coordinate
    # the function in the z component is already listed above. It does not need any transformation
    f_magnitude = np.sqrt((f_x**2) +(f_y**2) + (f_z**2)) # Computing the magnitude in cylindrical coordinates.
    # This should be the same thing as the magnitude in cartesian coordinates

    ## Adding all the parameters to a dictionary
    data_cylind = {
        'r': r,
        'phi': Phi,
        'z': data['Z'],
        'L': data['L'],
        force_type + '_r': f_r,
        force_type + '_phi': f_phi,
        force_type + '_z': f_z,
        force_type + '_magnitude': f_magnitude
    }

    ## just for python precession
    data_cylind['r'] = data_cylind['r'].astype(np.float64)
    data_cylind['phi'] = data_cylind['phi'].astype(np.float64)
    data_cylind['z'] = data_cylind['z'].astype(np.float64)
    data_cylind['L'] = data_cylind['L'].astype(np.float64)
    data_cylind[force_type + '_r'] = data_cylind[force_type + '_r'].astype(np.float64)
    data_cylind[force_type + '_phi'] = data_cylind[force_type + '_phi'].astype(np.float64)
    data_cylind[force_type + '_z'] = data_cylind[force_type + '_z'].astype(np.float64)
    data_cylind[force_type + '_magnitude'] = data_cylind[force_type + '_magnitude'].astype(np.float64)

    # Saving and creating a file in cylindrical coordinates
    df_cylind = pd.DataFrame(data_cylind)
    df_cylind.to_csv('cylindrical_' + force_type +'_data_cgs_L_' + str(Resolution_L) +'.txt', index=False, sep='\t')
    return df_cylind # The code is going to return a data frame that has functions in each direction along with the positions and L

