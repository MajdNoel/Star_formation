# This code is for calculating the gravitational force.
''' It will be used as one of the files to be imported to the main code to calculate the net force of the driven
outflows from the protostellar system. The simulation gets us the animation with some of the parameters as well.
For example, we have the velocity fields, the magnetic fields and so on. The main purpose is to find the
net forces on these outflows, so we can determine what forces are driving the outflows. In the simulation we
have many files for the resolution. The more the resolution, the smaller the grid cell, and the more we are zooming
in the center of the star '''

## Importing the modules and libraries:

import pandas as pd
import os
import numpy as np


# Reading and opening the file

#Resolution_L = 12 # This is the resolution. The highest resolution is 13.
#directory = "/Users/majd_noel/Documents/Western_University_Overgraduate_studies/Academics/Research/Python_project/Gravitational_force" # directory of the file
#file_name = "Density_s" + str(Resolution_L)+ "_excel_.xlsx" # the name of the file. The file has to be in the extension of xltx

# To run the above function, the directory of the file needs to be specified, along with the file name.

# Defining all the constants and parameters in cgs units
#G = 6.67 * 10 **(-8) # The gravitational constant in (dyne((cm^2)/(g^2)
#M_ps = 0.32542335 * 1.989 * 10**33 # mass of the protostar in g. This will depend on the stage of the simulation.
# The mass of the protostar grows in time

# This is a function that computes the gravitational force. IT returns the data frame where the gravitational forces
# along with the magnitude are calculated
def compute_gravitational_force(pandas_data_frame, Resolution_L, G, M_ps):
    ## Pulling the arrays as lists in cgs units
    df_Den = pandas_data_frame['Density']  # this is the list in g/cm^3
    df_X = pandas_data_frame['X']  # the X position in cm
    df_Y = pandas_data_frame['Y']  # the Y position in cm
    df_Z = pandas_data_frame['Z']  # the Z position in cm

    # Calculating the distance of each grid cell
    r_cell = np.sqrt(df_X ** 2 + df_Y ** 2 + df_Z ** 2)  # the radius in units of cm
    ## Calculating the forces in cgs units
    F_g_X = G * M_ps * (1 / ((r_cell) ** 3)) * (df_Den) * (- df_X)  # The force in the x component F_x in dyne
    F_g_Y = G * M_ps * (1 / ((r_cell) ** 3)) * (df_Den) * (- df_Y)  # The force in the y component F_y in dyne
    F_g_Z = G * M_ps * (1 / ((r_cell) ** 3)) * (df_Den) * (- df_Z)  # The force in the Z component F_z in dyne
    # Computing the magnitude of the forces
    F_g_magnitude = np.sqrt(F_g_X ** 2 + F_g_Y ** 2 + F_g_Z ** 2)
    # We put the negative because we sre interested in the force that the protostar exerts on the grid cell
    Gravity_forces_data_dictionary = {
        'X': pandas_data_frame['X'],
        'Y': pandas_data_frame['Y'],
        'Z': pandas_data_frame['Z'],
        'L': pandas_data_frame['L'],
        'F_g_x': F_g_X,
        'F_g_y': F_g_Y,
        'F_g_z': F_g_Z,
        'F_g_magnitude': F_g_magnitude
    }

    ## just for python precession
    Gravity_forces_data_dictionary['X'] = Gravity_forces_data_dictionary['X'].astype(np.float64)
    Gravity_forces_data_dictionary['Y'] = Gravity_forces_data_dictionary['Y'].astype(np.float64)
    Gravity_forces_data_dictionary['Z'] = Gravity_forces_data_dictionary['Z'].astype(np.float64)
    Gravity_forces_data_dictionary['L'] = Gravity_forces_data_dictionary['L'].astype(np.float64)
    Gravity_forces_data_dictionary['F_g_z'] = Gravity_forces_data_dictionary['F_g_z'].astype(np.float64)
    Gravity_forces_data_dictionary['F_g_y'] = Gravity_forces_data_dictionary['F_g_y'].astype(np.float64)
    Gravity_forces_data_dictionary['F_g_x'] = Gravity_forces_data_dictionary['F_g_x'].astype(np.float64)
    Gravity_forces_data_dictionary['F_g_magnitude'] = Gravity_forces_data_dictionary['F_g_magnitude'].astype(np.float64)


    df_forces = pd.DataFrame(Gravity_forces_data_dictionary)
    df_forces.to_excel('cartesian_Gravitational_forces_cgs_L_' + str(Resolution_L) + '.xlsx', index=False,
                       engine='openpyxl')
    return df_forces
