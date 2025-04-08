'''
                                                 MAIN CODE TO RUN
__________________________________
PHYSICS BEHIND THE CODE:
The main purpose of the code is find the driving forces for the outflows from protostellar systems.
It is still not very understood how these outflows gets ejected from the protostellar systems. However, for sometime
now, the main idea became that ejections have been associated with magnetic fields. the current understanding is that
magnetic fields can extract angular momentum form the system and be able to launch these outflows. In this research,
we are trying to find the net forces on these outflows, to determine the driving forces on these outflows, and explain
their origin and ejections. The magnetic fields, velocity fields, Temperature, density and the rest of the parameters
were and published in the simulation that is mentioned in this paper: https://arxiv.org/abs/2401.04260.
The forces that are considered are the gravitational force, the thermal pressure gradient, and the Lorentz force
__________________________________
INPUT AND PROCESS:
This main code imports several files that should be in the same directory as this file. All of these files have
functions that complete each other are used in this code. All the files should be in Excel '.xlsx' format.
There are files for:
    - Calculating gravitational force
    - Calculating pressure gradient force
    - Calculating Lorentz force
    - Computing the Net force
    - Transform the data from Cartesian to Cylindrical coordinates
    - Plotting the data
The simulation has a cube of grid cells 64x64x64. There are different set of resolution, which makes the grid cell
smaller, and it is like zooming into the center of the cube, where the protostar is.
the code uses the method of finite differences of the second order to calculate the gradient of the variables.
For plotting (for now), the code is using azimuthal angle averaged over the quantities, so we can plot this in the
r-z plane.
__________________________________
OUTPUT AND RESULTS:
The code could output plots for each individual force and the net force in total, along with saving files in the same
directory of the code that carries all the data for the forces, in cartesian and cylindrical coordinates
__________________________________
NEXT STEPS:
For the next steps as of now:
    - Plotting the left side of the Euler equation
    - Getting the data for thr previous times to see if the velocity fields are changing with time at this stage?
    - Plot all the forces for the previous time steps to see how the forces before this stage and compare it with current
      velocity fields (since velocity takes time to go along the acceleration (other than in the centripetal
      acceleration case))
    - Plot all the vector fields in 3D, to see the forces shape
__________________________________
'''

# Importing the libraries and modules

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from General__calculating_gravitational_Force import compute_gravitational_force
from File_load import load_file
from Cylindrical_cartesian import cartesian_to_cylindrical
from General_calculating_pressure_gradient import compute_pressure_gradient
from General_calculating_Lorentz_force import compute_Lorentz_force
from Summation_net_forces import summation_and_net_force
from Plot_cylindrical_file import plotting_vector_field
from removing_boundary_points import removing_boundaries


# Defining all the constants and parameters in cgs units
G = 6.67 * 10 **(-8) # The gravitational constant in (dyne((cm^2)/(g^2)
M_ps = 0.32542335 * 1.989 * 10**33 # mass of the protostar in g at this stage.
k = 1.3807 * (10**(-16)) # Boltzmann constant in cm^2, g s^-1
Background_Temperature = 10 # The temperature in Kelvin
choices = ['X', 'Y', 'Z'] # The choices for the coordinates. It is going to be used later in the code
AU_TO_CM = 1.496e13 # conversion from cm to AU
mean_molecular_weight = 2.3 # the mean molecular weight taken
mass_of_hydrogen = 1.67*(10**-24) # the mass of the hydrogen atom in g
Resolution = 11 # This is the parameter that will determine which files are we dealing with

'''_________________________________________________________________________________________________________
                                        CALCULATING THE GRAVITATIONAL FORCE
_________________________________________________________________________________________________________'''
## Importing the needed files:
file_directory = "/Users/majd_noel/Documents/Western_University_Overgraduate_studies/Academics/Research/Python_project/ALl_files" # directory of the file
f_Density_name = "Density_s" + str(Resolution)+ "_excel.xlsx" # the name of the file. It has to be in the extension of .xlsx
## Calling the function
Density_df = load_file(file_directory, f_Density_name) # Density_df is a file that has the density variable in it
## Calculating the gravity and creating a file of it:
Gravity_df = compute_gravitational_force(Density_df, Resolution, G, M_ps) # Gravity_df is a data frame that has positions, and F_g in cartesian coordinates

'''_________________________________________________________________________________________________________
                                        CALCULATING THE THERMAL PRESSURE GRADIENT
_________________________________________________________________________________________________________'''

## Importing the needed files:
file_directory = "/Users/majd_noel/Documents/Western_University_Overgraduate_studies/Academics/Research/Python_project/ALl_files" # directory of the file
file_Temperature_name = "Temperature_s" + str(Resolution)+ "_excel.xlsx" # the name of the file. It has to end in the extension of .xlsx
## Calling the function
Pressure_Temperature_density_df = load_file(file_directory, file_Temperature_name) # Density_df is a file that has the density variable in it

Density = np.array(Density_df['Density'])
n = (Density)/((mean_molecular_weight)*(mass_of_hydrogen)) # the number denisty
Pressure_Temperature_density_df['n'] = n # adding the number density to the temperature/density file

Pressure_Temperature_density_df['Pressure'] = Pressure_Temperature_density_df['Temp']* n * k
# Now, we have a file that has the temperature and number density for each grid
## Calculating the cell size:


df_X_SI = Pressure_Temperature_density_df['X'] # the X position in cm
for i in range(len(df_X_SI)): # This is a loop to calculate the grid size.
    cell_size = abs(df_X_SI.iloc[i+1] - df_X_SI.iloc[i]) # the size of the grid cell in m
    if cell_size != 0.0:
        break

Pressure_Gradient_df = Pressure_Temperature_density_df # define another data frame to put it as an argument in the function
# Now we will Calculate the pressure gradient for the multiple components we have:

Pressure_Gradient_df = compute_pressure_gradient(Pressure_Gradient_df, 'X', Background_Temperature, cell_size, np.array(Pressure_Gradient_df['n']), k)
Pressure_Gradient_df = compute_pressure_gradient(Pressure_Gradient_df, 'Y', Background_Temperature, cell_size, np.array(Pressure_Gradient_df['n']), k)
Pressure_Gradient_df = compute_pressure_gradient(Pressure_Gradient_df, 'Z', Background_Temperature, cell_size, np.array(Pressure_Gradient_df['n']), k)
Pressure_Gradient_df = Pressure_Gradient_df.sort_values(by=['Z', 'Y', 'X']) # sort them by to their original position
# Calculate the pressure gradient magnitude
Pressure_gradient_magnitude = np.sqrt(Pressure_Gradient_df['Grad_P_X'] ** 2 + Pressure_Gradient_df['Grad_P_Y'] ** 2 + Pressure_Gradient_df['Grad_P_Z'] ** 2)
# Add the pressure gradient magnitude to the data frame
Pressure_Gradient_df['Grad_pressure_magnitude'] = Pressure_gradient_magnitude


## Add a dictionary
Pressure_Gradient_data_dictionary = {
    'X': Pressure_Gradient_df['X'],
    'Y': Pressure_Gradient_df['Y'],
    'Z': Pressure_Gradient_df['Z'],
    'L': Pressure_Gradient_df['L'],
    'dP_dx': Pressure_Gradient_df['Grad_P_X'],
    'dP_dy': Pressure_Gradient_df['Grad_P_Y'],
    'dP_dz': Pressure_Gradient_df['Grad_P_Z'],
    'Gradient_pressure_magnitude': Pressure_Gradient_df['Grad_pressure_magnitude']
}
# Create the DataFrame
Pressure_gradient_data_frame = pd.DataFrame(Pressure_Gradient_data_dictionary)

# Just for the precession of python
Pressure_gradient_data_frame['X'] = Pressure_gradient_data_frame['X'].astype(np.float64) # for Python precision
Pressure_gradient_data_frame['Y'] = Pressure_gradient_data_frame['Y'].astype(np.float64) # for Python precision
Pressure_gradient_data_frame['Z'] = Pressure_gradient_data_frame['Z'].astype(np.float64) # for Python precision
Pressure_gradient_data_frame['L'] = Pressure_gradient_data_frame['L'].astype(np.float64) # for Python precision
Pressure_gradient_data_frame['dP_dx'] = Pressure_gradient_data_frame['dP_dx'].astype(np.float64) # for Python precision
Pressure_gradient_data_frame['dP_dy'] = Pressure_gradient_data_frame['dP_dy'].astype(np.float64) # for Python precision
Pressure_gradient_data_frame['dP_dz'] = Pressure_gradient_data_frame['dP_dz'].astype(np.float64) # for Python precision
Pressure_gradient_data_frame['Gradient_pressure_magnitude'] = Pressure_gradient_data_frame['Gradient_pressure_magnitude'].astype(np.float64) # for Python precision
# Saving the data frame and creating it as an Excel file
Pressure_gradient_data_frame.to_excel('cartesian_Pressure_gradient_cgs_L_' + str(Resolution) + '.xlsx', index=False,
                       engine='openpyxl')

'''_________________________________________________________________________________________________________
                                        CALCULATING THE LORENTZ FORCE
_________________________________________________________________________________________________________'''

# Loading the files
file_directory = "/Users/majd_noel/Documents/Western_University_Overgraduate_studies/Academics/Research/Python_project/ALl_files" # directory of the file
file_magnetic_field_name = "B_field_s" + str(Resolution)+ "_excel.xlsx" # the name of the file. It has to end in the extension of .xlsx
## Calling the function of loading the files
magnetic_field_df = load_file(file_directory, file_magnetic_field_name) # magnetic_field_df is a file that has the magnetic field variables in it

## Calling the function to calculate the Lorentz force
Lorentz_force_df = compute_Lorentz_force(magnetic_field_df, cell_size, Resolution)

## Add a dictionary for the forces only
Lorentz_force_data_dictionary = {
    'X': Lorentz_force_df['X'],
    'Y': Lorentz_force_df['Y'],
    'Z': Lorentz_force_df['Z'],
    'L': Lorentz_force_df['L'],
    'F_L_x': Lorentz_force_df['L_F_x_component'],
    'F_L_y': Lorentz_force_df['L_F_y_component'],
    'F_L_z': Lorentz_force_df['L_F_z_component'],
    'F_L_magnitude': Lorentz_force_df['L_F_magnitude']
}
## Creating a data frame from the dictionary
Lorentz_force_data_frame = pd.DataFrame(Lorentz_force_data_dictionary)
# Saving the data frame
Lorentz_force_data_frame.to_excel('cartesian_Lorentz_forces_cgs_L_' + str(Resolution) + '.xlsx', index=False, engine='openpyxl')

'''_________________________________________________________________________________________________________
                                        Calculating the net force
_________________________________________________________________________________________________________'''

## In this file, we will combine all the forces together and then compute the net force
All_forces_df = summation_and_net_force(Gravity_df, Pressure_gradient_data_frame, Lorentz_force_data_frame)

# Saving the data frame
All_forces_df.to_excel('ALL_and_net_forces_cgs_L_' + str(Resolution) + '.xlsx', index=False, engine='openpyxl')


'''_________________________________________________________________________________________________________
                                        Removing the boundaries
_________________________________________________________________________________________________________'''
# Since it is not very accurate to calculate the boundaries of the cube. The best thing is not to show them in the plot
# and remove them

## Calling the function:
df_update = removing_boundaries(All_forces_df, Resolution)

'''_________________________________________________________________________________________________________
                                        Plotting the vector fields
_________________________________________________________________________________________________________'''
# The dictionary for plotting

All_force_data_dictionary = {
    'X': df_update['X'],
    'Y': df_update['Y'],
    'Z': df_update['Z'],
    'L': df_update['L'],
    'F_x': df_update['F_net_x'],
    'F_y': df_update['F_net_y'],
    'F_z': df_update['F_net_z'],
}

All_force_data = pd.DataFrame(All_force_data_dictionary) # Transferring the dictionary into a data frame
# First, in order to plot the data, we need to transfer the valuables in cylindrical coordinates
All_force_data.to_excel('All_force_boundaries_Removed_cartesian_cgs_L_' + str(Resolution) + '.xlsx', index=False, engine='openpyxl')


## Calling the cartesian to cylindrical function to switch the cartesian data into cylindrical
Net_force_Cylindrical_df = cartesian_to_cylindrical(All_force_data, 'Net force', Resolution)
# Saving the data frame
Net_force_Cylindrical_df.to_excel('Net_forces_cylindrical_coordinates_cgs_L_' + str(Resolution) + '.xlsx', index=False, engine='openpyxl')

# Taking care of the boundaries. We will create a data frame to exclude all the boundary points, because it is hard
# to calculate

##In order for use to plot the function, we need take the data and average it on the angle
Averaged_data_frame = Net_force_Cylindrical_df.groupby(['r', 'z']).mean().reset_index() # Averaging the quantities over the azimuthal direction

# Saving the data frame of averaged data
Averaged_data_frame.to_excel('Net_force_averaged_cylindrical_coordinates_L_' + str(Resolution) + '.xlsx', index=False, engine='openpyxl')

# Recalling the data frame into another one
Averaged_df = Averaged_data_frame

#Plotting the function.
# To have good visualization, we need to transfer the r and z to AU, and take scale the forces by their magnitude

Plotting_variables = plotting_vector_field(Averaged_df, AU_TO_CM, Resolution)
fig1, fig2 = Plotting_variables  # Now both figures are stored

