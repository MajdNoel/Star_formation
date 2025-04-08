## This function is plot the data in cylindrical coordinates as r-z planes

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
def plotting_vector_field(Averaged_data_frame, AU_TO_CM, Resolution_L, force_type):
    # Extract arrays for plotting. Transferring the padnas data frame to numpy arrays
    R_list = Averaged_data_frame['r'].values
    Z_list = Averaged_data_frame['z'].values
    L_list = Averaged_data_frame['L'].values


    lists_appending = []
    # The order of the data has to be F_r -> F_phi -> F_z as the last coloumn


    for col in Averaged_data_frame.columns:
        if col not in ['r', 'phi', 'z', 'L']:  # Exclude these columns
            lists_appending.append(np.array(Averaged_data_frame[col]))  # Append the remaining columns
    f_r_list_to_append_perm = lists_appending[0]
    f_phi_list_to_append_perm = lists_appending[1]
    f_z_list_to_append_perm = lists_appending[2]
    R_AU = R_list / AU_TO_CM  # Converting the radial coordinate form cm to AU
    Z_AU = Z_list / AU_TO_CM  # Converting the z coordinate from cm to AU
    ## Creating a dictionary to plot the important data
    data_cylind_N_AU = {
        'F_r': f_r_list_to_append_perm,
        'F_z': f_z_list_to_append_perm,
        'F_phi': f_phi_list_to_append_perm,
        'r': R_AU,
        'z': Z_AU, }

    # Creating and saving a file with
    f_r_cylind_in_cgs_AU = pd.DataFrame(data_cylind_N_AU)
    f_r_cylind_in_cgs_AU.to_csv('Net_force_AU_cgs_L_' + str(Resolution_L) + '.txt', index=False,
                                     sep='\t')
    ## Transferring the pandas data frame to numpy arrays

    R_AU = f_r_cylind_in_cgs_AU["r"].values  # Radial positions in AU
    Z_AU = f_r_cylind_in_cgs_AU["z"].values  # Z positions in AU
    function_r = f_r_cylind_in_cgs_AU["F_r"].values  # Radial pressure_gradient in N
    function_z = f_r_cylind_in_cgs_AU["F_z"].values  # Vertical pressure_gradient in N
    function_phi = f_r_cylind_in_cgs_AU["F_phi"].values

    # Computing the pressure_Gradient magnitudes to scale the vectors to one size

    function_magnitude = np.sqrt(function_r ** 2 + function_z ** 2 + function_phi ** 2)
    df_r_norm = function_r / function_magnitude  # Dividing the pressure in the radial component by the magnitude
    df_z_norm = function_z / function_magnitude  # Dividing the pressure in the z component by the magnitude

    # Using a logarithmic scale for color mapping
    function_magnitude_log = np.log10(function_magnitude)

    # Set a color range
    vmin = np.min(function_magnitude_log)  # Minimum log pressure
    vmax = np.max(function_magnitude_log)  # Maximum log pressure

    # Creating the figure
    plt.figure(figsize=(8, 6))

    # Plot the quiver plot with colors based on log pressure_gradient magnitude
    quiver = plt.quiver(R_AU, Z_AU, df_r_norm, df_z_norm, function_magnitude_log,
                        cmap='plasma', scale=50, width=0.002, headwidth=3)

    # Add color bar
    cbar = plt.colorbar(quiver)
    cbar.set_label(r'$\log_{10}$ Net force field (Dyne/cm^3)')
    cbar.mappable.set_clim(vmin, vmax)
    ###### Plotting the first figure
    fig_1 = plt.figure(1)

    # Labels and title
    plt.xlabel(r'$r$ (AU)')
    plt.ylabel(r'$z$ (AU)')
    plt.title(force_type + ' (Log Color Scale)_L_' + str(Resolution_L))
    plt.show()    # Show the plot
    ###### Plotting the next figure
    fig_2 = plt.figure(2)

    step = 7 # This is to change how crowded you want the vector to be

    # Plot fewer vectors by slicing arrays with step size
    quiver = plt.quiver(R_AU[::step], Z_AU[::step],
                        df_r_norm[::step], df_z_norm[::step],
                        function_magnitude_log[::step], cmap='plasma',
                        scale=50, width=0.002, headwidth=3)
    cbar = plt.colorbar(quiver)
    cbar.set_label(r'$\log_{10}$ Net force field (Dyne/cm^3)')
    plt.xlabel(r'$r$ (AU)')
    plt.ylabel(r'$z$ (AU)')
    plt.title(force_type + ' (Log Color Scale)_L_' + str(Resolution_L))
    plt.show() # show the plot
    return fig_1, fig_2
