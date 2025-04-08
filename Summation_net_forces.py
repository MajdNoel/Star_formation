## In this file, the summation of all the forces will be calculated

import numpy as np
import pandas as pd

def summation_and_net_force(Gravity_file, Pressure_file, Lorentz_force):
    df = Gravity_file # making the gravity file as dataframe (we will add the pressure and Lorentz force)
    df["dP_dx"] = Pressure_file["dP_dx"]
    df["dP_dy"] = Pressure_file["dP_dy"]
    df["dP_dz"] = Pressure_file["dP_dz"]
    df["Grad_pressure_magnitude"] = Pressure_file['Gradient_pressure_magnitude']
    df["F_L_x"] = Lorentz_force['F_L_x']
    df["F_L_y"] = Lorentz_force['F_L_y']
    df["F_L_z"] = Lorentz_force['F_L_z']
    df["F_L_magnitude"] = Lorentz_force['F_L_magnitude']

    # Calculating the components of the net force
    df["F_net_x"] =  df["F_g_x"] + df["dP_dx"] + df["F_L_x"]
    df["F_net_y"] = df["F_g_y"] + df["dP_dy"] + df["F_L_y"]
    df["F_net_z"] = df["F_g_z"] + df["dP_dz"] + df["F_L_z"]

    ## In this for loop, which ever force is strongest, it will be written in the same row. This is to give us an estimate
    # when is each force strongest
    df_list = []
    for i in range(len(df["F_L_magnitude"])):
        if df["F_g_magnitude"][i] > df["Grad_pressure_magnitude"][i] and df["F_g_magnitude"][i] > df["F_L_magnitude"][i]:
            df_list.append("Gravity")
        elif df["Grad_pressure_magnitude"][i] > df["F_g_magnitude"][i] and df["Grad_pressure_magnitude"][i] > df["F_L_magnitude"][i]:
            df_list.append("Pressure Gradient")
        else :
            df_list.append("Lorentz force")

    # Adding the list to the data frame:
    df["Strongest_force"] = df_list


    # Calculating the net force of the magnitude
    df["F_net_magnitude"] = np.sqrt((df["F_net_x"]**2)+(df["F_net_y"]**2)+(df["F_net_z"]**2))

    return df



