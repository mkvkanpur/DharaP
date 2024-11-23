# --------------------------------------------------- Libraries ----------------------------------------------------- #
import numpy as np
import sys
import h5py
import pandas as pd
# ------------------------------------------------------------------------------------------------------------------- #

# ------------------------------------------------ Importing files -------------------------------------------------- #
from set_path import *
# ------------------------------------------------------------------------------------------------------------------- #

# ----------------------------------------- Reading and defining parameters ----------------------------------------- #

with open(output_folder_path + parameters_filename) as file:
    exec(file.read())

r_b = 0                                             # Making sure that r_b = z(0) = 0 for Cartesian coordinates
alpha = 1/(gamma - 1)                               # Adiabatic index
m = (D*(1+alpha)/(D+epsilon)) - 1                   # Polytropic index

C1 = (1-r_b)/(epsilon*D*(alpha+1))                  # Normalized gas constant, p = C1*rho*T
C2 = (1-r_b)**(3/2)*np.sqrt(Pr/Ra)                  # Normalized dynamic viscoisty
C3 = 1/epsilon                                      # Normalized acceleration due to gravity
C4 = (1-r_b)**(5/2)/(epsilon*D*np.sqrt(Ra*Pr))      # Normalized thermal diffusivity

if dimension == 2:
    Lx = A*Lz
    Ly = 0
    Volume = Lx*Lz                              # Volume of the box
elif dimension == 3:
    Lx = A*Lz
    Ly = B*Lz
    Volume = Lx*Ly*Lz                           # Volume of the box
    
# ------------------------------------------------------------------------------------------------------------------- #

# -------------------------------------- Reading and defining global quantities ------------------------------------- #

with h5py.File(output_folder_path + global_filename, "r") as f:
    dt_c = np.array(f['dt_c'])
    t = np.array(f['t'])
    M_T = np.array(f['M_T'])
    Ue = np.array(f['Ue'])
    Ke = np.array(f['Ke'])
    Ie = np.array(f['Ie'])
    E_T = np.array(f['E_T'])
    Ma = np.array(f['Ma'])
    V_rms = np.array(f['V_rms'])
    Re = np.array(f['Re'])
    Nu = np.array(f['Nu'])
f.close()

M_0 = M_T[0]           # Initial total mass
E_0 = E_T[0]           # Initial energy

# ------------------------------------------------------------------------------------------------------------------- #

# ------------------------------------------- Moving average of an array -------------------------------------------- #

def moving_average(data, window_size):
    return np.convolve(data, np.ones(window_size)/window_size, mode='same')

# ------------------------------------------------------------------------------------------------------------------- #
