# --------------------------------------------------- Libraries ----------------------------------------------------- #
import os
import sys
import time 
# ------------------------------------------------------------------------------------------------------------------- #

# ------------------------------------------------ Importing files -------------------------------------------------- #
import para
import DHARA.exceptions
from DHARA.config.config import ncp
from DHARA.src.lib.compressible import Compressible
from DHARA.src.evolution import time_advance_euler, time_advance_RK2
from DHARA.src.univ.boundary import boundary
from DHARA.src.univ.data_io import gen_path, parameters_save
from DHARA.init_cond.init import init_fields
# ------------------------------------------------------------------------------------------------------------------- #

# --------------------------------------------- Output files and print ---------------------------------------------- #

print('\n# Starting the simulation... DO NOT CANCEL THE PROCESS!')
if para.dimension == 2:
    print('# Simulation setup: 2D', ', Grid size: ', para.Nx, 'x', para.Nz)
elif para.dimension == 3:
    print('# Simulation setup: 3D', ', Grid size: ', para.Nx, 'x', para.Ny, 'x', para.Nz)

parent_directory = os.path.join(os.getcwd())
output_directory = os.path.join(parent_directory, para.output_dir)
print('\nPlease check the output directory for the results.\nOutput directory:', output_directory)

gen_path()                          # Creating output directories
parameters_save()                   # Saving parameters file

# ------------------------------------------------------------------------------------------------------------------- #

# ----------------------------------------------- Main simulation --------------------------------------------------- #

compress = Compressible()           # Making of Compressible object
compress.set_arrays()               # Setting of all arrays
init_fields(compress)               # Initiating initial field profiles at initial time
boundary(compress)                  # Boundary condition for initial field profiles                

ti = time.time()                    # Initial time

if para.scheme == 'EULER':
    time_advance_euler(compress)    # Time advance steps using Euler method

elif para.scheme == 'RK2':
    time_advance_RK2(compress)      # Time advance steps using RK2 method

tf = time.time()                    # Final time

# ------------------------------------------------------------------------------------------------------------------- #

print('\n# ----------------------------------------------------------------------------------------- #')
print('\n# Total time of simulation =', tf-ti, 'seconds.')        # Time taken to run the code 
print('\n# Finished the simulation. Please check the output folder for the results.')
