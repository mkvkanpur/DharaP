# --------------------------------------------------- Libraries ----------------------------------------------------- #
import os
# ------------------------------------------------------------------------------------------------------------------- #

# ------------------------------------------------ Importing files -------------------------------------------------- #
from set_path import *
# ------------------------------------------------------------------------------------------------------------------- #


# --------------------------------- Creating the plots folder if it does not exist ---------------------------------- #

if not os.path.exists(output_folder_path + 'plots'):
    os.makedirs(output_folder_path + 'plots')
if not os.path.exists(output_folder_path + 'plots/density_plots'):
    os.makedirs(output_folder_path + 'plots/density_plots')

# ------------------------------------------------------------------------------------------------------------------- #

# -------------------------------------------------- Types of plot -------------------------------------------------- #

plot_type = 'Nusselt'
t_save = 100

# For fields: field
# For energies and fluxes: energies, Nusselt
# For other quantities: error_total_mass, V_rms, Mach_number, criteria, ratio 

# ------------------------------------------------------------------------------------------------------------------- #

# ----------------------------------------------- Importing plot files ---------------------------------------------- #

import plot_scripts.field_plots
import plot_scripts.output_plots 

# ------------------------------------------------------------------------------------------------------------------- #

# ------------------------------------------------ Printing output -------------------------------------------------- #

def print_output():
    print('\n# --------------------------------------------------- Making plots ---------------------------------------------------- #')
    print('\nPlot type: ', plot_type)
    print('\nPlease check the output plot folder for plots!')
    print('Output plot folder path: ', output_folder_path + 'plots')
    print('\n# --------------------------------------------------------------------------------------------------------------------- #')

# ------------------------------------------------------------------------------------------------------------------- #

