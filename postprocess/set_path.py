# --------------------------------------------------- Libraries ----------------------------------------------------- #
from pathlib import Path
# ------------------------------------------------------------------------------------------------------------------- #

# ------------------------------------------------- Defining paths -------------------------------------------------- #

main_path = Path(__file__).parents[1].absolute()
DHARA_path = str(main_path)                                             # Main program path        

# Provide a custom path of the program if it is not in the same folder as the postprocessing scripts
# DHARA_path = '/home/phyguest/harshit/DHARA_V_0.2/'

# Output folder path
output_folder_path = '/mnt/disk3_4TB/harshit/dhara_seq/2druns/Ra1e5/'                                  

parameters_filename = 'parameters_0.py'                                 # Output file name
global_filename = 'global_0.00.h5'                                      # Output file name

# ------------------------------------------------------------------------------------------------------------------- #

# ------------------------------------------------- Printing paths -------------------------------------------------- #

print('\n# ---------------------------- Check the output path and filenames ---------------------------- #')
print('\nDHARA program path: ', DHARA_path)
print('Output folder path: ', output_folder_path)
print('\nParameters filename : ', parameters_filename)
print('Global filename : ', global_filename)
print('\n# --------------------------------------------------------------------------------------------- #')

# ------------------------------------------------------------------------------------------------------------------- #
