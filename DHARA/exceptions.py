# ------------------------------------------------ Importing files -------------------------------------------------- #
from para import *
# ------------------------------------------------------------------------------------------------------------------- #

# ---------------------------------------------- Raising exceptions ------------------------------------------------- #

if device not in ['CPU', 'GPU']:
    raise ValueError('Unsupported device: {}. Please type CPU or GPU clearly!'.format(device))

if dimension not in [2,3]:
    raise ValueError('Unsupported dimension: {}. Only 2D and 3D are allowed!'.format(dimension))

if dt_far < dt:
    raise Exception('dt_far must be less than or equal to dt!')

if integration_scheme == 'simpsons13':
    if Nz%2 == 0 or Nx%2 == 0 or Ny%2 == 0:
        raise Exception('Keep the number of grids, Nx, Ny and Nz as odd numbers! (This make sure that Simpson 1/3 is computed correctly.)')

if save_global_quantities_dt < dt:
    raise Exception('Keep save_global_quantities_dt equal to or less than dt!')

# ------------------------------------------------------------------------------------------------------------------- #
