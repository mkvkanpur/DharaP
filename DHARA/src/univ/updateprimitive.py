# ------------------------------------------------ Importing files -------------------------------------------------- #
import para
from DHARA.config.config import ncp
import DHARA.src.univ.grid as grid
from DHARA.src.lib.compressible import Compressible
from DHARA.src.univ.boundary import *
# ------------------------------------------------------------------------------------------------------------------- #

# ---------------------------- Update primitive variables with proper boundary conditions --------------------------- #

def update_primitive(compress):
    # Update the primitive variables using Q's

    compress.rho = ncp.copy(compress.Q[0])
    imposeBC_rho(compress)
    
    compress.ux = compress.Q[1]/compress.rho

    if para.dimension == 2:
        compress.uz = compress.Q[2]/compress.rho
        imposeBC_u(compress)

        compress.T = ((compress.Q[3]/compress.rho) - ((1/2)*(compress.ux**2 + compress.uz**2)) - grid.C3*grid.Z_mesh)*(para.gamma-1)/grid.C1
        imposeBC_T(compress)
    elif para.dimension == 3:
        compress.uy = compress.Q[2]/compress.rho
        compress.uz = compress.Q[3]/compress.rho
        imposeBC_u(compress)

        compress.T = ((compress.Q[4]/compress.rho) - ((1/2)*(compress.ux**2 + compress.uy**2 + compress.uz**2)) - grid.C3*grid.Z_mesh)*(para.gamma-1)/grid.C1
        imposeBC_T(compress)

# ------------------------------------------------------------------------------------------------------------------- #
