# ------------------------------------------------ Importing files -------------------------------------------------- #
import para
from DHARA.config.config import ncp
import DHARA.src.univ.grid as grid
from DHARA.src.lib.compressible import Compressible
from DHARA.src.lib.derivative import *
from DHARA.src.univ.integrate import *
# ------------------------------------------------------------------------------------------------------------------- #

# ------------------------------------------- Computing global quantites -------------------------------------------- #

def total_mass(compress):
    # Total mass
    if para.dimension == 2:
        return integrate_2d(compress.rho)
    elif para.dimension == 3:
        return integrate_3d(compress.rho)

def total_potential_energy(compress):
    # Total gravitational potential energy integral
    if para.dimension == 2:
        return grid.C3*integrate_2d((compress.rho*grid.Z_mesh))
    elif para.dimension == 3:
        return grid.C3*integrate_3d((compress.rho*grid.Z_mesh))
        
def total_kinetic_energy(compress):
    # Total kinetic energy integral
    if para.dimension == 2:
        return integrate_2d((compress.rho)*(1/2)*(compress.ux**2 + compress.uz**2))
    elif para.dimension == 3:
        return integrate_3d((compress.rho)*(1/2)*(compress.ux**2 + compress.uy**2 + compress.uz**2))

def total_internal_energy(compress):
    # Total kinetic energy integral
    if para.dimension == 2:
        return (grid.C1/(para.gamma-1))*integrate_2d((compress.rho)*(compress.T))
    elif para.dimension == 3:
        return (grid.C1/(para.gamma-1))*integrate_3d((compress.rho)*(compress.T))

def V_rms(compress):
    # Root mean square velocity
    if para.dimension == 2:
        return ncp.sqrt(integrate_2d((compress.ux**2 + compress.uz**2))/grid.Volume)
    elif para.dimension == 3:
        return ncp.sqrt(integrate_3d((compress.ux**2 + compress.uy**2+ compress.uz**2))/grid.Volume)

def Max_Mach(compress):
    # Maximum Mach number
    if para.dimension == 2:
        return ncp.max(ncp.sqrt((compress.ux**2+compress.uz**2)/(para.gamma*grid.C1*compress.T)))
    elif para.dimension == 3:
        return ncp.max(ncp.sqrt((compress.ux**2+compress.uy**2+compress.uz**2)/(para.gamma*grid.C1*compress.T)))
    
def Nusselt_number(compress):
    # Nusselt number at boundaries
    compress.temp = compress.T - (1-para.D*grid.Z_mesh)
    if para.dimension == 2:    
        return - integrate_1d((h0*compress.temp[:,0] + h1*compress.temp[:,1] + h2*compress.temp[:,2] + h_0*compress.temp[:,-1] + h_1*compress.temp[:,-2] + h_2*compress.temp[:,-3])/2)/(para.epsilon*grid.Area_flux[0])
    elif para.dimension == 3:
        return - integrate_2d((h0*compress.temp[:,:,0] + h1*compress.temp[:,:,1] + h2*compress.temp[:,:,2] + h_0*compress.temp[:,:,-1] + h_1*compress.temp[:,:,-2] + h_2*compress.temp[:,:,-3])/2)/(grid.Area_flux*para.epsilon)
    
# ------------------------------------------------------------------------------------------------------------------- #