# --------------------------------------------------- Libraries ----------------------------------------------------- #
import h5py
import sys
# ------------------------------------------------------------------------------------------------------------------- #

# ------------------------------------------------ Importing files -------------------------------------------------- #
import para
from DHARA.config.config import ncp
import DHARA.src.univ.grid as grid
from DHARA.src.lib.compressible import Compressible
# ------------------------------------------------------------------------------------------------------------------- #

# -------------------------------------------- Initial field profiles  ---------------------------------------------- #

def init_fields(compress):
    # Initial field profiles for density, temperature, velocity in x and z directions
    
    if para.profile == 'static':
        # Static equilibrium profiles for temperature and density 
        compress.T = 1-(para.epsilon+para.D)*((grid.Z_mesh-grid.r_b)/(1-grid.r_b))
        compress.rho = compress.T**grid.m

        if para.dimension == 2:
            # Initial velocity profiles for 2D simulations
            if para.boundary_u_x == 'PD' or para.boundary_T_x == 'PD':
                compress.ux = 1e-3*ncp.sin(para.A*grid.X_mesh*2*ncp.pi/grid.Lx)*ncp.sin(((grid.Z_mesh-grid.r_b)/(1-grid.r_b))*2*ncp.pi)
                compress.uz = - 1e-3*ncp.cos(para.A*grid.X_mesh*2*ncp.pi/grid.Lx)*ncp.sin(((grid.Z_mesh-grid.r_b)/(1-grid.r_b))*ncp.pi)

        elif para.dimension == 3:
            # Initial velocity profiles for 3D simulations
            if para.boundary_u_x == 'PD' or para.boundary_T_x == 'PD':
                compress.ux = 1e-3*ncp.sin((1/para.A)*grid.X_mesh*2*ncp.pi)*ncp.sin((1/para.B)*grid.Y_mesh*2*ncp.pi)*ncp.sin(grid.Z_mesh*2*ncp.pi) 
                compress.uy = 1e-3*ncp.sin((1/para.A)*grid.X_mesh*2*ncp.pi)*ncp.sin((1/para.B)*grid.Y_mesh*2*ncp.pi)*ncp.sin(grid.Z_mesh*2*ncp.pi) 
                compress.uz = - 1e-3*ncp.cos((1/para.A)*grid.X_mesh*2*ncp.pi)*ncp.cos((1/para.B)*grid.Y_mesh*2*ncp.pi)*ncp.sin(grid.Z_mesh*ncp.pi)

    elif para.profile == 'continue':
        # Continue previous run from some time, put it in para.tinit 
        if para.dimension == 2:
            # For 2D simulations
            with h5py.File(para.output_dir + "fields/2D_%.2f.h5" %(grid.tinit), "r") as f:
                # Get the data for the fields
                compress.rho = ncp.array(f['rho'])
                compress.T = ncp.array(f['theta']) + (1-para.D*((grid.Z_mesh-grid.r_b)/(1-grid.r_b)))
                compress.ux = ncp.array(f['ux'])
                compress.uz = ncp.array(f['uz'])
        
        elif para.dimension == 3:
            # For 3D simulations
            with h5py.File(para.output_dir + "fields/3D_%.2f.h5" %(grid.tinit), "r") as f:
                # Get the data
                compress.rho = ncp.array(f['rho'])
                compress.T = ncp.array(f['theta']) + (1-para.D*((grid.Z_mesh-grid.r_b)/(1-grid.r_b)))
                compress.ux = ncp.array(f['ux'])
                compress.uy = ncp.array(f['uy'])
                compress.uz = ncp.array(f['uz'])

    elif para.profile == 'custom':
        # Continue previous run from some time, put it in para.tinit 
        if para.dimension == 2:
            # For 2D simulations
            with h5py.File('input/2D_init.h5', "r") as f:
                # Get the data for the fields
                compress.rho = ncp.array(f['rho'])
                compress.T = ncp.array(f['theta']) + (1-para.D*((grid.Z_mesh-grid.r_b)/(1-grid.r_b)))
                compress.ux = ncp.array(f['ux'])
                compress.uz = ncp.array(f['uz'])
        
        elif para.dimension == 3:
            # For 3D simulations
            with h5py.File("input/3D_init.h5", "r") as f:
                # Get the data
                compress.rho = ncp.array(f['rho'])
                compress.T = ncp.array(f['theta']) + (1-para.D*((grid.Z_mesh-grid.r_b)/(1-grid.r_b)))
                compress.ux = ncp.array(f['ux'])
                compress.uy = ncp.array(f['uy'])
                compress.uz = ncp.array(f['uz'])    

# ------------------------------------------------------------------------------------------------------------------- #

