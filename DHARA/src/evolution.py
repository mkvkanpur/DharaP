# ------------------------------------------------ Importing files -------------------------------------------------- #
import sys
import para
from DHARA.config.config import ncp
import DHARA.src.univ.grid as grid
from DHARA.src.lib.compressible import Compressible
from DHARA.src.lib.derivative import *
from DHARA.src.univ.updateprimitive import *
from DHARA.src.univ.data_io import print_output, saving_per_timestep
# ------------------------------------------------------------------------------------------------------------------- #

# -------------------------------------------- Time-stepping functions ---------------------------------------------- #

def CFL_condition(compress):
    # CFL condition
    if para.dimension == 2:
        dt_C = para.cfl_cons/ncp.max(ncp.abs(compress.ux)/grid.dx + ncp.abs(compress.uz)/(grid.hz))
    elif para.dimension == 3:
        dt_C = para.cfl_cons/ncp.max(ncp.abs(compress.ux)/grid.dx + ncp.abs(compress.uy)/grid.dy + ncp.abs(compress.uz)/(grid.hz))
    grid.dt_c = min(para.dt,dt_C)

    # Stop the code if dt_c becomes too small
    if grid.dt_c < 1e-10:
        print ('# Try different parameters.. Simulation blew up!')    
        sys.exit(1)


def time_advance_single_step(dt, d_type, compress):
    # Q for single time step dt
    
    compress.compute_convective_flux_x()
    compress.compute_viscous_flux_x(d_type)
    compress.flux_derivative_x(d_type)

    compress.Q -= dt*(compress.F)
    
    if para.dimension == 3:
        compress.compute_convective_flux_y()
        compress.compute_viscous_flux_y(d_type)
        compress.flux_derivative_y(d_type)

        compress.Q -= dt*(compress.F)

    compress.compute_convective_flux_z()
    compress.compute_viscous_flux_z(d_type)
    compress.flux_derivative_z(d_type)

    compress.Q -= dt*(compress.F)
        
    # Source term
    if para.dimension == 2:
        compress.Q[2] -= dt*(compress.rho*grid.C3)
    elif para.dimension == 3:
        compress.Q[3] -= dt*(compress.rho*grid.C3)
    
    return compress


def time_advance_euler(compress):
    # Time advance using Euler method, 1st order

    t = para.tinit

    if para.print_output_terminal == True:
        print_output(t,compress)               # Print output

    saving_per_timestep(t,compress)            # Saving output files

    while t <= para.tfinal+para.dt:
                
        CFL_condition(compress)

        compress.update_conserved()
        time_advance_single_step(grid.dt_c, 0, compress)
        update_primitive(compress)

        t = t+grid.dt_c

        if t >= (para.t_far):
            para.dt = para.dt_far
        
        if para.print_output_terminal == True:
            print_output(t,compress)               # Print output

        saving_per_timestep(t,compress)            # Saving output files
    

def time_advance_RK2(compress):
    # Time advance using RK2 method, 2nd order

    t = para.tinit

    if para.print_output_terminal == True:
        print_output(t,compress)               # Print output

    saving_per_timestep(t,compress)            # Saving output files

    while t <= para.tfinal+para.dt:
                
        CFL_condition(compress)

        # Saving Q of first point
        compress.update_conserved()
        compress.Q_copy1 = ncp.copy(compress.Q)

        # Finding the Q's of middle point 
        compress=time_advance_single_step(grid.dt_c/2, 0, compress)
        update_primitive(compress)

        # Finding the Q's of next point
        compress.Q = ncp.copy(compress.Q_copy1)
        compress=time_advance_single_step(grid.dt_c, 0, compress)
        update_primitive(compress)

        t = t+grid.dt_c

        if t >= (para.t_far):
            para.dt = para.dt_far
        
        if para.print_output_terminal == True:
            print_output(t,compress)               # Print output

        saving_per_timestep(t,compress)            # Saving output files
        
# ------------------------------------------------------------------------------------------------------------------- #
