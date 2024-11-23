# --------------------------------------------------- Libraries ----------------------------------------------------- #
import h5py
import math as math
import os
import sys
import shutil
# ------------------------------------------------------------------------------------------------------------------- #

# ------------------------------------------------ Importing files -------------------------------------------------- #
import para
from DHARA.config.config import ncp
import DHARA.src.univ.grid as grid
from DHARA.src.lib.compressible import Compressible
from DHARA.src.univ.compute.glob import *
# ------------------------------------------------------------------------------------------------------------------- #

# ----------------------------------------------- Generate directories ---------------------------------------------- #

def gen_path():
    # Creating the output directories if they do not exist
    if not os.path.exists(para.output_dir):
        os.makedirs(para.output_dir)
    if not os.path.exists(para.output_dir + '/fields/'):
        os.makedirs(para.output_dir + '/fields/')

    # Copy para.py file
    counter = 0
    filename = "parameters_{}.py"
    while os.path.isfile(para.output_dir + filename.format(counter)):
        counter += 1
    filename = filename.format(counter)
    shutil.copy2('para.py', para.output_dir + filename)

# ------------------------------------------------------------------------------------------------------------------- #

# ------------------------------------------------ Printing functions ----------------------------------------------- #

def print_output(t, compress):
    # Print output info at the start of the simulation
    if t == para.tinit:
        print('\n# ------------------------------------- Time instance ------------------------------------- #') 
        print('\n# The following columns contains: time step (dt), time (t), Root Mean Square Velocity (Vrms)')

    # Print output at every time step
    print(grid.dt_c, t, V_rms(compress))

# ------------------------------------------------------------------------------------------------------------------- #

# -------------------------------------------- Computing global quantities ------------------------------------------ #

def global_quantites(j, t, compress):

    # Save global quantities at every time step, t_g
    grid.global_quan[0,j] = grid.dt_c
    grid.global_quan[1,j] = t
    grid.global_quan[2,j] = total_mass(compress)
    grid.global_quan[3,j] = total_potential_energy(compress)
    grid.global_quan[4,j] = total_kinetic_energy(compress)
    grid.global_quan[5,j] = total_internal_energy(compress)
    grid.global_quan[6,j] = grid.global_quan[3,j] + grid.global_quan[4,j] + grid.global_quan[5,j]
    grid.global_quan[7,j] = V_rms(compress)
    grid.global_quan[8,j] = Max_Mach(compress)
    grid.global_quan[9,j] = ncp.sqrt(para.Ra/para.Pr)*grid.global_quan[7,j]
    grid.global_quan[10,j] = Nusselt_number(compress)

    # Stop the simulation if kinetic energy is NaN
    if math.isnan(grid.global_quan[4,j]) == True:
        sys.exit('Error: Kinetic energy is NaN at time t = %f' %t)

# ------------------------------------------------------------------------------------------------------------------- #

# ---------------------------------------------- Saving functions --------------------------------------------------- #

def parameters_save():
    # Creating the parameters file, if it already exists, then it will create a new file with a number at the end
    counter = 0
    filename = "parameters_{}.txt"
    while os.path.isfile(para.output_dir + filename.format(counter)):
        counter += 1
    filename = filename.format(counter)

    with open(para.output_dir + filename, "w") as file:
        file.write('# --------------------------------------------------------- Simulation details -------------------------------------------------------- #') 
        file.write('\n\n# This run belongs to the sequential version of DHARA.')
        file.write('\n\n# Scheme details (dimension, grid_type, profile, scheme):\n')
        file.write(f'{para.dimension} {para.grid_type} {para.profile} 0 0 0 0 0')
        file.write('\n\n# Boundary conditions (boundary_u_x, boundary_u_y, boundary_u_z, boundary_T_x, boundary_T_y, boundary_T_z):\n')
        file.write(f'{para.boundary_u_x} {para.boundary_u_y} {para.boundary_u_z} {para.boundary_T_x} {para.boundary_T_y} {para.boundary_T_z} 0 0')
        file.write('\n\n# Grid parameters (Lx, Ly, Lz, Nx, Ny, Nz, beta):\n')
        file.write(f'{grid.Lx} {grid.Ly} {para.Lz} {para.Nx} {para.Ny} {para.Nz} {para.beta} 0')
        file.write('\n\n# Time-stepping parameters (dt, tinit, tfinal, dt_far, t_far, cfl_cons):\n')
        file.write(f'{para.dt} {para.tinit} {para.tfinal} {para.dt_far} {para.t_far} {para.cfl_cons} 0 0')
        file.write('\n\n# Control parameters (A, B, epsilon, D, gamma, Pr, Ra):\n')
        file.write(f'{para.A} {para.B} {para.epsilon} {para.D} {para.gamma} {para.Pr} {"{:e}".format(para.Ra)} 0')
        file.write('\n\n# Constants in governing equations (C1, C2, C3, C4):\n')
        file.write(f'{grid.C1} {grid.C2} {grid.C3} {grid.C4} 0 0 0 0')
        file.write('\n\n# Related parameters (alpha, m, dz_min, t_nu, t_kappa, Chi_0, Chi_a, Ma):\n')
        file.write(f'{grid.alpha} {grid.m} {grid.dz_min} {grid.t_nu} {grid.t_kappa} {grid.Chi_0} {grid.Chi_a} {grid.Ma}')
        file.write('\n\n# ------------------------------------------------------------------------------------------------------------------------------------- #')   
    file.close()


def global_quan_save(compress):
    # Saving the global quantities 
    filename = 'global_%.2f.h5' %(para.tinit)
    hf2 = h5py.File(para.output_dir + filename, 'w')
    if para.device == 'CPU':
        hf2.create_dataset('dt_c', data=grid.global_quan[0])
        hf2.create_dataset('t', data=grid.global_quan[1])
        hf2.create_dataset('M_T', data=grid.global_quan[2])
        hf2.create_dataset('Ue', data=grid.global_quan[3])
        hf2.create_dataset('Ke', data=grid.global_quan[4])
        hf2.create_dataset('Ie', data=grid.global_quan[5])
        hf2.create_dataset('E_T', data=grid.global_quan[6])
        hf2.create_dataset('V_rms', data=grid.global_quan[7])
        hf2.create_dataset('Ma', data=grid.global_quan[8])
        hf2.create_dataset('Re', data=grid.global_quan[9])
        hf2.create_dataset('Nu', data=grid.global_quan[10])
    elif para.device == 'GPU':
        hf2.create_dataset('dt_c', data=grid.global_quan[0].get())
        hf2.create_dataset('t', data=grid.global_quan[1].get())
        hf2.create_dataset('M_T', data=grid.global_quan[2].get())
        hf2.create_dataset('Ue', data=grid.global_quan[3].get())
        hf2.create_dataset('Ke', data=grid.global_quan[4].get())
        hf2.create_dataset('Ie', data=grid.global_quan[5].get())
        hf2.create_dataset('E_T', data=grid.global_quan[6].get())
        hf2.create_dataset('V_rms', data=grid.global_quan[7].get())
        hf2.create_dataset('Ma', data=grid.global_quan[8].get())
        hf2.create_dataset('Re', data=grid.global_quan[9].get())
        hf2.create_dataset('Nu', data=grid.global_quan[10].get())
    hf2.close()


def save_fields(t,compress):
    # Saving the fields
    field_path = para.output_dir + '/fields/'
    grid.parameters[0] = t                          # Saving the exact time of field file
    if para.dimension == 2:
        hf1 = h5py.File(field_path + "2D_%.2f.h5" % (t), 'w')
        if para.device == 'CPU':
            hf1.create_dataset('parameters', data=grid.parameters)
            hf1.create_dataset('rho', data=(compress.rho))
            hf1.create_dataset('theta', data=(compress.T - (1-para.D*((grid.Z_mesh-grid.r_b)/(1-grid.r_b)))))
            hf1.create_dataset('ux', data=compress.ux)
            hf1.create_dataset('uz', data=compress.uz)
            hf1.create_dataset('x', data=grid.X)
            hf1.create_dataset('z', data=grid.Z)
        elif para.device == 'GPU':
            hf1.create_dataset('parameters', data=grid.parameters.get())
            hf1.create_dataset('rho', data=(compress.rho).get())
            hf1.create_dataset('theta', data=(compress.T - (1-para.D*((grid.Z_mesh-grid.r_b)/(1-grid.r_b)))).get())
            hf1.create_dataset('ux', data=compress.ux.get())
            hf1.create_dataset('uz', data=compress.uz.get())
            hf1.create_dataset('x', data=grid.X.get())
            hf1.create_dataset('z', data=grid.Z.get())
        hf1.close()
    elif para.dimension == 3:
        hf1 = h5py.File(field_path + "3D_%.2f.h5" % (t), 'w')
        if para.device == 'CPU':
            hf1.create_dataset('parameters', data=grid.parameters)
            hf1.create_dataset('rho', data=(compress.rho))
            hf1.create_dataset('theta', data=(compress.T - (1-para.D*grid.Z_mesh)))
            hf1.create_dataset('ux', data=compress.ux)
            hf1.create_dataset('uy', data=compress.uy)
            hf1.create_dataset('uz', data=compress.uz)
            hf1.create_dataset('x', data=grid.X)
            hf1.create_dataset('y', data=grid.Y)
            hf1.create_dataset('z', data=grid.Z)
        elif para.device == 'GPU':
            hf1.create_dataset('parameters', data=grid.parameters.get())
            hf1.create_dataset('rho', data=(compress.rho).get())
            hf1.create_dataset('theta', data=(compress.T - (1-para.D*grid.Z_mesh)).get())
            hf1.create_dataset('ux', data=compress.ux.get())
            hf1.create_dataset('uy', data=compress.uy.get())
            hf1.create_dataset('uz', data=compress.uz.get())
            hf1.create_dataset('x', data=grid.X.get())
            hf1.create_dataset('y', data=grid.Y.get())
            hf1.create_dataset('z', data=grid.Z.get())
        hf1.close()


def saving_per_timestep(t,compress):
    if para.save_global_quantities == True and grid.s_j <= grid.Ng and ((grid.t_g_step[grid.s_j]-t)/grid.dt_c) <= grid.dt_c:
        global_quantites(grid.s_j,t,compress)    # Computing the global quantities 
        grid.s_j=grid.s_j+1
    if para.save_fields == True and grid.s_l <= grid.Nf and ((grid.t_f_step[grid.s_l]-t)/grid.dt_c) <= grid.dt_c:
        save_fields(t,compress)                  # Saving the fields
        grid.s_l=grid.s_l+1
    if grid.s_m <= grid.Ngt and ((grid.t_gt_step[grid.s_m]-t)/grid.dt_c) <= grid.dt_c:
        if para.save_global_quantities == True:
            global_quan_save(compress)           # Saving global quantities
        grid.s_m=grid.s_m+1

# ------------------------------------------------------------------------------------------------------------------- #
