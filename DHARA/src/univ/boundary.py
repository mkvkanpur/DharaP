# ------------------------------------------------ Importing files -------------------------------------------------- #
import para
import DHARA.src.univ.grid as grid
from DHARA.src.lib.compressible import Compressible
from DHARA.src.lib.derivative import *
# ------------------------------------------------------------------------------------------------------------------- #

# -------------------------------------------- Boundary conditions -------------------------------------------------- #

def imposeBC_u(compress):
    # Boundary condition for velocity fields

    if para.dimension == 2:
        if para.boundary_u_x == 'PD' and para.boundary_u_z == 'NS':
            # Periodic boundary condition for velocity in x-direction and no-slip in z-direction
            compress.ux[:,0] = compress.ux[:,-1] = 0
            compress.uz[:,0] = compress.uz[:,-1] = 0

    elif para.dimension == 3:
        if para.boundary_u_x == 'PD' and para.boundary_u_y == 'PD' and para.boundary_u_z == 'NS':
            # Periodic boundary condition for velocity in x and y-direction and no-slip in z-direction
            compress.ux[:,:,0] = compress.ux[:,:,-1] = 0
            compress.uy[:,:,0] = compress.uy[:,:,-1] = 0
            compress.uz[:,:,0] = compress.uz[:,:,-1] = 0


def imposeBC_T(compress):
    # Boundary condition for temperature

    if para.dimension == 2:
        if para.boundary_T_x == 'PD' and para.boundary_T_z == 'FX':
            # Fixed boundary on vertical walls, periodic on horizontal walls for temperature
            compress.T[:,0] = grid.T_b
            compress.T[:,-1] = grid.T_t

    elif para.dimension == 3:
        if para.boundary_T_x == 'PD' and para.boundary_T_y == 'PD' and para.boundary_T_z == 'FX':
            # Fixed boundary on vertical walls, periodic on horizontal walls for temperature
            compress.T[:,:,0] = grid.T_b
            compress.T[:,:,-1] = grid.T_t


def imposeBC_rho(compress):
    # Boundary condition for density

    if para.dimension == 2:
        if para.scheme == 'MC' or para.scheme == 'TVD-MC' and para.boundary_u_x == 'PD' and para.boundary_u_z == 'NS':
            if para.coordinates == 'cartesian':
                # Using the equation, del_z(p) = - rho g + del_z(tau_zz). Note that boundary terms are second order accurate.      
                compress.rho[:,0] = ((grid.C2/(grid.C1))*((4/3)*dfzz_b(compress.uz,0)) - grid.T_b*(h1*compress.rho[:,1] + h2*compress.rho[:,2]))/(2*h0*grid.T_b + h1*compress.T[:,1] + h2*compress.T[:,2] + grid.C3/grid.C1)
                compress.rho[:,-1] = ((grid.C2/(grid.C1))*((4/3)*dfzz_b(compress.uz,-1)) - grid.T_t*(h_1*compress.rho[:,-2] + h_2*compress.rho[:,-3]))/(2*h_0*grid.T_t + h_1*compress.T[:,-2] + h_2*compress.T[:,-3] + grid.C3/grid.C1)
            elif para.coordinates == 'polar':
                # Using the equation, del_z(p) = - rho g + del_z(tau_zz). Note that boundary terms are second order accurate.      
                compress.rho[:,0] = ((grid.C2/(grid.C1))*((4/3)*dfzz_b(compress.uz,0) - (8/3)*dfz_s(compress.uz,0)/para.radius_b) - grid.T_b*(h1*compress.rho[:,1] + h2*compress.rho[:,2]))/(2*h0*grid.T_b + h1*compress.T[:,1] + h2*compress.T[:,2] + grid.C3/grid.C1)
                compress.rho[:,-1] = ((grid.C2/(grid.C1))*((4/3)*dfzz_b(compress.uz,-1) - (8/3)*dfz_s(compress.uz,para.Nz-1)) - grid.T_t*(h_1*compress.rho[:,-2] + h_2*compress.rho[:,-3]))/(2*h_0*grid.T_t + h_1*compress.T[:,-2] + h_2*compress.T[:,-3] + grid.C3/grid.C1)
        pass
    elif para.dimension == 3:
        if para.scheme == 'MC' or para.scheme == 'TVD-MC' and para.boundary_u_x == 'PD' and para.boundary_u_y == 'PD' and para.boundary_u_z == 'NS':
            # Using the equation, del_z(p) = - rho g + del_z(tau_zz). Note that boundary terms are second order accurate.      
            compress.rho[:,:,0] = ((grid.C2/(grid.C1))*((4/3)*dfzz_b(compress.uz,0)) - grid.T_b*(h1*compress.rho[:,:,1] + h2*compress.rho[:,:,2]))/(2*h0*grid.T_b + h1*compress.T[:,:,1] + h2*compress.T[:,:,2] + grid.C3/grid.C1)
            compress.rho[:,:,-1] = ((grid.C2/(grid.C1))*((4/3)*dfzz_b(compress.uz,-1)) - grid.T_t*(h_1*compress.rho[:,:,-2] + h_2*compress.rho[:,:,-3]))/(2*h_0*grid.T_t + h_1*compress.T[:,:,-2] + h_2*compress.T[:,:,-3] + grid.C3/grid.C1)
        pass

def boundary(compress):
    # Function to impose boundary conditions
    imposeBC_rho(compress)
    imposeBC_u(compress)
    imposeBC_T(compress)

# ------------------------------------------------------------------------------------------------------------------- #
