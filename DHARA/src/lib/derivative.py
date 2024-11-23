# ------------------------------------------------ Importing files -------------------------------------------------- #
import para
import DHARA.src.univ.grid as grid
# ------------------------------------------------------------------------------------------------------------------- #

# ------------------------------------- Coefficients for non-uniform grid ------------------------------------------- #

# Coefficients for 2nd order forward difference at z = 0 for non-uniform grid
h0,h1,h2=0,0,0
h0 = (2*grid.Z[0] - grid.Z[1] -  grid.Z[2])/((grid.Z[1] - grid.Z[0])*(grid.Z[2] - grid.Z[0]))
h1 = (grid.Z[2] -  grid.Z[0])/((grid.Z[1] - grid.Z[0])*(grid.Z[2] - grid.Z[1]))
h2 = (grid.Z[0] -  grid.Z[1])/((grid.Z[2] - grid.Z[1])*(grid.Z[2] - grid.Z[0]))

# Coefficients for 2nd order backward difference at z = Lz for non-uniform grid
h_0,h_1,h_2=0,0,0
h_0 = (2*grid.Z[-1] - grid.Z[-2] -  grid.Z[-3])/((grid.Z[-1] - grid.Z[-2])*(grid.Z[-1] - grid.Z[-3]))
h_1 = (grid.Z[-3] -  grid.Z[-1])/((grid.Z[-2] - grid.Z[-3])*(grid.Z[-1] - grid.Z[-2]))
h_2 = (grid.Z[-1] -  grid.Z[-2])/((grid.Z[-1] - grid.Z[-3])*(grid.Z[-2] - grid.Z[-3]))

# Coefficients for 2nd order 2nd derivative forward difference at z = 0 for non-uniform grid
g0,g1,g2,g3=0,0,0,0
g0 = 2*(grid.Z[1]+grid.Z[2]+grid.Z[3]-3*grid.Z[0])/((grid.Z[1]-grid.Z[0])*(grid.Z[2]-grid.Z[0])*(grid.Z[3]-grid.Z[0]))
g1 = 2*(grid.Z[2]+grid.Z[3]-2*grid.Z[0])/((grid.Z[1]-grid.Z[0])*(grid.Z[2]-grid.Z[1])*(grid.Z[1]-grid.Z[3]))
g2 = -2*(grid.Z[1]+grid.Z[3]-2*grid.Z[0])/((grid.Z[2]-grid.Z[0])*(grid.Z[2]-grid.Z[1])*(grid.Z[2]-grid.Z[3]))
g3 = 2*(grid.Z[1]+grid.Z[2]-2*grid.Z[0])/((grid.Z[3]-grid.Z[0])*(grid.Z[2]-grid.Z[3])*(grid.Z[3]-grid.Z[1]))

# Coefficients for 2nd order 2nd derivative backward difference at z = Lz for non-uniform grid
g_0,g_1,g_2,g_3=0,0,0,0
g_0 = -2*(grid.Z[-2]+grid.Z[-3]+grid.Z[-4]-3*grid.Z[-1])/((grid.Z[-2]-grid.Z[-1])*(grid.Z[-3]-grid.Z[-1])*(grid.Z[-4]-grid.Z[-1]))
g_1 = -2*(grid.Z[-3]+grid.Z[-4]-2*grid.Z[-1])/((grid.Z[-2]-grid.Z[-1])*(grid.Z[-3]-grid.Z[-2])*(grid.Z[-2]-grid.Z[-4]))
g_2 = 2*(grid.Z[-2]+grid.Z[-4]-2*grid.Z[-1])/((grid.Z[-3]-grid.Z[-1])*(grid.Z[-3]-grid.Z[-2])*(grid.Z[-3]-grid.Z[-4]))
g_3 = -2*(grid.Z[-2]+grid.Z[-3]-2*grid.Z[-1])/((grid.Z[-4]-grid.Z[-1])*(grid.Z[-3]-grid.Z[-4])*(grid.Z[-4]-grid.Z[-2]))

# ------------------------------------------------------------------------------------------------------------------- #

# ----------------------------------------------- First derivatives ------------------------------------------------- #

# d_type = 0 : Central Difference, d_type = 1 : Forward Difference, d_type = 2 : Backward Difference

def dfx(f,d_type):
    # Derivative of array 'f' w.r.t. x

    if d_type == 0:
        # Central difference of array 'f' w.r.t. x
        if para.boundary_u_x == 'PD' and para.boundary_T_x == 'PD':
            # For peridoic f
            grid.temp_dx[1:-1] = (f[2:] - f[:-2])/(2*grid.dx)
            # Boundary terms when f is periodic in x
            grid.temp_dx[0] = (f[1] - f[-2])/(2*grid.dx)
            grid.temp_dx[-1] = (grid.temp_dx[0])
    elif d_type == 1:
        # Forward difference of array 'f' w.r.t. x
        if para.boundary_u_x == 'PD' and para.boundary_T_x == 'PD':
            # For peridoic f
            grid.temp_dx[:-1] = (f[1:] - f[:-1])/(grid.dx)
            # Boundary terms when f is periodic in x
            grid.temp_dx[-1] = (f[1] - f[-1])/(grid.dx)
    elif d_type == 2:
        # Backward difference of array 'f' w.r.t. x
        if para.boundary_u_x == 'PD' and para.boundary_T_x == 'PD':
            # For peridoic f
            grid.temp_dx[1:] = (f[1:] - f[:-1])/(grid.dx)
            # Boundary terms when f is periodic in x
            grid.temp_dx[0] = (f[0] - f[-2])/(grid.dx)

    return grid.temp_dx


def dfy(f,d_type):
    if d_type == 0:
        # Central difference of array 'f' w.r.t. x
        if para.boundary_u_x == 'PD' and para.boundary_T_x == 'PD':
            # For peridoic f
            grid.temp_dy[:,1:-1,:] = (f[:,2:,:] - f[:,:-2,:])/(2*grid.dy)
            # Boundary terms when f is periodic in x
            grid.temp_dy[:,0,:] = (f[:,1,:] - f[:,-2,:])/(2*grid.dy)
            grid.temp_dy[:,-1,:] = (grid.temp_dy[:,0,:])
    elif d_type == 1:
        # Forward difference of array 'f' w.r.t. x
        if para.boundary_u_x == 'PD' and para.boundary_T_x == 'PD':
            # For peridoic f
            grid.temp_dy[:,:-1,:] = (f[:,1:,:] - f[:,:-1,:])/(grid.dy)
            # Boundary terms when f is periodic in x
            grid.temp_dy[:,-1,:] = (f[:,1,:] - f[:,-1,:])/(grid.dy)
    elif d_type == 2:
        # Backward difference of array 'f' w.r.t. x
        if para.boundary_u_x == 'PD' and para.boundary_T_x == 'PD':
            # For peridoic f
            grid.temp_dy[:,1:,:] = (f[:,1:,:] - f[:,:-1,:])/(grid.dy)
            # Boundary terms when f is periodic in x
            grid.temp_dy[:,0,:] = (f[:,0,:] - f[:,-2,:])/(grid.dy)
    return grid.temp_dy


def dfz(f,d_type):
    # Derivative of array 'f' w.r.t. z

    if para.dimension == 2:
        if d_type == 0: 
            # Central difference of array 'f' w.r.t. z
            grid.temp_dz[:,1:-1] = (f[:,2:] - (grid.hz[1:-1]/grid.hz[:-2])**2*f[:,:-2] - (1-(grid.hz[1:-1]/grid.hz[:-2])**2)*f[:,1:-1])/(grid.hz[1:-1]*(1+(grid.hz[1:-1]/grid.hz[:-2])))
            # Boundary terms using forward and backward differences
            grid.temp_dz[:,0] = (h0*f[:,0] + h1*f[:,1] + h2*f[:,2])
            grid.temp_dz[:,-1] = (h_0*f[:,-1] + h_1*f[:,-2] + h_2*f[:,-3])           
        elif d_type == 1:
            # Forward difference of array 'f' w.r.t. z
            grid.temp_dz[:,:-1] = (f[:,1:] - f[:,:-1])/(grid.Z[1:]-grid.Z[:-1])
            # Boundary terms using backward difference
            grid.temp_dz[:,-1] = (f[:,-1] - f[:,-2])/(grid.Z[-1]-grid.Z[-2])
        elif d_type == 2:
            # Backward difference of array 'f' w.r.t. z
            grid.temp_dz[:,1:] = (f[:,1:] - f[:,:-1])/(grid.Z[1:]-grid.Z[:-1])
            # Boundary terms using forward difference
            grid.temp_dz[:,0] = (f[:,1] - f[:,0])/(grid.Z[1]-grid.Z[0])

    elif para.dimension == 3:
        if d_type == 0: 
            # Central difference of array 'f' w.r.t. z
            grid.temp_dz[:,:,1:-1] = (f[:,:,2:] - (grid.hz[1:-1]/grid.hz[:-2])**2*f[:,:,:-2] - (1-(grid.hz[1:-1]/grid.hz[:-2])**2)*f[:,:,1:-1])/(grid.hz[1:-1]*(1+(grid.hz[1:-1]/grid.hz[:-2])))
            # Boundary terms using forward and backward differences
            grid.temp_dz[:,:,0] = (h0*f[:,:,0] + h1*f[:,:,1] + h2*f[:,:,2])
            grid.temp_dz[:,:,-1] = (h_0*f[:,:,-1] + h_1*f[:,:,-2] + h_2*f[:,:,-3])           
        elif d_type == 1:
            # Forward difference of array 'f' w.r.t. z
            grid.temp_dz[:,:,:-1] = (f[:,:,1:] - f[:,:,:-1])/(grid.Z[1:]-grid.Z[:-1])
            # Boundary terms using backward difference
            grid.temp_dz[:,:,-1] = (f[:,:,-1] - f[:,:,-2])/(grid.Z[-1]-grid.Z[-2])
        elif d_type == 2:
            # Backward difference of array 'f' w.r.t. z
            grid.temp_dz[:,:,1:] = (f[:,:,1:] - f[:,:,:-1])/(grid.Z[1:]-grid.Z[:-1])
            # Boundary terms using forward difference
            grid.temp_dz[:,:,0] = (f[:,:,1] - f[:,:,0])/(grid.Z[1]-grid.Z[0])

    return grid.temp_dz

# ------------------------------------------------------------------------------------------------------------------- #

# ------------------------------------ Second derivative w.r.t 'z' at boundaries ------------------------------------ #

def dfzz_b(f,p):
    if para.dimension == 2:
        if p == 0:
            return g0*f[:,0]+g1*f[:,1]+g2*f[:,2]+g3*f[:,3]
        elif p == -1:
            return g_0*f[:,-1]+g_1*f[:,-2]+g_2*f[:,-3]+g_3*f[:,-4]
    elif para.dimension == 3:
        if p == 0:
            return g0*f[:,:,0]+g1*f[:,:,1]+g2*f[:,:,2]+g3*f[:,:,3]
        elif p == -1:
            return g_0*f[:,:,-1]+g_1*f[:,:,-2]+g_2*f[:,:,-3]+g_3*f[:,:,-4]
    
# ------------------------------------------------------------------------------------------------------------------- #

