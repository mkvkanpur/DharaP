# ------------------------------------------------ Importing files -------------------------------------------------- #
import para
from DHARA.config.config import ncp
import DHARA.src.univ.grid as grid
# ------------------------------------------------------------------------------------------------------------------- #

# ------------------------------------------------- Integration ----------------------------------------------------- #

# Trapezoidal and Simpson's 1/3 rule for 1D, 2D and 3D integrals

if para.dimension == 2:
    def trap1d(f):
        return (ncp.sum(f) - (f[0] + f[-1])/2)*grid.dx   
    def trap2d_u(f):
        return trap1d((ncp.sum(f,axis=1) - (f[:,0] + f[:,-1])/2)*grid.dz)   
    def trap2d_nu(f):
        return trap1d(ncp.sum((f[:,1:]+f[:,:-1])*(grid.Z[1:]-grid.Z[:-1])/2,axis=1))

    def simps1d(f):
        return ncp.sum((f[0:-2:2]+4*f[1:-1:2]+f[2::2])*(grid.dx)/3)
    def simps2d_u(f):
        return simps1d(ncp.sum((f[:,0:-2:2]+4*f[:,1:-1:2]+f[:,2::2])*(grid.dz)/3,axis=1))   
    def simps2d_nu(f):
        return simps1d(ncp.sum(((2-(grid.Z[2::2]-grid.Z[1:-1:2])/(grid.Z[1:-1:2]-grid.Z[:-2:2]))*f[:,0:-2:2] \
                    + (2+(grid.Z[2::2]-grid.Z[1:-1:2])/(grid.Z[1:-1:2]-grid.Z[:-2:2])+(grid.Z[1:-1:2]-grid.Z[:-2:2])/(grid.Z[2::2]-grid.Z[1:-1:2]))*f[:,1:-1:2] \
                    + (2-(grid.Z[1:-1:2]-grid.Z[:-2:2])/(grid.Z[1:-1:2]-grid.Z[:-2:2]))*f[:,2::2])*(grid.Z[2::2]-grid.Z[:-2:2])/6,axis=1))

elif para.dimension == 3:
    def trap1d(f):
        return (ncp.sum(f) - (f[0] + f[-1])/2)*grid.dx
    def trap2d(f):
        return trap1d((ncp.sum(f,axis=1) - (f[:,0] + f[:,-1])/2)*grid.dy)
    def trap3d_u(f):
        return trap2d((ncp.sum(f,axis=2) - (f[:,:,0] + f[:,:,-1])/2)*grid.dz)
    def trap3d_nu(f):
        return trap2d(ncp.sum((f[:,:,1:]+f[:,:,:-1])*(grid.Z[1:]-grid.Z[:-1])/2,axis=2))
    
    def simps1d(f):
        return ncp.sum((f[0:-2:2]+4*f[1:-1:2]+f[2::2])*(grid.dx)/3)
    def simps2d(f):
        return simps1d(ncp.sum((f[:,0:-2:2]+4*f[:,1:-1:2]+f[:,2::2])*(grid.dy)/3,axis=1))
    def simps3d_u(f):
        return simps2d(ncp.sum((f[:,:,0:-2:2]+4*f[:,:,1:-1:2]+f[:,:,2::2])*grid.dz/3,axis=2))
    def simps3d_nu(f):
        return simps2d(ncp.sum(((2-(grid.Z[2::2]-grid.Z[1:-1:2])/(grid.Z[1:-1:2]-grid.Z[:-2:2]))*f[:,:,0:-2:2] \
                    + (2+(grid.Z[2::2]-grid.Z[1:-1:2])/(grid.Z[1:-1:2]-grid.Z[:-2:2])+(grid.Z[1:-1:2]-grid.Z[:-2:2])/(grid.Z[2::2]-grid.Z[1:-1:2]))*f[:,:,1:-1:2] \
                    + (2-(grid.Z[1:-1:2]-grid.Z[:-2:2])/(grid.Z[1:-1:2]-grid.Z[:-2:2]))*f[:,:,2::2])*(grid.Z[2::2]-grid.Z[:-2:2])/6,axis=2))


# Initialize integration functions
if para.integration_scheme == 'trapezoidal':
    integrate_1d = trap1d
    if para.dimension == 2:
        if para.grid_type == 'uniform':
            integrate_2d = trap2d_u
        elif para.grid_type == 'non-uniform':
            integrate_2d = trap2d_nu
    elif para.dimension == 3:
        integrate_2d = trap2d
        if para.grid_type == 'uniform':
            integrate_3d = trap3d_u
        elif para.grid_type == 'non-uniform':
            integrate_3d = trap3d_nu

if para.integration_scheme == 'simpsons13':
    integrate_1d = simps1d
    if para.dimension == 2:
        if para.grid_type == 'uniform':
            integrate_2d = simps2d_u
        elif para.grid_type == 'non-uniform':
            integrate_2d = simps2d_nu
    elif para.dimension == 3:
        integrate_2d = simps2d
        if para.grid_type == 'uniform':
            integrate_3d = simps3d_u
        elif para.grid_type == 'non-uniform':
            integrate_3d = simps3d_nu

# ------------------------------------------------------------------------------------------------------------------- #