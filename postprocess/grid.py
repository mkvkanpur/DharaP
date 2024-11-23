# ------------------------------------------------ Importing files -------------------------------------------------- #
from read_output import * 
# ------------------------------------------------------------------------------------------------------------------- #

# ------------------------------------------------- Grid variables -------------------------------------------------- #

# For Cartesian coordinates   
dx = Lx/(Nx-1)                                  # Length between two consecutive grid points in x-direction
dy = Ly/(Ny-1)                                  # Length between two consecutive grid points in y-direction
dz = Lz/(Nz-1)                                  # Length between two consecutive grid points in z-direction         

X = np.arange(0,Nx)*dx                          # Array consisting of points in x-direction

if dimension == 3:
    Y = np.arange(0,Ny)*dy                      # Array consisting of points in y-direction

if grid_type == 'uniform':
    beta_s = 0
    Z = np.arange(0,Nz)*dz                      # Array consisting of points in z-direction

elif grid_type == 'non-uniform':
    beta_s = beta
    # Tangent-Hyperbolic grid for z-direction
    Z = (Lz/2)*(1-(np.tanh(beta*(1-2*(np.arange(0,Nz)*dz)))/(np.tanh(beta))))

if dimension == 2:
    Volume = Lx*Lz                              # Volume of the box
    Area_flux = Lx*np.ones_like(Z)              # Area of the flux surface perpendicular to z-axis
elif dimension == 3:
    Volume = Lx*Ly*Lz                           # Volume of the box
    Area_flux = Lx*Ly                           # Area of the flux surface perpendicular to z-axis

if dimension == 2:
    X_mesh, Z_mesh = np.meshgrid(X, Z,indexing = 'ij')                  # Meshgrids

elif dimension == 3:
    X_mesh, Y_mesh, Z_mesh = np.meshgrid(X, Y, Z,indexing = 'ij')       # Meshgrids  
    
# ------------------------------------------------------------------------------------------------------------------- #

# ------------------------------------- Simpson 1/3 rule to calculate integrals ------------------------------------- #

if dimension == 2:
    def simps1d(f,X):
        return np.sum((f[0:-2:2]+4*f[1:-1:2]+f[2::2])*(dx)/3)
        
    def simps2d(f,X,Z):
        if grid_type == 'uniform':
            return simps1d(np.sum((f[:,0:-2:2]+4*f[:,1:-1:2]+f[:,2::2])*(dz)/3,axis=1),X)
        elif grid_type == 'non-uniform':
            return simps1d(np.sum(((2-(Z[2::2]-Z[1:-1:2])/(Z[1:-1:2]-Z[:-2:2]))*f[:,0:-2:2] \
                                    + (2+(Z[2::2]-Z[1:-1:2])/(Z[1:-1:2]-Z[:-2:2])+(Z[1:-1:2]-Z[:-2:2])/(Z[2::2]-Z[1:-1:2]))*f[:,1:-1:2] \
                                    + (2-(Z[1:-1:2]-Z[:-2:2])/(Z[1:-1:2]-Z[:-2:2]))*f[:,2::2])*(Z[2::2]-Z[:-2:2])/6,axis=1),X)

elif dimension == 3:
    def simps1d(f,X):
        return np.sum((f[0:-2:2]+4*f[1:-1:2]+f[2::2])*(dx)/3)
        
    def simps2d(f,X,Y):
        return simps1d(np.sum((f[:,0:-2:2]+4*f[:,1:-1:2]+f[:,2::2])*(dy)/3,axis=1),X)
    
    def simps3d(f,X,Y,Z):
        if grid_type == 'uniform':
            return simps2d(np.sum((f[:,:,0:-2:2]+4*f[:,:,1:-1:2]+f[:,:,2::2])*dz/3,axis=2),X,Y)
        elif grid_type == 'non-uniform':
            return simps2d(np.sum(((2-(Z[2::2]-Z[1:-1:2])/(Z[1:-1:2]-Z[:-2:2]))*f[:,:,0:-2:2] \
                                    + (2+(Z[2::2]-Z[1:-1:2])/(Z[1:-1:2]-Z[:-2:2])+(Z[1:-1:2]-Z[:-2:2])/(Z[2::2]-Z[1:-1:2]))*f[:,:,1:-1:2] \
                                    + (2-(Z[1:-1:2]-Z[:-2:2])/(Z[1:-1:2]-Z[:-2:2]))*f[:,:,2::2])*(Z[2::2]-Z[:-2:2])/6,axis=2),X,Y)
        
# ------------------------------------------------------------------------------------------------------------------- #