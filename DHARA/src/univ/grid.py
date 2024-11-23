# ------------------------------------------------ Importing files -------------------------------------------------- #
from para import *
from DHARA.config.config import ncp
# ------------------------------------------------------------------------------------------------------------------- #

# ------------------------------------------------- Time variables -------------------------------------------------- #

dt_c = dt                                                       # Courant time step
t = ncp.arange(tinit,tfinal+dt,dt_c)                            # Time axis

Nf = int((tfinal-tinit)/save_fields_dt)                         # Number of field saving times
t_f_step = ncp.linspace(tinit,tfinal,Nf+1)                      # Field saving times

Ng = int((tfinal-tinit)/save_global_quantities_dt)              # Number of global quantities saving times
t_g_step = ncp.linspace(tinit,tfinal,Ng+1)                      # Global quantities saving times

Ngt = int((tfinal-tinit)/save_global_file)                      # Number of global file saving times
t_gt_step = ncp.linspace(tinit,tfinal,Ngt+1)                    # Global file saving times

s_j,s_l,s_m = 0,0,0                                             # For saving per time-step

# ------------------------------------------------------------------------------------------------------------------- #

# ------------------------------------------------- Grid variables -------------------------------------------------- #
 
Lx = A*Lz                                       # Length of box in x-direction
Ly = B*Lz                                       # Length of box in y-direction

r_b = 0                                         # Making sure that r_b = z(0) = 0 for Cartesian coordinates

dx = Lx/(Nx-1)                                  # Length between two consecutive grid points in x-direction
dy = Ly/(Ny-1)                                  # Length between two consecutive grid points in y-direction
dz = Lz/(Nz-1)                                  # Length between two consecutive grid points in z-direction         

X = ncp.arange(0,Nx)*dx                         # Array consisting of points in x-direction

if dimension == 3:
    Y = ncp.arange(0,Ny)*dy                     # Array consisting of points in y-direction

if grid_type == 'uniform':
    beta_s = 0
    Z = ncp.arange(0,Nz)*dz                     # Array consisting of points in z-direction

elif grid_type == 'non-uniform':
    beta_s = beta
    # Tangent-Hyperbolic grid for z-direction
    Z = (Lz/2)*(1-(ncp.tanh(beta*(1-2*(ncp.arange(0,Nz)*dz)))/(ncp.tanh(beta))))

if dimension == 2:
    Volume = Lx*Lz                              # Volume of the box
    Area_flux = Lx*ncp.ones_like(Z)             # Area of the flux surface perpendicular to z-axis
elif dimension == 3:
    Volume = Lx*Ly*Lz                           # Volume of the box
    Area_flux = Lx*Ly                           # Area of the flux surface perpendicular to z-axis

hz = ncp.zeros(Nz)                                  # Grid spacing in z-direction
hz[:-1] = Z[1:] - Z[:-1]                    
hz[-1] = hz[-2]

if dimension == 2:
    X_mesh, Z_mesh = ncp.meshgrid(X, Z,indexing = 'ij')                  # Meshgrids

elif dimension == 3:
    X_mesh, Y_mesh, Z_mesh = ncp.meshgrid(X, Y, Z,indexing = 'ij')       # Meshgrids  
    
# ------------------------------------------------------------------------------------------------------------------- #

# ----------------------------------------------- Boundary parameters ----------------------------------------------- #

T_t = 1 - (D+epsilon)                               # Temperature at top boundary
T_b = 1                                             # Temperature at bottom boundary

# ------------------------------------------------------------------------------------------------------------------- #

# ----------------------------------------------- Equilibrium profile ----------------------------------------------- #

alpha = 1/(gamma - 1)                               # Adiabatic index
m = (D*(1+alpha)/(D+epsilon)) - 1                   # Polytropic index

# ------------------------------------------------------------------------------------------------------------------- #

# --------------------------------------------- Constants in equations ---------------------------------------------- #

C1 = (1-r_b)/(epsilon*D*(alpha+1))                  # Normalized gas constant, p = C1*rho*T
C2 = (1-r_b)**(3/2)*ncp.sqrt(Pr/Ra)                 # Normalized dynamic viscoisty
C3 = 1/epsilon                                      # Normalized acceleration due to gravity
C4 = (1-r_b)**(5/2)/(epsilon*D*ncp.sqrt(Ra*Pr))     # Normalized thermal diffusivity

# ------------------------------------------------------------------------------------------------------------------- #

# ------------------------------------------------ Related parameters ----------------------------------------------- #

# Minimum grid spacing in z-direction
if grid_type == 'uniform':
    dz_min = dz
elif grid_type == 'non-uniform':
    dz_min = (1/2)*(1-ncp.tanh(beta*(1-2*dz))/(ncp.tanh(beta))) 
    
t_nu = ncp.sqrt(Ra/Pr)                                         # Viscous time scale
t_kappa = ncp.sqrt(Pr*Ra)                                      # Thermal diffusive time scale
Chi_0 = 1/(1-epsilon-D)**m                                     # Initial static density constrast (rho_b/rho_t)
Chi_a = 1/(1-D)**alpha                                         # Initial adiabatic density constrast (rho_b/rho_t)
Ma = ncp.sqrt(epsilon*D/(gamma-1))                             # Mach number, U/sqrt(gamma*R_gas*T_b) ; U = sqrt(epsilon*g*l)

# ------------------------------------------------------------------------------------------------------------------- #

# ----------------------------------------- Temporary and time series arrays ---------------------------------------- #

temp_dx = ncp.zeros_like(X_mesh)                    # Temporary array for derivatives calcualtion and other purposes
if dimension == 3:
    temp_dy = ncp.zeros_like(Y_mesh)                # Temporary array for derivatives calcualtion and other purposes
temp_dz = ncp.zeros_like(Z_mesh)                    # Temporary array for derivatives calcualtion and other purposes

# For storing parameters inside the field saving file
parameters = ncp.array([0, gamma, D, epsilon, Pr, Ra, beta_s])  

global_quan = ncp.zeros([11,len(t_g_step)])         # Global quantities at respective time       

eps = 1e-20                                         # Small number for avoiding division by zero

# ------------------------------------------------------------------------------------------------------------------- #