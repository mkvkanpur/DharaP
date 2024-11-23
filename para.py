# ---------------------------------------------------- Simulation details ---------------------------------------------------- #
 
output_dir = '2druns/Ra1e5/'     # Output directory

device = "CPU"            # Use Target to set device as CPU or GPU                                        
device_rank = 0           # Set GPU device rank     

dimension = 2             # Dimension of the simulation: 2, 3

grid_type = 'uniform'     # Grid type, uniform, non-uniform

# NS: No-slip, FS: Free-slip, PD: Periodic
boundary_u_x = 'PD'       # Boundary condtions for velocity in horizontal or x-direction
boundary_u_y = 'PD'       # Boundary condtions for velocity in horizontal or y-direction
boundary_u_z = 'NS'       # Boundary condtions for velocity in veritical or z-direction   

# FX: Fixed, AD: Adiabatic
boundary_T_x = 'PD'       # Boundary condtions for temperature in horizontal or x-direction
boundary_T_y = 'PD'       # Boundary condtions for temperature in horizontal or y-direction
boundary_T_z = 'FX'       # Boundary condtions for temperature in veritical or z-direction

# Note: For polar coordinates, keep theta or x-direction always periodic.

# Initial profiles; static: equilibrium profiles, custom: custom init.h5 in ../input/, continue: continue some run at time placed in tinit
profile = 'static'        

scheme = 'RK2'            # Time integration scheme, EULER, RK2

integration_scheme = 'trapezoidal'      # Integration schemes: trapezoidal or simpsons13

# --------------------------------------------------------------------------------------------------------------------------- #

# ------------------------------------------------ Time-stepping parameters ------------------------------------------------- #

tinit = 0                 # Initial time
tfinal = 100              # Final time
dt = 1e-3                 # Single time step

dt_far = 1e-3             # dt after transient state
t_far = 0                 # Time after transient state

# Note: After t > t_far, dt = dt_far. Use it to decrease dt after transient state.

cfl_cons = 0.5            # CFL constant for time step calculation

# --------------------------------------------------------------------------------------------------------------------------- #

# ----------------------------------------------------- Grid parameters ----------------------------------------------------- #

Lz = 1                    # Length of box in z-direction

# Note: Keep all number of grid points odd (For accurate Simpson's 1/3 rule) 

Nz = 65                   # Number of grid points in z-direction
Nx = 65                   # Number of grid points in x-direction
# Ny is for 3D simulations
Ny = 65                   # Number of grid points in y-direction

# For non-uniform grid
beta = 1                  # Stretching parameter, beta = 0 for uniform grid

# --------------------------------------------------------------------------------------------------------------------------- #

# --------------------------------------------------- Control parameters ---------------------------------------------------- #

A = 1                     # Aspect ratio in x-direction
B = 1                     # Aspect ratio in y-direction

epsilon = 0.1             # Superadiabaticity, Delta T/T_r, T_r is reference value
D = 0.5                   # Dissipation number
gamma = 1.4               # gamma = C_p/C_v 
Pr = 0.7                  # Prandtl number at reference boundary
Ra = 1e5                  # Rayleigh number at reference boundary

# --------------------------------------------------------------------------------------------------------------------------- #

# ---------------------------------------------------- Saving parameters ---------------------------------------------------- #

print_output_terminal = True        # Print running time-instance of the simulation in terminal

save_fields = True                  # Saving fields files
save_fields_dt = 5                  # Field files saving time-step

save_global_quantities = True       # Saving global quantitites
save_global_quantities_dt = 1e-3    # Global quantities saving time-step, keep it greater or equal to dt

save_global_file = 1                # Saving time of global file

# --------------------------------------------------------------------------------------------------------------------------- #