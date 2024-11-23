# ------------------------------------------------ Importing files -------------------------------------------------- #
from DHARA.config.config import ncp
import DHARA.src.univ.grid as grid
from DHARA.src.lib.derivative import *
# ------------------------------------------------------------------------------------------------------------------- #

# --------------------------------------------- Compressible class -------------------------------------------------- #

class Compressible:
    # Class for all selfible variables

    def __init__(self):
        # Initialising all variables

        # Velocity vector components
        self.ux = []
        self.uy = []
        self.uz = []

        # Scalar fields
        self.rho = []
        self.T = []

        # Conserved variables
        self.Q = []
        self.F = []
        self.E = []

        # For time-stepping
        self.Q_copy1 = []
        
        # Temporary array
        self.temp = []

        self.F1 = []
        self.F2 = []
        self.F3 = []
        

    def set_arrays(self):
        # Setting arrays initially in the memory

        if para.dimension == 2:
            # Velocity vector components
            self.ux = ncp.zeros([grid.Nx,grid.Nz])
            self.uz = ncp.zeros([grid.Nx,grid.Nz])

            # Scalar fields
            self.rho = ncp.zeros([grid.Nx,grid.Nz])
            self.T = ncp.zeros([grid.Nx,grid.Nz])
        
            # Conserved variables
            self.Q = ncp.zeros([4,grid.Nx,grid.Nz])
            self.F = ncp.zeros([4,grid.Nx,grid.Nz])

            # For time-stepping
            self.Q_copy1 = ncp.zeros([4,grid.Nx,grid.Nz])
            
            # Temporary array
            self.temp = ncp.zeros([grid.Nx,grid.Nz])

            self.F1 = ncp.zeros([grid.Nx,grid.Nz])
            self.F2 = ncp.zeros([grid.Nx,grid.Nz])

        elif para.dimension == 3:
            # Velocity vector components
            self.ux = ncp.zeros([grid.Nx,grid.Ny,grid.Nz])
            self.uy = ncp.zeros([grid.Nx,grid.Ny,grid.Nz])
            self.uz = ncp.zeros([grid.Nx,grid.Ny,grid.Nz])

            # Scalar fields
            self.rho = ncp.zeros([grid.Nx,grid.Ny,grid.Nz])
            self.T = ncp.zeros([grid.Nx,grid.Ny,grid.Nz])
        
            # Conserved variables
            self.Q = ncp.zeros([5,grid.Nx,grid.Ny,grid.Nz])
            self.F = ncp.zeros([5,grid.Nx,grid.Ny,grid.Nz])
            
            self.Q_copy1 = ncp.zeros([5,grid.Nx,grid.Ny,grid.Nz])
            
            self.temp = ncp.zeros([grid.Nx,grid.Ny,grid.Nz])

            self.F1 = ncp.zeros([grid.Nx,grid.Ny,grid.Nz])
            self.F2 = ncp.zeros([grid.Nx,grid.Ny,grid.Nz])
            self.F3 = ncp.zeros([grid.Nx,grid.Ny,grid.Nz])


    def update_conserved(self):
        # Put the primitive variables inside Q's
        if para.dimension == 2:
            self.Q[0] = ncp.copy(self.rho)
            self.Q[1] = self.rho*self.ux
            self.Q[2] = self.rho*self.uz
            self.Q[3] = self.rho*(0.5*(self.ux*self.ux + self.uz*self.uz) + grid.C1/(para.gamma-1)*self.T + grid.C3*grid.Z_mesh) 
            
        elif para.dimension == 3:
            self.Q[0] = ncp.copy(self.rho)
            self.Q[1] = self.rho*self.ux
            self.Q[2] = self.rho*self.uy
            self.Q[3] = self.rho*self.uz
            self.Q[4] = self.rho*(0.5*(self.ux*self.ux + self.uy*self.uy + self.uz*self.uz) + grid.C1/(para.gamma-1)*self.T + grid.C3*grid.Z_mesh) 
            

    def compute_convective_flux_x(self):
        # Convective flux terms in x-direction
        if para.dimension == 2:
            self.F[0] = self.rho*self.ux
            self.F[1] = self.rho*(self.ux*self.ux + grid.C1*self.T)     # Pressure is given by p = C1*rho*T
            self.F[2] = self.rho*self.ux*self.uz
            self.F[3] = self.rho*self.ux*(0.5*(self.ux*self.ux + self.uz*self.uz) + (para.gamma/(para.gamma - 1))*grid.C1*self.T + grid.C3*grid.Z_mesh)
            
        elif para.dimension == 3:
            self.F[0] = self.rho*self.ux
            self.F[1] = self.rho*(self.ux*self.ux + grid.C1*self.T)     # Pressure is given by p = C1*rho*T
            self.F[2] = self.rho*self.ux*self.uy
            self.F[3] = self.rho*self.ux*self.uz
            self.F[4] = self.rho*self.ux*(0.5*(self.ux*self.ux + self.uy*self.uy + self.uz*self.uz) + (para.gamma/(para.gamma - 1))*grid.C1*self.T + grid.C3*grid.Z_mesh)
            

    def compute_convective_flux_y(self):
        # Convective flux terms in y-direction
        self.F[0] = self.rho*self.uy
        self.F[1] = self.rho*self.uy*self.ux
        self.F[2] = self.rho*(self.uy*self.uy + grid.C1*self.T)
        self.F[3] = self.rho*self.uy*self.uz
        self.F[4] = self.rho*self.uy*(0.5*(self.ux*self.ux + self.uy*self.uy + self.uz*self.uz) + (para.gamma/(para.gamma - 1))*grid.C1*self.T + grid.C3*grid.Z_mesh)
        

    def compute_convective_flux_z(self):
        # Convective flux terms in z-direction
        if para.dimension == 2:
            self.F[0] = self.rho*self.uz
            self.F[1] = self.rho*self.ux*self.uz
            self.F[2] = self.rho*(self.uz*self.uz + grid.C1*self.T)     # Pressure is given by p = C1*rho*T
            self.F[3] = self.rho*self.uz*(0.5*(self.ux*self.ux + self.uz*self.uz) + (para.gamma/(para.gamma - 1))*grid.C1*self.T + grid.C3*grid.Z_mesh)
            
        elif para.dimension == 3:
            self.F[0] = self.rho*self.uz
            self.F[1] = self.rho*self.uz*self.ux
            self.F[2] = self.rho*self.uz*self.uy
            self.F[3] = self.rho*(self.uz*self.uz + grid.C1*self.T)     # Pressure is given by p = C1*rho*T
            self.F[4] = self.rho*self.uz*(0.5*(self.ux*self.ux + self.uy*self.uy + self.uz*self.uz) + (para.gamma/(para.gamma - 1))*grid.C1*self.T + grid.C3*grid.Z_mesh)
            

    def compute_viscous_flux_x(self,d_type):
        # Viscous flux terms in x-direction

        if para.dimension == 2:
            dfx(self.ux,d_type)
            dfz(self.uz,0)
            self.F1 = grid.C2*(-4*grid.temp_dx + 2*grid.temp_dz)/3
            dfx(self.uz,d_type)
            dfz(self.ux,0)
            self.F2 = -grid.C2*(grid.temp_dx + grid.temp_dz)
            dfx(self.T,d_type)

            self.F[1] += self.F1
            self.F[2] += self.F2
            self.F[3] += self.ux*self.F1 + self.uz*self.F2 - grid.C4*grid.temp_dx
            
        elif para.dimension == 3:
            dfx(self.ux,d_type)
            dfy(self.uy,0)
            dfz(self.uz,0)
            self.F1 = grid.C2*(-4*grid.temp_dx + 2*grid.temp_dy + 2*grid.temp_dz)/3
            dfx(self.uy,d_type)
            dfy(self.ux,0)
            self.F2 = -grid.C2*(grid.temp_dx + grid.temp_dy)
            dfx(self.uz,d_type)
            dfz(self.ux,0)
            self.F3 = -grid.C2*(grid.temp_dx + grid.temp_dz)
            dfx(self.T,d_type)

            self.F[1] += self.F1
            self.F[2] += self.F2
            self.F[3] += self.F3
            self.F[4] += self.ux*self.F1 + self.uy*self.F2 + self.uz*self.F3 - grid.C4*grid.temp_dx
            
    
    def compute_viscous_flux_y(self,d_type):
        # Viscous flux terms in y-direction

        dfx(self.uy,0)
        dfy(self.ux,d_type)
        self.F1 = -grid.C2*(grid.temp_dx + grid.temp_dy)
        dfx(self.ux,0)
        dfy(self.uy,d_type)
        dfz(self.uz,0)
        self.F2 = grid.C2*(2*grid.temp_dx - 4*grid.temp_dy + 2*grid.temp_dz)/3 
        dfy(self.uz,d_type)
        dfz(self.uy,0)
        self.F3 = -grid.C2*(grid.temp_dy+ grid.temp_dz)
        dfy(self.T,d_type)

        self.F[1] += self.F1
        self.F[2] += self.F2
        self.F[3] += self.F3
        self.F[4] += self.ux*self.F1 + self.uy*self.F2 + self.uz*self.F3 - grid.C4*grid.temp_dy
        

    def compute_viscous_flux_z(self,d_type):
        # Viscous flux terms in z-direction

        if para.dimension == 2:
            dfx(self.uz,0)
            dfz(self.ux,d_type)
            self.F1 = -grid.C2*(grid.temp_dx + grid.temp_dz)
            dfx(self.ux,0)
            dfz(self.uz,d_type)
            self.F2 = grid.C2*(2*grid.temp_dx - 4*grid.temp_dz)/3
            dfz(self.T,d_type)

            self.F[1] += self.F1
            self.F[2] += self.F2
            self.F[3] += self.ux*self.F1 + self.uz*self.F2 - grid.C4*grid.temp_dz
            
        elif para.dimension == 3:
            dfx(self.uz,0)
            dfz(self.ux,d_type)
            self.F1 = -grid.C2*(grid.temp_dx + grid.temp_dz)
            dfy(self.uz,0)
            dfz(self.uy,d_type)
            self.F2 = -grid.C2*(grid.temp_dy + grid.temp_dz)
            dfx(self.ux,0)
            dfy(self.uy,0)
            dfz(self.uz,d_type)
            self.F3 = grid.C2*(2*grid.temp_dx + 2*grid.temp_dy - 4*grid.temp_dz)/3
            dfz(self.T,d_type)

            self.F[1] += self.F1
            self.F[2] += self.F2
            self.F[3] += self.F3
            self.F[4] += self.ux*self.F1 + self.uy*self.F2 + self.uz*self.F3 - grid.C4*grid.temp_dz
        

    def flux_derivative_x(self,d_type):
        # Derivatives of total flux terms in x-direction
        self.F[0] = ncp.copy(dfx(self.F[0],d_type))
        self.F[1] = ncp.copy(dfx(self.F[1],d_type))
        self.F[2] = ncp.copy(dfx(self.F[2],d_type))
        self.F[3] = ncp.copy(dfx(self.F[3],d_type))
        if para.dimension == 3:
            self.F[4] = ncp.copy(dfx(self.F[4],d_type))


    def flux_derivative_y(self,d_type):
        # Derivatives of total flux terms in y-direction
        self.F[0] = ncp.copy(dfy(self.F[0],d_type))
        self.F[1] = ncp.copy(dfy(self.F[1],d_type))
        self.F[2] = ncp.copy(dfy(self.F[2],d_type))
        self.F[3] = ncp.copy(dfy(self.F[3],d_type))
        self.F[4] = ncp.copy(dfy(self.F[4],d_type))   


    def flux_derivative_z(self,d_type):
        # Derivatives of total flux terms in z-direction
        self.F[0] = ncp.copy(dfz(self.F[0],d_type))
        self.F[1] = ncp.copy(dfz(self.F[1],d_type))
        self.F[2] = ncp.copy(dfz(self.F[2],d_type))
        self.F[3] = ncp.copy(dfz(self.F[3],d_type))
        if para.dimension == 3:
            self.F[4] = ncp.copy(dfz(self.F[4],d_type))


# ------------------------------------------------------------------------------------------------------------------- #

