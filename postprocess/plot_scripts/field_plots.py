#-------------------------------------- Libraries -----------------------------------------#
import h5py
import math as math
import os
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from pylab import rcParams
mpl.rcParams.update(mpl.rcParamsDefault)
import matplotlib.ticker as ticker
#------------------------------------------------------------------------------------------#

#----------------------------------- Import other files -----------------------------------#
from plot import *
from read_output import * 
import grid 
#------------------------------------------------------------------------------------------#

if Nz <= 256:
    plot_grid_space = 1                          # Grid space for density plots 
else:
    plot_grid_space = int((Nz-1)/1024)                               
quiver_grid_space = int((Nz-1)/16)               # Grid space for quiver plots

Ra_n = int(np.log10(Ra))

bound_p = 2

mpl.style.use('classic')
plt.style.use('dark_background')

if plot_type == "field":
    with h5py.File(output_folder_path + "/fields/2D_%d.00.h5" %(t_save), "r") as f:
        # Get the data
        T = np.array(f['theta'])
        rho = np.array(f['rho'])
        ux = np.array(f['ux'])
        uz = np.array(f['uz'])
        param = np.array(f['parameters'])
        f.close()

    fig = plt.figure(figsize=(16,9), tight_layout = True)
    fig.suptitle(r'$\gamma = %.1f$,  ' %(gamma) + r'$\epsilon = %.1f$,  ' %(epsilon) + r'$D = %.1f$,  ' %(D) + r'$\mathrm{Pr} = %.1f$,  ' %(Pr) + r'$\mathrm{Ra} = 10^{%d}$' %(Ra_n) + '\n\n' + r'$t = %.2f$' %(t_save), fontsize=24, y=0.96)  
    
    # For temperature
    ax1 = fig.add_subplot(121)
    c1 = ax1.pcolor(grid.X_mesh[::plot_grid_space,::plot_grid_space], grid.Z_mesh[::plot_grid_space,::plot_grid_space], T[::plot_grid_space,::plot_grid_space], cmap=cm.inferno, vmin=-0.1, vmax=0)
    ax1.set_aspect(aspect = 1)
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('right', size='2%', pad=0.1)
    cb1 = fig.colorbar(c1, ticks=[0,-0.05,-0.1], cax=cax)
    cb1.ax.tick_params(labelsize=18)
    cb1.set_ticklabels([r'$0$',r'$-0.05$',r'$-0.1$'])
    ax1.quiver(grid.X_mesh[::quiver_grid_space,::quiver_grid_space],grid.Z_mesh[::quiver_grid_space,::quiver_grid_space],ux[::quiver_grid_space,::quiver_grid_space],uz[::quiver_grid_space,::quiver_grid_space],color = 'k',label='$\vec{u}$', scale=5, scale_units='xy', angles='xy', width=0.004)
    ax1.set_xlim(0,1)
    ax1.set_ylim(0,1)
    ax1.set_xticks([1],[r'$1$'], fontsize=20)
    ax1.set_yticks([0,1],[r'$0$',r'$1$'], fontsize=20)
    plt.gca().xaxis.set_minor_locator(ticker.FixedLocator([0,1,2]))
    plt.gca().xaxis.set_minor_formatter(ticker.NullFormatter())
    ax1.tick_params(axis='both', which='major', labelsize=18)
    ax1.set_xlabel(r'$x$', fontsize=20, labelpad=-1)
    ax1.set_ylabel(r'$z$', fontsize=20)
    ax1.set_title(r'$T-T_A(z)$', fontsize=22, pad=10)

    # For density
    ax2 = fig.add_subplot(122)
    temp = rho - (1-D*grid.Z_mesh)**grid.alpha
    c2 = ax2.pcolor(grid.X_mesh[::plot_grid_space,::plot_grid_space], grid.Z_mesh[::plot_grid_space,::plot_grid_space], temp[::plot_grid_space,::plot_grid_space], cmap=cm.seismic, vmin=-0.1, vmax=0.1)
    ax2.set_aspect(aspect = 1)
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes('right', size='2%', pad=0.1)
    cb2 = fig.colorbar(c2, ticks=[0.1,0,-0.1], cax=cax)
    cb2.ax.tick_params(labelsize=18)
    cb2.set_ticklabels([r'$0.1$',r'$0$',r'$-0.1$'])
    ax2.quiver(grid.X_mesh[::quiver_grid_space,::quiver_grid_space],grid.Z_mesh[::quiver_grid_space,::quiver_grid_space],ux[::quiver_grid_space,::quiver_grid_space],uz[::quiver_grid_space,::quiver_grid_space],color = 'k',label='$\vec{u}$', scale=5, scale_units='xy', angles='xy', width=0.004)
    ax2.set_xlim(0,1)
    ax2.set_ylim(0,1)
    ax2.set_xticks([1],[r'$1$'], fontsize=20)
    ax2.set_yticks([0,1],[r'$0$',r'$1$'], fontsize=20)
    ax2.tick_params(axis='both', which='major', labelsize=18)
    ax2.set_xlabel(r'$x$', fontsize=20, labelpad=-1)
    ax2.set_ylabel(r'$z$', fontsize=20)
    ax2.set_title(r'$\rho-\rho_A(z)$', fontsize=22, pad=10)

    plt.rc('text', usetex=True) 

    plt.savefig(output_folder_path + '/plots/density_plots/%.2f.png' %(t_save), dpi=200)
    plt.close(fig)
