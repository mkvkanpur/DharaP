#-------------------------------------- Libraries -----------------------------------------#
import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from pylab import rcParams
import pandas as pd
#------------------------------------------------------------------------------------------#

#----------------------------------- Import other files -----------------------------------#
from plot import *
from read_output import * 
#------------------------------------------------------------------------------------------#

mpl.style.use('classic')

if plot_type == 'error_total_mass':
    plt.figure(figsize=(5,3))
    plt.plot(t, np.abs(M_T-M_0)/M_0, label = r'$\frac{|M_T-M_0|}{M_0}$', color = 'r')
    plt.xlabel(r'$t$', fontsize=13)          
    plt.xlim(tinit,tfinal)
    # plt.ylim(1.6,1.8)      
    plt.legend(loc='lower right', prop={'size': 15}, frameon=False) 
    plt.tight_layout()
    plt.rc('text', usetex=True)          
    plt.savefig(output_folder_path + 'plots/error_total_mass.png', dpi =200)

elif plot_type == 'energy':
    plt.figure(figsize=(5,3))
    plt.loglog(t, Ke, label = r'$K_e = \frac{1}{2}\rho V^2$', color = 'r')
    plt.loglog(t, Ie, label = r'$I_e = c_v\rho T$', color = 'b')
    # plt.loglog(t, Ue, label = r'$U_e = \rho gz$', color = 'g')
    plt.xlabel(r'$t$', fontsize=13)        
    # plt.ylabel(r'$E$', fontsize=13)     
    plt.xlim(tinit,tfinal)      
    # plt.ylim(1e-13,1e1)
    plt.xticks([1,1e1,1e2],fontsize=10)
    plt.yticks([1e-9,1e-6,1e-3,1],fontsize=10)
    plt.legend(loc='lower right', prop={'size': 15}, frameon=False) 
    plt.tight_layout()
    plt.rc('text', usetex=True)          
    plt.savefig(output_folder_path + 'plots/energy.png', dpi =200)

elif plot_type == 'energies':
    fig  = plt.figure(figsize=(12,6), tight_layout = True)

    ax1 = fig.add_subplot(221)
    ax1.plot(t, Ke, label = r'$K_e = \frac{1}{2}\rho V^2$', color = 'r')
    ax1.set_xlabel(r'$t$', fontsize=13)          
    ax1.set_xlim(tinit,tfinal)
    plt.legend(loc='best', prop={'size': 12}, frameon=False, ncol=3) 

    ax2 = fig.add_subplot(222)
    ax2.plot(t, Ie, label = r'$I_e = c_v\rho T$', color = 'b')
    ax2.set_xlabel(r'$t$', fontsize=13)          
    ax2.set_xlim(tinit,tfinal) 
    # ax2.set_ylim(61.4,62)       
    plt.legend(loc='best', prop={'size': 12}, frameon=False, ncol=3) 

    ax3 = fig.add_subplot(223)
    ax3.plot(t, Ue, label = r'$U = \rho gz$', color = 'g')
    ax3.set_xlabel(r'$t$', fontsize=13)          
    ax3.set_xlim(tinit,tfinal)
    # ax3.set_ylim(4.15,4.24)       
    plt.legend(loc='best', prop={'size': 12}, frameon=False, ncol=3) 

    ax4 = fig.add_subplot(224)
    ax4.plot(t, E_T, label = r'$E_T = K_e+I_e+U$', color = 'k')
    ax4.set_xlabel(r'$t$', fontsize=13)          
    ax4.set_xlim(tinit,tfinal)
    # ax4.set_ylim(65.6,66.2)       
    plt.legend(loc='best', prop={'size': 12}, frameon=False, ncol=3) 
    
    fig.suptitle(r'$\mathrm{Energies \ vs \ time \ plots }$', fontsize=16)
    plt.savefig(output_folder_path + 'plots/energies.png', dpi =200)

elif plot_type == 'V_rms':
    plt.figure(figsize=(5,3))
    plt.plot(t, V_rms, color = 'k')
    plt.xlabel(r'$t$', fontsize=13)        
    plt.ylabel(r'$\langle V_{rms} \rangle$', fontsize=13)     
    plt.xlim(tinit,tfinal)      
    plt.plot([], [], ' ', label=r'$\mathrm{Re} = %.1f$' %(V_rms[500]/C2))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.ylim(0,4e-1)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.legend(loc=4, prop={'size': 15}, frameon=False) 
    plt.tight_layout()
    plt.rc('text', usetex=True)          
    plt.savefig(output_folder_path + 'plots/V_rms.png', dpi =200)

elif plot_type == 'Mach_number':
    plt.figure(figsize=(5,3))
    plt.plot(t, Ma, color = 'k')
    plt.xlabel(r'$t$', fontsize=13)        
    plt.ylabel(r'$Max \ Mach \ number$', fontsize=13)     
    plt.xlim(tinit, tfinal)      
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    # plt.ylim(0,1e-1)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.legend(loc=4, prop={'size': 15}, frameon=False) 
    plt.tight_layout()
    plt.rc('text', usetex=True)          
    plt.savefig(output_folder_path + 'plots/Mach_number.png', dpi =200)

elif plot_type == 'Nusselt':
    # Plotting Nusselt number    
    plt.figure(figsize=(5,3))
    plt.plot(t, Nu, label = r'$\mathrm{\overline{Nu}}$', color = 'b')
    plt.xlabel(r'$t$', fontsize=13)        
    plt.ylabel(r'$\mathrm{Nu}$', fontsize=13)     
    plt.xlim(0,100)      
    plt.ylim(1,5)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.legend(loc='best', prop={'size': 15}, frameon=False) 
    plt.tight_layout()
    plt.rc('text', usetex=True)          
    plt.savefig(output_folder_path + 'plots/Nusselt.png', dpi =200)