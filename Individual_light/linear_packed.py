#%% NEQ imports
import os
N_cores = 8
os.environ["OMP_NUM_THREADS"] = str(N_cores)        # export OMP_NUM_THREADS=4
os.environ["OPENBLAS_NUM_THREADS"] = str(N_cores)   # export OPENBLAS_NUM_THREADS=4
os.environ["MKL_NUM_THREADS"] = str(N_cores)        # export MKL_NUM_THREADS=6
os.environ["VECLIB_MAXIMUM_THREADS"] = str(N_cores) # export VECLIB_MAXIMUM_THREADS=4
os.environ["NUMEXPR_NUM_THREADS"] = str(N_cores)    # export NUMEXPR_NUM_THREADS=6

import numpy as np
import matplotlib.pyplot as plt

import jclsquant as jcl
import sys
sys.path.append('Code/Neq_TFM')

from ham_creation import create_graphene_ham
from lat_creation import get_positions_graphene
from core import DOS_sparse
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
# Benchmarking
from time import time
from datetime import timedelta
#%% NEQ Calcs
# ------------------------------------------------------------------------------
####  Defining a hamiltonian (own)
N = 2**20
N1 = N2 = int(np.sqrt(N))//2
S = get_positions_graphene(N1, N2)

Ham = create_graphene_ham(S, N1, N2, out_format='ELL')

dE = (Ham.bounds[1] - Ham.bounds[0])/2
# ------------------------------------------------------------------------------
#### PARAMETERS   
# Energy in pulse                        
E = 1            
# Temperature                       
Temp = 0.0
# Chemical potential
mu = 0.01
# Amount of periods to be simulated
n_periods = 20
# Amount of half multiples of E where the occupation is obtained
hE_reps = 2
# ??
tau = 0.0 
# Amount of measures per period
N_measures = 4*n_periods
# Parameters of the laser
w = E/jcl.hbar_fs 
T = 2*np.pi/w   
# Parameter of the wave package
Tp = T                         
# Time of  (fs)
t_vec = np.linspace(0,n_periods*T , 1000*n_periods)      
# Type of light               
modifier_id = 'linear_packed'
# Intensity param     (no units)
gamma = 0.025                   

# Intensity of the laser
Phi0 = jcl.hbar_fs*2*np.pi
A0 = gamma*Phi0/(2*3**0.5*jcl.a_cc)

# Amount of random vectors used in calculation
N_random_vector = 1

# Momenta
multM = 1.5
M = int(multM*np.sqrt(N))
# Broadening in the energies
broad = dE*np.pi/M
# Parameters of the light modifier
modifier_params = (A0, w, Tp)
obs_list = [['n', N_measures, M]]

# Names of file and info on graphs
folder_name = f'linear_packed/N={N}_E={E}_Temp={Temp}_mu={mu:.2f}_G={gamma:.3f}_M={multM:.1f}'
fig_title_info = f'N={N}, E={E} eV, T={Temp} K, $\\mu$={mu} eV, $\\Gamma$={gamma}'
pulse_suptitle = fr'$E={E}, \omega = {w:.3f}$ fs$^{{-1}}, T={T:.3f}$ fs'
extra_text = f'Light: linear\n$\\delta E={broad:.3f}$'

# Create said folder
os.makedirs(f'Out/{folder_name}', exist_ok=True)

# ------------------------------------------------------------------------------
## PRE-CALL: calculating the pulse in the given timeframe
Pulse = np.sin(t_vec*w)


cmap = plt.cm.plasma 
norm = Normalize(vmin=t_vec.min(), vmax=t_vec.max())

t_vec_measures = np.linspace(0, n_periods*T, N_measures)

fig, ax = plt.subplots()
ax.plot(t_vec/T, Pulse, color='blue')
for i in range(len(t_vec_measures)):
    ax.vlines(t_vec_measures[i]/T, -1, 1, color=cmap(norm(t_vec_measures[i])))
ax.set_xlabel('Time (Periods)')
ax.set_ylabel('Pulse amplitude')
fig.suptitle(pulse_suptitle)
ax.set_title('Pulse shape')
fig.savefig(f'Out/{folder_name}/PULSE.png', bbox_inches='tight')

taux = time()
print('Calculation starts!')
trep = time()
n_mat, dos_n_mat, t_vec_measures = jcl.kpm_rho_neq(Ham,t_vec,tau,modifier_id,modifier_params,Temp,mu,obs_list,M)
n_mat_total, dos_n_mat_total = np.zeros_like(n_mat), np.zeros_like(dos_n_mat)
print('Calc #1 done!')
print(f'Time elapsed: {timedelta(seconds=time() - trep)}')
for i in range(2, N_random_vector+1):
    trep = time()
    n_mat, dos_n_mat, t_vec_measures = jcl.kpm_rho_neq(Ham,t_vec,tau,modifier_id,modifier_params,Temp,mu,obs_list,M)
    n_mat_total += n_mat 
    dos_n_mat_total += dos_n_mat
    print(f'Calc #{i} done!')
    print(f'Time elapsed: {timedelta(seconds=time() - trep)}')
n_mat_total /= N_random_vector
dos_n_mat_total /= N_random_vector
print('Calculation finished!')
# Saving results
np.save(f'Out/{folder_name}/E.npy', n_mat_total[:,:,0])
np.save(f'Out/{folder_name}/n_E.npy', n_mat_total[:,:,1])
np.save(f'Out/{folder_name}/dosn_E.npy', dos_n_mat_total[:,:,1])
print(f'Time elapsed: {timedelta(seconds=time() - taux)}')


# ------------------------------------------------------------------------------
#%% Graph results
# Energy window and lines
min_e, max_e = -2.5, 2.5
hw_lines_step = 0.5
hw_hlines = [i*E for i in [hw_lines_step*i for i in range(1, int(max_e/hw_lines_step+1))]]
hw_hlines += [i*E for i in [-hw_lines_step*i for i in range(1, int(max_e/hw_lines_step+1))]]

# Text box
props = dict(boxstyle='round', facecolor='white', alpha=1.0)

fig, ax = plt.subplots()
for i in range(N_measures):
    ax.plot(n_mat[i,:,0], n_mat[i,:,1], color=cmap(norm(t_vec_measures[i])), label=f't={round(t_vec_measures[i]/T, 3)}T')
#ax.legend()
cbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax, orientation='vertical', label='Time (periods)')
ax.set_xlabel('Energy (eV)')
ax.set_ylabel('$n(\\varepsilon)$')
ax.set_xlim(min_e, max_e)
ax.vlines(hw_hlines, 0, 1, color='grey', ls='--', alpha=0.5, zorder=1)
fig.suptitle(fig_title_info)
ax.set_title('Occupation Number')
ax.text(min_e + 0.5, 0.2, extra_text, bbox=props)
fig.savefig(f'Out/{folder_name}/N(E).png', bbox_inches='tight')
plt.show()


#%% Extra graphs for the n for different times
hE_list = [hE*E/2 for hE in range(-hE_reps, hE_reps+1)]
total_nhE = 2*hE_reps + 1 
for (i, hE) in enumerate(hE_list):
    ind_time = np.where(n_mat[i,:,0] > hE)[0][0]
    n_timed = [n_mat[i,ind_time,1] for i in range(N_measures)]

    fig, ax = plt.subplots()
    ax.plot(np.array(t_vec_measures)/T, n_timed, c='blue', marker='.', ls='--')
    ax.set_xlabel('Time (Periods)')
    ax.set_ylabel('$n(\\varepsilon)$')
    ax.set_title(fr'Occupation in $E={hE/E:.1f}E$')
    fig.suptitle(fig_title_info)
    fig.savefig(f'Out/{folder_name}/N(T)_{hE:.1f}.png', bbox_inches='tight')
