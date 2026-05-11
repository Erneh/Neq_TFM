"""
Script made to analyze the output of the calculations. Will fail if the given
calcs do not exist!
"""
#%% First part: imports and analysis of just one specific case
import os
N_cores = 1
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


# ------------------------------------------------------------------------------
####  Defining a hamiltonian (own)
N_pot = 23
N = 2**N_pot
N1 = N2 = int(np.sqrt(N))//2
S = get_positions_graphene(N1, N2)

Ham = create_graphene_ham(S, N1, N2, out_format='ELL')

dE = (Ham.bounds[1] - Ham.bounds[0])/2
# ------------------------------------------------------------------------------
#### PARAMETERS   
# Energy in pulse                        
E = 1            
# Temperature                       
Temp = 1e-9
# Chemical potential
mu = 0.01
# Amount of periods to be simulated
n_periods = 100
# Amount of half multiples of E where the occupation is obtained
hE_reps = 2
# ??
tau = 0.0 
# Amount of measures per period
meas_per_T = 4
N_measures = meas_per_T*n_periods
# Parameters of the laser
w = E/jcl.hbar_fs 
T = 2*np.pi/w                            
# Time of  (fs)

steps_per_T = 1000
t_vec = np.linspace(0,n_periods*T , steps_per_T*n_periods)      
# Type of light               
modifier_id = 'circle'
# Intensity param     (no units)
gamma = 0.025            


# Intensity of the laser
Phi0 = jcl.hbar_fs*2*np.pi
A0 = gamma*Phi0/(2*3**0.5*jcl.a_cc)

# Amount of random vectors used in calculation
N_random_vector = 1

# Momenta
M = int(np.sqrt(N))
# Broadening in the energies
broad_neq = dE*np.pi/M

t_vec_measures = np.linspace(0, n_periods*T, N_measures)
# Print of parameters to check results
print('PARAMETERS OF ANALYSIS')
print()
print(f'Type of light is {modifier_id}')
print(f'# of atoms: {N}')
print(f'Energy: {E} eV')
print(f'Intensity param: {gamma}')
print(f'Temperature: {Temp} K')
print(f'Chem potential: {mu} eV')
print(f'# of periods: {n_periods}')
print(f'steps/period: {steps_per_T}')
print(f'# measures/T: {meas_per_T}')


# Names of file and info on graphs
folder_name = f'{modifier_id}/G={gamma:.3f}_E={float(E)}_Temp={Temp}_mu={mu:.2f}/N={N_pot}_M={M}_R={N_random_vector}_nT={n_periods}_measT={meas_per_T}_stT={steps_per_T}'
fig_title_info = f'N={N}, E={E} eV, T={Temp} K, $\\mu$={mu} eV, $\\Gamma$={gamma}'
pulse_suptitle = fr'$E={E}, \omega = {w:.3f}$ fs$^{{-1}}, T={T:.3f}$ fs'
extra_text = f'$\\delta E={broad_neq:.3f}$\n# Rand Vecs: {N_random_vector}'

# Loading the info in the .npy files
EF_list = np.load(f'Out/{folder_name}/E.npy')
dos_E_list = np.load(f'Out/{folder_name}/dosn_E.npy')

# Normalizing the DOS to put it in equal terms to the other
C_norm = np.trapezoid(dos_E_list, EF_list)
dos_E_norm = dos_E_list / C_norm[:,None]


cmap = plt.cm.plasma 
norm = Normalize(vmin=t_vec.min()/T, vmax=t_vec.max()/T)
min_e, max_e = -1.5, 1.5
# Find the max of the DOS in the selected energy range
inds_min_e = EF_list > min_e
inds_max_e = EF_list < max_e

inds_range_e = inds_min_e * inds_max_e
DOS_indmax = np.max(dos_E_norm[inds_range_e])

hw_lines_step = 0.5
hw_hlines = [i*E for i in [hw_lines_step*i for i in range(1, int(max_e/hw_lines_step+1))]]
hw_hlines += [i*E for i in [-hw_lines_step*i for i in range(1, int(max_e/hw_lines_step+1))]]

fig, ax = plt.subplots()
for i in range(N_measures):
    ax.plot(EF_list[i,:], dos_E_norm[i,:], color=cmap(norm(t_vec_measures[i]/T)), label=f't={round(t_vec_measures[i]/T, 3)}T')
#ax.legend()
cbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax, orientation='vertical', label='Time (periods)')
ax.set_xlabel('Energy (eV)')
ax.set_ylabel('$n(\\varepsilon)$')
ax.set_xlim(min_e, max_e)
ax.vlines(hw_hlines, 0, np.max(dos_E_norm), color='grey', ls='--', alpha=0.5, zorder=1)
ax.vlines([0], 0, np.max(dos_E_norm), color='grey', ls='--', alpha=0.8, zorder=1)
fig.suptitle(fig_title_info)
ax.set_title('Density of States')
ax.set_ylim((0, DOS_indmax))
#fig.savefig(f'Out/{folder_name}/dos(E).png', bbox_inches='tight')

# We get the mean and represent 
N_measurements, n_E = dos_E_norm.shape
dos_E_mean = np.zeros(n_E)
count = 0
for i in range(N_measurements):
    #if i % 4 == 0:
        count += 1 
        dos_E_mean += dos_E_norm[i,:]
dos_E_mean /= count

# Text box
props = dict(boxstyle='round', facecolor='white', alpha=1.0)

fig, ax = plt.subplots()
for i in range(N_measures):
    ax.plot(EF_list[0], dos_E_mean, color='blue')
#ax.legend()
ax.set_xlabel('Energy (eV)')
ax.set_ylabel('$n(\\varepsilon)$')
ax.set_xlim(min_e, max_e)
ax.vlines(hw_hlines, 0, np.max(dos_E_norm), color='grey', ls='--', alpha=0.5, zorder=1)
ax.vlines([0], 0, np.max(dos_E_norm), color='grey', ls='--', alpha=0.8, zorder=1)
ax.text(min_e + 0.5, 0.1, extra_text, bbox=props)
fig.suptitle(fig_title_info)
ax.set_title('Density of States')
ax.set_ylim((0, DOS_indmax))
#fig.savefig(f'Out/{folder_name}/dos(E).png', bbox_inches='tight')


#%% ----------------------------------------------------------------------------
# COMPARATION OF THE DOS WITH MY OWN FLOQUET
# Floquet 
model_select = 'real'
broad = 0.005
min_e, max_e = -8.5, 8.5
dE_factor = 1.5
rep = 0
hw = 1.0
n_rep = 4

dE = broad / dE_factor
EF_fl = np.arange(min_e, max_e, dE)

folder_name = f'br={float(broad)}_ran={min_e}_{max_e}_dEf={dE_factor}_hw={float(hw)}'
subname = f'G={gamma}_rep={rep}'
out_folder = f'{folder_name}/{subname}_{model_select}'

DOS_fl = np.load(rf'/home/eperez/Code/Floquet_tfm/Outr/{out_folder}/dos_nrep={n_rep}.npy')

# We do the mean between all the different taken measurements and graph again

dos_E_mean = np.sum(dos_E_norm, axis=0) / N_measurements
EF_mean = np.sum(EF_list, axis=0) / N_measurements

# Take the maximum values for the density of states 
DOS_flmax = np.max(DOS_fl)
DOS_neqmax = np.max(dos_E_norm)
DOS_max = np.max([DOS_flmax, DOS_neqmax])

fig, ax = plt.subplots()
ax.plot(EF_mean, dos_E_mean, color='blue', label=f'NEQ ($\\delta E = {broad_neq:.3f}$)')
ax.plot(EF_fl, DOS_fl, color='red', label=f'Floquet ($\\eta = {broad:.3f}$)')
ax.set_xlabel('Energy (eV)')
ax.set_ylabel('$n(\\varepsilon)$')
ax.set_xlim(min_e, max_e)
ax.vlines(hw_hlines, 0, DOS_max, color='grey', ls='--', alpha=0.5, zorder=1)
ax.vlines([0], 0, DOS_max, color='grey', ls='--', alpha=0.8, zorder=1)
fig.suptitle(fig_title_info)
ax.set_title('DOS comparison')
ax.legend()
#fig.savefig(f'Out/{folder_name}/dos(E).png', bbox_inches='tight')