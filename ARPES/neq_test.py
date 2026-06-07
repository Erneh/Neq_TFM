#%%
import os
N_cores = 8
os.environ["OMP_NUM_THREADS"] = str(N_cores)        # export OMP_NUM_THREADS=4
os.environ["OPENBLAS_NUM_THREADS"] = str(N_cores)   # export OPENBLAS_NUM_THREADS=4
os.environ["MKL_NUM_THREADS"] = str(N_cores)        # export MKL_NUM_THREADS=6
os.environ["VECLIB_MAXIMUM_THREADS"] = str(N_cores) # export VECLIB_MAXIMUM_THREADS=4
os.environ["NUMEXPR_NUM_THREADS"] = str(N_cores)    # export NUMEXPR_NUM_THREADS=6

import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('/home/eperez/Code/Neq_TFM')

import jclsquant as jcl

from ham_creation import create_hex_ham
from lat_creation import get_positions_graphene
from ARPES.kpath_stuff import get_path
from ARPES.auxiliar_functions import plot_pulse
from core import random_vector
# ------------------------------------------------------------------------------
###### PARAMETERS OF THE SYSTEM
### PHYSICAL
# Hamiltonian construction
mass = 0.0
t = -2.7
a_l = 0.24595
type_ham = 'basic'
# kpath n
path_type = 'full'
nk = 100
# Type of light               
modifier_id = 'circle'
# Energy in pulse                        
E = 1.0
# Temperature                       
Temp = 1e-9
# Chemical potential
mu = 0.00
# Intensity param     (no units)
gamma = 0.025

### SIMULATION
# Size of hamiltonian (2**N_pot)
N_pot = 18
N = 2**N_pot
# Amount of periods to be simulated
n_periods = 4
# Simulation steps per period
steps_per_T = 500
# Amount of measures per period
meas_per_T = 4
N_measures = meas_per_T*n_periods
# Amount of random vectors used in calculation
N_random_vector = 1
# Momenta
M = int(np.sqrt(N))

# ------------------------------------------------------------------------------
# CREATION OF THE HAMILTONIAN

# Print of parameters to check results
print('PARAMETERS OF CALCULATION')
print(f'Type of hamiltonian is {type_ham}')
print(f'Hamiltonian mass is {mass}')
print(f'Type of light is {modifier_id}')
print(f'# of atoms: {N}')
print(f'# of moments: {M}')
print(f'Energy: {E} eV')
print(f'Intensity param: {gamma}')
print(f'Temperature: {Temp} K')
print(f'Chem potential: {mu} eV')
print(f'# of Random Vectors: {N_random_vector}')
print(f'# of periods: {n_periods}')
print(f'steps/period: {steps_per_T}')
print(f'# measures/T: {meas_per_T}')

if type_ham == 'basic':
    type_ham = ''
    N1 = N2 = int(np.sqrt(N))//2
    S = get_positions_graphene(N1, N2)
    Ham = create_hex_ham(S, N1, N2, out_format='ELL')
    folder_name = f'{modifier_id}{type_ham}/{path_type}/G={gamma:.3f}_E={float(E)}_Temp={Temp}_mu={mu:.2f}/N={N_pot}_M={M}_nT={n_periods}_measT={meas_per_T}_stT={steps_per_T}_nk={nk}'

elif type_ham == 'jcl':
    positions = jcl.lattice_hexagonal(N)
    Ham = jcl.H_graphene(positions, -2.7 + 0j, periodic=True, type_H='ELL')
    folder_name = f'{modifier_id}{type_ham}/G={gamma:.3f}_E={float(E)}_Temp={Temp}_mu={mu:.2f}/N={N_pot}_M={M}_nT={n_periods}_measT={meas_per_T}_stT={steps_per_T}_nk={nk}'

elif type_ham == 'hbn':
    if mass is None:
        mass = 0.5

    N1 = N2 = int(np.sqrt(N))//2
    S = get_positions_graphene(N1, N2, a_l = a_l)
    Ham = create_hex_ham(S, N1, N2, t=-2.7, M=mass, a_l=0.25, out_format='ELL')
    folder_name = f'{modifier_id}{type_ham}/{path_type}/G={gamma:.3f}_E={float(E)}_Temp={Temp}_mu={mu:.2f}_m={mass:.2f}/N={N_pot}_M={M}_nT={n_periods}_measT={meas_per_T}_stT={steps_per_T}_nk={nk}'

dE = (Ham.bounds[1] - Ham.bounds[0])/2
# ------------------------------------------------------------------------------
#### PARAMETERS   
hE_reps = 2
# ??
tau = 0.0 

N_measures = meas_per_T*n_periods
# Parameters of the laser
w = E/jcl.hbar_fs 
T = 2*np.pi/w
# Time of  (fs)

t_vec = np.linspace(0,n_periods*T , steps_per_T*n_periods)      
t_vec_meass = np.linspace(0, n_periods*T, meas_per_T*n_periods) 
# Intensity of the laser
#Phi0 = jcl.hbar_fs*2*np.pi
#A0 = gamma*Phi0/(2*3**0.5*jcl.a_cc)
A0 = np.pi*gamma/(3**0.5*jcl.a_cc)
# Broadening in the energies
broad = dE*np.pi/M

# Different possibilities depending on type of light
if modifier_id == 'circle':
    # Polarization (right or left)
    pol = 'r'
    modifier_params = (A0, w, pol)

elif modifier_id == 'linear':
    modifier_params = (A0, w)

elif modifier_id == 'linear_packed':
    Tp = T
    modifier_params = (A0, w, Tp)

elif modifier_id == 'circle_packed':
    Tp = T
    pol = 'r'
    modifier_params = (A0, w, pol, Tp)


#M = 1000

# Selecting the indices accordingly
ar1 = np.array([1, 0, 1, 0], dtype=bool)
ind1 = np.kron(np.ones(N//4, dtype=bool), ar1)
ind2 = np.kron(np.ones(N//4, dtype=bool), np.bool(1-ar1))
total = np.arange(N)
index_list = np.array([total[ind1], total[ind2]])


a1 = a_l*np.array([1/2, 3**0.5/2])
a2 = a_l*np.array([-1/2, 3**0.5/2])
r1 = np.array([0.0, 0.0])
r2 = np.array([a_l/np.sqrt(3), 0.0])
Rat = np.array([r1, r2])

# Reference hamiltonian
def H_og(k):
    """
    Basic hamiltonian of a graphene model
    """
    k_shape = np.shape(k)[1:]
    d1 = r2
    d2 = -a1 + r2 
    d3 = -a2 + r2
    delta = np.exp(1j*(d1@k)) + np.exp(1j*(d2@k)) + np.exp(1j*(d3@k))
    H = np.array([[np.ones(k_shape)*mass, t*delta], 
                  [t*np.conj(delta), -np.ones(k_shape)*mass]])
    return np.einsum('...i->i...', H)

plot_pulse(modifier_id, modifier_params, t_vec, t_vec_meass)
# ------------------------------------------------------------------------------
# CREATION OF A K-PATH (in the original cell, I suppose?)
path_type = 'papr'
nk = 100
rLat, Rat, kpath, kind, kdist, klabs = get_path(path_type, nk)



obs_list = [['n', N_measures, M, kpath.T],
    ['n_f', N_measures, M, kpath.T, S, index_list]]

#%%
n_mat,dos_n_mat,t_vec_meass_n,n_mat_f,dos_n_mat_f,t_vec_meass_n_f = jcl.kpm_rho_neq_f(Ham,t_vec,0.0,modifier_id=modifier_id,modifier_params=modifier_params,Temp=Temp,mu=mu,observale_list=obs_list,M=M)

#%% All Times
EF_list = dos_n_mat_f[0,0,:,0]
#dosn_f_mean = np.mean(dos_n_mat_f[:meas_per_T], axis=0)


H_bounds = Ham.bounds
H_shape = Ham.shape

dos_n_mat_f_norm = dos_n_mat_f[:,:,:,1]/(H_bounds[1]*N**2)


fig_title = f'{modifier_id}, $N={N}$, $\\Gamma={gamma}$, $M={M}$'

for i in range(len(t_vec_meass_n_f)):
    col_min = 0
    col_max = 1
    levels = 400

    col_levels = np.linspace(col_min, col_max, levels)
    fig, ax = plt.subplots()
    contour = ax.contourf(kdist, EF_list, dos_n_mat_f_norm[i,:,:].T, col_levels, extend='both')
    ax.set_title(f't={t_vec_meass_n_f[i]/T}')
    #ax.imshow((dosn_f_mean.T)[::-1], aspect='auto')
    cbar = plt.colorbar(contour)
    #ax.imshow(dos_n_mat_f[i,:,:,1].T/(H_bounds[1]*N**2), extent=[kdist[0], kdist[-1], dos_n_mat_f[i,:,:,0].min(), dos_n_mat_f[i,:,:,0].max()],origin='lower', aspect='auto', cmap='viridis')
    ax.set_xticks(kdist[kind], labels=klabs)
    ax.set_ylabel('Energy')
    fig.suptitle(fig_title)
    ax.set_ylim(-1.5, 1.5)
    plt.show()
    #ax.plot(kdist[:,None], autV, c='pink')


#%% Mean Time
EF_list = dos_n_mat_f[0,0,:,0]
#dosn_f_mean = np.mean(dos_n_mat_f[:meas_per_T], axis=0)

H_bounds = Ham.bounds
H_shape = Ham.shape

dos_n_mat_f_norm = dos_n_mat_f[:,:,:,1]/(H_bounds[1]*N**2)




for i in range(1, n_periods+1):
    dosn_f_mean = np.mean(dos_n_mat_f_norm[(i-1)*meas_per_T:i*meas_per_T,:,:], axis=0)

    col_min = 0
    col_max = 1
    levels = 400
    col_levels = np.linspace(col_min, col_max, levels)

    fig_title = f'{modifier_id}, $N={N}$, $\\Gamma={gamma}$, $M={M}$'
    fig, ax = plt.subplots()
    contour = ax.contourf(kdist, EF_list, dosn_f_mean.T, col_levels, extend='both')
    #ax.imshow((dosn_f_mean.T)[::-1], aspect='auto')
    cbar = plt.colorbar(contour)
    #ax.imshow(dos_n_mat_f[i,:,:,1].T/(H_bounds[1]*N**2), extent=[kdist[0], kdist[-1], dos_n_mat_f[i,:,:,0].min(), dos_n_mat_f[i,:,:,0].max()],origin='lower', aspect='auto', cmap='viridis')
    ax.set_xticks(kdist[kind], labels=klabs)
    ax.set_ylabel('Energy')
    fig.suptitle(fig_title)
    ax.set_title(f'Mean between periods {i-1} and {i}')
    ax.set_ylim(-1.5, 1.5)
    plt.show()
    #ax.plot(kdist[:,None], autV, c='pink')

## Saving results
#print('Results are being saved...')
#save_path = f'ARPES/Out/{folder_name}/dosn'
#os.makedirs(save_path, exist_ok=True)
## Read existing files in the folder
#
## Only one is needed for each saved stuff
#np.save(f'{save_path}/E.npy', EF_list)
#np.save(f'{save_path}/dosn_full.npy', dos_n_mat_f)



# %%
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

cmap = plt.cm.plasma 
norm = Normalize(vmin=t_vec.min()/T, vmax=t_vec.max()/T)

min_e, max_e = -2.5, 2.5
hw_lines_step = 0.5*E
hw_hlines = [i for i in [hw_lines_step*i for i in range(1, int(max_e/hw_lines_step+1))]]
hw_hlines += [i for i in [-hw_lines_step*i for i in range(1, int(max_e/hw_lines_step+1))]]
n_E_list = n_mat[:,:, 1]
dosn_list = dos_n_mat[:,:,1]
# Text box
props = dict(boxstyle='round', facecolor='white', alpha=0.8)
# ------------------------------------------------------------------------------
# Occupation (E)
# Text box
fig, ax = plt.subplots()
for i in range(N_measures):
    ax.plot(EF_list, n_E_list[i,:], color=cmap(norm(t_vec_meass_n[i]/T)), label=f't={round(t_vec_meass_n[i]/T, 3)}T')
#ax.legend()
cbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax, orientation='vertical', label='Time (periods)')
ax.set_xlabel('Energy (eV)')
ax.set_ylabel('$n(\\varepsilon)$')
ax.set_xlim(min_e, max_e)
ax.vlines(hw_hlines, 0, 1, color='grey', ls='--', alpha=0.5, zorder=1)
ax.set_title('Occupation Number')

# ------------------------------------------------------------------------------
# DOS*Occ (Energy)
fig, ax = plt.subplots()
for i in range(N_measures):
    ax.plot(EF_list, dosn_list[i,:], color=cmap(norm(t_vec_meass_n[i])), label=f't={round(t_vec_meass_n[i]/T, 3)}T')
#ax.legend()
cbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax, orientation='vertical', label='Time (periods)')
ax.set_xlabel('Energy (eV)')
ax.set_ylabel('$n(\\varepsilon)$')
#ax.vlines(hw_hlines, 0, np.max(dos_list), color='grey', ls='--', alpha=0.5, zorder=1)
ax.vlines([0], 0, np.max(dosn_list), color='grey', ls='--', alpha=0.8, zorder=1)
ax.set_title('Density of States')


# %%
