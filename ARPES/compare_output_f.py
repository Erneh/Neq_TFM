#%%
import numpy as np
import matplotlib.pyplot as plt
from ARPES.kpath_stuff import get_path
from ARPES.auxiliar_functions import plot_pulse, load_data
import jclsquant as jcl

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
modifier_id = 'linear_packed'
# Energy in pulse                        
E = 1.0
# Temperature                       
Temp = 1e-9
# Chemical potential
mu = 0.01
# Intensity param     (no units)
gamma = 0.020

### SIMULATION
# Size of hamiltonian (2**N_pot)
N_pot = 16
N = 2**N_pot
# Amount of periods to be simulated
n_periods = 1
# Simulation steps per period
steps_per_T = 1000
# Amount of measures per period
meas_per_T = 40
N_measures = meas_per_T*n_periods
# Amount of random vectors used in calculation
N_random_vector = 1
# Momenta
M = int(np.sqrt(N))

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
print(f'Path chosen is: {path_type}')
print(f'nk is: {nk}')

# CREATION OF A K-PATH (in the original cell, I suppose?)


rLat, Rat, kpath, kind, kdist, klabs = get_path(path_type, nk)

# Reference hamiltonian
def H_og(k):
    """
    Basic hamiltonian of a graphene model
    """
    k_shape = np.shape(k)[1:]
    d1 = Rat[1]
    d2 = -rLat[0] + Rat[1]
    d3 = -rLat[1] + Rat[1]
    delta = np.exp(1j*(d1@k)) + np.exp(1j*(d2@k)) + np.exp(1j*(d3@k))
    H = np.array([[np.ones(k_shape)*mass, t*delta], 
                [t*np.conj(delta), -np.ones(k_shape)*mass]])
    return np.einsum('...i->i...', H)
#%%
# ------------------------------------------------------------------------------
# DATASET 1
# Things different from the general parameters
modifier_id = 'linear_packed'
mu = 0.00

label1 = 'linear_packed'
# Loading data
EF_list1, t_vec_meas1, n_f1, dosn_f1 = load_data(modifier_id, N_pot, E, Temp, mu, gamma, M, N_random_vector, 
              path_type, nk, n_periods, meas_per_T, steps_per_T, type_ham, mass)


# ------------------------------------------------------------------------------
# DATASET 2
# Things different from the general parameters
modifier_id = 'linear_packed'
mu = 5.00

# Loading data
EF_list2, t_vec_meas2, n_f2, dosn_f2 = load_data(modifier_id, N_pot, E, Temp, mu, gamma, M, N_random_vector, 
              path_type, nk, n_periods, meas_per_T, steps_per_T, type_ham, mass)


label2 = 'linear'

# Preparing data to be graphed
w = E/jcl.hbar_fs 
T = 2*np.pi/w
# Time of  (fs)

t_vec = np.linspace(0,n_periods*T , steps_per_T*n_periods)

dosn_f_mean1 = np.mean(dosn_f1, axis=0)
dosn_f_mean2 = np.mean(dosn_f2, axis=0)
"""
# Seeing light pulse
if modifier_id == 'circle':
    # Polarization (right or left)
    pol = 'r'
    modifier_stuff = (modifier_id, w, pol)

elif modifier_id == 'linear':
    modifier_stuff = (modifier_id, w)

elif modifier_id == 'linear_packed':
    Tp = T
    modifier_stuff = (modifier_id, w, Tp)

#fig, ax = plot_pulse(modifier_stuff, t_vec, t_vec_meas1)
"""

#%%
# Readying reference lines
autV, autE = np.linalg.eigh(H_og(kpath))
dosf_difference = dosn_f_mean1 - dosn_f_mean2

cmap = 'seismic'
col_min = -1e10
col_max = 1e10
levels = 400
col_levels = np.linspace(col_min, col_max, levels)
fig, ax = plt.subplots()
contour = ax.contourf(kdist, EF_list1, dosf_difference.T, col_levels, extend='both', cmap=cmap)
cbar = plt.colorbar(contour)
ax.set_xticks(kdist[kind], labels=klabs)
ax.set_ylabel('Energy')
#ax.plot(kdist[:,None], autV, c='pink')

# %%
