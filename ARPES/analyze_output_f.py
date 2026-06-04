#%%
import numpy as np
import matplotlib.pyplot as plt
from ARPES.kpath_stuff import get_path
from ARPES.auxiliar_functions import plot_pulse, load_data, string_to_parameters
import jclsquant as jcl

# ------------------------------------------------------------------------------
###### PARAMETERS OF THE  (manual)
### PHYSICAL
# Hamiltonian construction
mass = 0.50
t = -2.7
a_l = 0.24595
type_ham = 'hbn'
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
mu = 0.00
# Intensity param     (no units)
gamma = 0.020

### SIMULATION
# Size of hamiltonian (2**N_pot)
N_pot = 18
N = 2**N_pot
# Amount of periods to be simulated
n_periods = 1
# Simulation steps per period
steps_per_T = 1000
# Amount of measures per period
meas_per_T = 40
N_measures = meas_per_T*n_periods
# Amount of random vectors used in calculation
N_random_vector = 2
# Momenta
M = int(np.sqrt(N))

str_params = 'linear 18 1.0 1e-09 0.00 0.020 512 5 full 100 1 40 1000 basic 0.50'
modifier_id, N_pot, E, Temp, mu, gamma, M, N_random_vector, path_type, nk, n_periods, meas_per_T, steps_per_T, type_ham, mass = string_to_parameters(str_params)
fig_title = f'{modifier_id}, {type_ham}, $N={2**N_pot}$, $\\mu={mu}$ eV, $\\Gamma={gamma}$, $M={M}$'

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


# Loading data
EF_list, t_vec_meas, n_f, dosn_f = load_data(modifier_id, N_pot, E, Temp, mu, gamma, M, N_random_vector, 
              path_type, nk, n_periods, meas_per_T, steps_per_T, type_ham, mass)


# Preparing data to be graphed
w = E/jcl.hbar_fs 
T = 2*np.pi/w
# Time of  (fs)

t_vec = np.linspace(0,n_periods*T , steps_per_T*n_periods)

#%% Graphing stuff 
H_bounds = Ham.bounds
H_shape = Ham.shape
dosn_f_mean = np.mean(dosn_f[:meas_per_T], axis=0)/(H_bounds[1]*N**2)
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
plot_pulse(modifier_stuff, t_vec, t_vec_meas)

# Readying reference lines
autV, autE = np.linalg.eigh(H_og(kpath))

col_min = dosn_f_mean.min()
col_max = dosn_f_mean.max()/64
levels = 400
col_levels = np.linspace(col_min, col_max, levels)
fig, ax = plt.subplots()
contour = ax.contourf(kdist, EF_list, dosn_f_mean.T, col_levels, extend='max')
cbar = plt.colorbar(contour)
ax.set_xticks(kdist[kind], labels=klabs)
ax.set_ylabel('Energy')
fig.suptitle(fig_title)
#ax.plot(kdist[:,None], autV, c='pink')

# %% Postprocessing to obtain the actual bands instead of just the weird intensity part
thresh_dosn_f = 1e6

band_out = dosn_f_mean > thresh_dosn_f

# %%
