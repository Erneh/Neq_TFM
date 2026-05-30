#%%
import numpy as np
import matplotlib.pyplot as plt
from ARPES.kpath_stuff import path_chart, rec_lattice, load_data

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
modifier_id = 'linear'
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
n_periods = 5
# Simulation steps per period
steps_per_T = 1000
# Amount of measures per period
meas_per_T = 4
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
a1 = a_l*np.array([1/2, 3**0.5/2])
a2 = a_l*np.array([-1/2, 3**0.5/2])
r1 = np.array([0.0, 0.0])
r2 = np.array([a_l/np.sqrt(3), 0.0])
Rat = np.array([r1, r2])

rLat = np.array([a1, a2])
recLat, BZ_points = rec_lattice(rLat)
K = BZ_points[4]
Kp = BZ_points[5]
M_point = (K + Kp)/2
Gamma = np.array([0.0, 0.0])
if path_type == 'full':
    kpoints = [Gamma, K, M_point, Kp, Gamma]
    kpath, kind, kdist = path_chart(kpoints, nk, recLat)
    klabs = ['$\\Gamma$', '$K$', '$M$', "$K'$", '$\\Gamma$']
elif path_type == 'part':
    kpoints = [Gamma, K, M_point, Gamma]
    kpath, kind, kdist = path_chart(kpoints, nk, recLat)
    klabs = ['$\\Gamma$', '$K$', '$M$', '$\\Gamma$']
elif path_type == 'vall':
    kpoints = [Gamma, K, Gamma]
    kpath, kind, kdist = path_chart(kpoints, nk, recLat)
    klabs = ['$\\Gamma$', '$K$', '$\\Gamma$']

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


# Loading data
EF_list, t_vec_meas, n_f, dosn_f = load_data(modifier_id, N_pot, E, Temp, mu, gamma, M, N_random_vector, 
              path_type, nk, n_periods, meas_per_T, steps_per_T, type_ham, mass)


# Preparing data to be graphed
dosn_f_mean = np.mean(dosn_f, axis=0)
#%% Graphing stuff 
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
#ax.plot(kdist[:,None], autV, c='pink')

# %%
