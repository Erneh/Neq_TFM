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
from ARPES.kpath_stuff import rec_lattice, plot_1BZ, path_chart
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
nk = 10
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
N_pot = 16
N = 2**N_pot
# Amount of periods to be simulated
n_periods = 1
# Simulation steps per period
steps_per_T = 1000
# Amount of measures per period
meas_per_T = 16
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


# ------------------------------------------------------------------------------
# CREATION OF A K-PATH (in the original cell, I suppose?)
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


obs_list = [['n_f', N_measures, M, kpath.T, S, index_list]]
n_mat_f,dos_n_mat_f,t_vec_meass_n_f = jcl.kpm_rho_neq_f(Ham,t_vec,0.0,modifier_id=modifier_id,modifier_params=modifier_params,Temp=Temp,mu=mu,observale_list=obs_list,M=M)

# Processing the data to reduce saved stuff
EF_list = dos_n_mat_f[0,0,:,0]
dosn_mean = np.mean(dos_n_mat_f, axis=0)

## Saving results
print('Results are being saved...')
save_path = f'ARPES/Out/{folder_name}/dosn'
os.makedirs(save_path, exist_ok=True)
# Read existing files in the folder

# Only one is needed for each saved stuff
np.save(f'{save_path}/E.npy', EF_list)
np.save(f'{save_path}/dosn_full.npy', dos_n_mat_f)



# %%
