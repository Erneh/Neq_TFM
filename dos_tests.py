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
from time import time

# Defining a hamiltonian (jcl)
N = 2**20
S_jcl = jcl.lattice_hexagonal(N)
H_jcl = jcl.H_graphene(S_jcl, -2.7 + 0j, 0.0 + 0j, 0.0 + 0j, True, 'ELL')
dE = (H_jcl.bounds[1] - H_jcl.bounds[0])/2

# Calculating DOS
# res = jcl.kpm_dos(H_jcl)
# EF_list = res[:,0]
# DOS_list = res[:,1]/dE/N

fig, ax = plt.subplots()
ax.plot(EF_list, DOS_list, c='blue')
ax.set_xlabel('Energy (eV)')
ax.set_ylabel('DOS (not normalized)')

np.trapezoid(DOS_list, EF_list)


# Defining a hamiltonian (own)
N = 2**17
N1 = N2 = int(np.sqrt(N))//2
S = get_positions_graphene(N1, N2)

Ham = create_graphene_ham(S, N1, N2, out_format='ELL')

dE = (Ham.bounds[1] - Ham.bounds[0])/2
res = jcl.kpm_dos(Ham)
EF_list = res[:,0]
DOS_list = res[:,1]/dE/N

#DOS_mine = DOS_sparse(jcl.ell_to_csr(Ham), Ham.bounds, EF_list, int(np.sqrt(N)))
fig, ax = plt.subplots()
ax.plot(EF_list, DOS_list, c='blue')
ax.set_xlabel('Energy (eV)')
ax.set_ylabel('DOS (not normalized)')
#%% MAIN NEQ CALCULATIONS
#### FIXED 'FREE' PARAMETERS                           
E = 1                                   # Energy in pulse? (eV)
Temp = 0.0
mu = 0.0      

#### Partially fixed parameters
w = E/jcl.hbar_fs 
T = 2*np.pi/w
Tp = T                                  # Pulse length parameter
n_periods = 4

t_vec = np.linspace(0,n_periods*T , 1000*n_periods)      # Timeframe of sim    (fs)
tau = 0.0                                       # ??
modifier_id = 'linear_packed'                   # Chosen light pol
gamma = 0.025                                 # Intensity param     (no units)
Phi0 = jcl.hbar_fs*2*np.pi
A0 = gamma*Phi0/(2*3**0.5*jcl.a_cc)

M = 2*int(np.sqrt(N))
N_measures = 4


modifier_params = (A0, w, Tp)
obs_list = [['n', N_measures, M]]

# Name of saved graph
graph_title = f'N={N}_E={E}_Temp={Temp}_mu={mu:.2f}_G={gamma:.3f}'
fig_title = f'N={N}, E={E} eV, T={Temp} K, $\\mu$={mu} eV, $\\Gamma$={gamma}'

##Pre-call: calculating the pulse in the given timeframe
Pulse = np.sin(t_vec*w)/np.cosh((t_vec - 2*Tp)/0.5673/Tp)
Pulse_wrap = 1/np.cosh((t_vec - 2*Tp)/0.5673/Tp)

cmap = plt.cm.plasma 
norm = Normalize(vmin=t_vec.min(), vmax=t_vec.max())

t_vec_measures = np.linspace(0, n_periods*T, N_measures)

fig, ax = plt.subplots()
ax.plot(t_vec, Pulse, color='blue')
# Graph the limits of the wave packet involved
ax.plot(t_vec, Pulse_wrap, ls='--', color='gray')
ax.plot(t_vec, -Pulse_wrap, ls='--', color='gray')
for i in range(len(t_vec_measures)):
    ax.vlines(t_vec_measures[i], -1, 1, color=cmap(norm(t_vec_measures[i])))
ax.set_xlabel('Time (fs)')
ax.set_ylabel('Pulse amplitude')
ax.set_title('Pulse shape')
fig.savefig(f'Out/PULSE_{graph_title}.png', bbox_inches='tight')

taux = time()
print('Calculation starts!')
n_mat, dos_n_mat, t_vec_measures = jcl.kpm_rho_neq(Ham,t_vec,tau,modifier_id,modifier_params,Temp,mu,obs_list,M)


# Graph results
min_e, max_e = -2.5, 2.5
hw_lines_step = 0.5
hw_hlines = [i*E for i in [hw_lines_step*i for i in range(1, int(max_e/hw_lines_step+1))]]
hw_hlines += [i*E for i in [-hw_lines_step*i for i in range(1, int(max_e/hw_lines_step+1))]]

fig, ax = plt.subplots()
for i in range(N_measures):
    ax.plot(n_mat[i,:,0], n_mat[i,:,1], color=cmap(norm(t_vec_measures[i])), label=f't={round(t_vec_measures[i]/T, 3)}')
ax.legend()
ax.set_xlabel('Energy (eV)')
ax.set_ylabel('$n(\\varepsilon)$')
ax.set_xlim(min_e, max_e)
ax.vlines(hw_hlines, 0, 1, color='grey', ls='--', alpha=0.5, zorder=1)
ax.set_title(fig_title)
fig.savefig(f'Out/N_{graph_title}.png', bbox_inches='tight')
plt.show()
print('Calculation finished and plot graphed!')
print(f'Time elapsed: {time() - taux}')