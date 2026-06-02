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


import jclsquant as jcl

from ham_creation import create_hex_ham
from lat_creation import get_positions_graphene
from ARPES.kpath_stuff import rec_lattice, plot_1BZ, path_chart
from core import random_vector

mass = 0.5
t = -2.7

a_l = 0.24595
a1 = a_l*np.array([1/2, 3**0.5/2])
a2 = a_l*np.array([-1/2, 3**0.5/2])
r1 = np.array([0.0, 0.0])
r2 = np.array([a_l/np.sqrt(3), 0.0])
Rat = np.array([r1, r2])

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
# CREATION OF THE HAMILTONIAN
N_pot = 18
ham_type = 'jcl'
N = 2**N_pot
if ham_type == 'basic':
    N1 = N2 = int(np.sqrt(N))//2
    S = get_positions_graphene(N1, N2)
    Ham = create_hex_ham(S, N1, N2, t=t, M=mass, out_format='ELL')
elif ham_type == 'jcl':
    S = jcl.lattice_hexagonal(N)
    Ham = jcl.H_graphene(S, -2.7 + 0j, mass + 0j,  periodic=True, type_H='ELL')
    M = int(N**0.5)
#M = 1000
# Selecting the indices accordingly
ar1 = np.array([1, 0, 1, 0], dtype=bool)
ind1 = np.kron(np.ones(N//4, dtype=bool), ar1)
ind2 = np.kron(np.ones(N//4, dtype=bool), np.bool(1-ar1))
total = np.arange(N)
index_list = np.array([total[ind1], total[ind2]])

# ------------------------------------------------------------------------------
# CREATION OF A K-PATH (in the original cell, I suppose?)
rLat = np.array([a1, a2])
recLat, BZ_points = rec_lattice(rLat)
K = BZ_points[4]
Kp = BZ_points[5]
M_point = (K + Kp)/2
Gamma = np.array([0.0, 0.0])
nk = 100
kpoints = [Gamma, K, M_point, Kp, Gamma]
kpath, kind, kdist = path_chart(kpoints, nk, recLat)
klabs = ['$\\Gamma$', '$K$', '$M$', "$K'$", '$\\Gamma$']

# %% Performing the actual calculations
rnd_vec = random_vector(N, 1)[:,0]
DOS_f = jcl.kpm_dos_f(Ham, M, kpath.T, S, index_list, rnd_vec)


# %%
autV, autE = np.linalg.eigh(H_og(kpath))

col_min = DOS_f[:,:,1].T.min()
col_max = DOS_f[:,:,1].T.max()/64
levels = 400
col_levels = np.linspace(col_min, col_max, levels)

fig_title = f'{ham_type},N={N}, $m={mass}$, $M={M}$'
EF_list = DOS_f[0,:,0]
fig, ax = plt.subplots()
contour = ax.contourf(kdist, EF_list, DOS_f[:,:,1].T, col_levels, extend='max')
cbar = plt.colorbar(contour)
ax.set_xticks(kdist[kind], labels=klabs)
ax.plot(kdist[:,None], autV, c='gray', ls='--')
fig.suptitle(fig_title)

# %%
