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
from tqdm import tqdm

import jclsquant as jcl

from ham_creation import create_hex_ham
from lat_creation import get_positions_graphene
from ARPES.kpath_stuff import rec_lattice, plot_1BZ, path_chart
from core import random_vector

mass = 0.5
t = -2.7

a_l = 0.24595
a1 = a_l*np.array([3**0.5/2, 1/2])
a2 = a_l*np.array([3**0.5/2, -1/2])
r1 = np.array([0.0, 0.0])
r2 = np.array([a_l/np.sqrt(3), 0.0])

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
N_pot = 16
N = 2**N_pot
N1 = N2 = int(np.sqrt(N))//2


S = get_positions_graphene(N1, N2)

Ham = create_hex_ham(S, N1, N2, t=t, M=mass, out_format='ELL')
positions = jcl.lattice_hexagonal(N)
Ham = jcl.H_graphene(positions, -2.7 + 0j, 0.5 + 0j,  periodic=True, type_H='ELL')
M = int(N**0.5)
#M = 1000
# Selecting the indices accordingly
index_list = [np.arange(0, N, 2), np.arange(1, N, 2)]

# ------------------------------------------------------------------------------
# CREATION OF A K-PATH (in the original cell, I suppose?)
rLat = np.array([a1, a2])
recLat, BZ_points = rec_lattice(rLat)
K = BZ_points[4]
Kp = BZ_points[5]
M_point = (K + Kp)/2
#plot_1BZ(recLat, BZ_points[[5],:])
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

EF_list = DOS_f[0,:,0]
fig, ax = plt.subplots()
ax.contourf(kdist, EF_list, np.log(DOS_f[:,:,1].T))
ax.set_xticks(kdist[kind], labels=klabs)
ax.plot(kdist[:,None], autV, c='pink')


# %%
