import numpy as np
import matplotlib.pyplot as plt



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

from ham_creation import create_hex_ham
from lat_creation import get_positions_graphene
from core import DOS_sparse, check_if_calculated
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
# Benchmarking
from time import time
from datetime import timedelta

#%% NEQ Calcs in Full
# ------------------------------------------------------------------------------
####  Defining a hamiltonian (own)
N = 2**(3+11)
positions = jcl.lattice_hexagonal(N)
Ham = jcl.H_graphene(positions, -2.7 + 0j, 0 + 0j, periodic=True, type_H='ELL')

print(Ham.len_row)


#%% Test random stuff
from core import random_vector
N = 2**16
dt = 0.24595
R = random_vector(N, n_rand=1)[:,0]

Rk = np.fft.fft(R)

freqs = np.fft.fftfreq(N, dt)


fig, ax = plt.subplots()
ax.plot(freqs, np.abs(Rk))