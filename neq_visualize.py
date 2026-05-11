<<<<<<< HEAD
#%% NEQ Calcs
=======
#%% NEQ imports
>>>>>>> Massive amount of changes and reorganizing of code. New results
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

from neq import neq_sim
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
# Benchmarking
from time import time
from datetime import timedelta
<<<<<<< HEAD

=======
#%% NEQ Calcs
>>>>>>> Massive amount of changes and reorganizing of code. New results
# ------------------------------------------------------------------------------
#### PARAMETERS   
# Type of Light
modifier_id = 'circle'
# Power to which the number of atoms is 'powered'
<<<<<<< HEAD
N_pot = 19
#M = int(np.sqrt(2**N_pot))
M = 362
=======
N_pot = 17
>>>>>>> Massive amount of changes and reorganizing of code. New results
# Energy in pulse                        
E = 1.0
# Temperature                       
Temp = 1e-9
# Chemical potential
mu = 0.01
# Intensity param     (no units)
gamma = 0.010
# Amount of random vectors used in calculation
<<<<<<< HEAD
N_random_vector = 1
=======
N_random_vector = 10
>>>>>>> Massive amount of changes and reorganizing of code. New results

# # periods included in sims
n_periods = 20
# Amount of measures per period
<<<<<<< HEAD
meas_per_T = 4
# steps/T
steps_per_T = 1000
# Force recalculation
force_recalc = True
show_figs = True


neq_sim(modifier_id, N_pot, E, Temp, mu, gamma, M, N_random_vector, 
=======
meas_per_T = 32
# steps/T
steps_per_T = 1000
# Force recalculation
force_recalc = False
show_figs = True


neq_sim(modifier_id, N_pot, E, Temp, mu, gamma, N_random_vector, 
>>>>>>> Massive amount of changes and reorganizing of code. New results
            n_periods, meas_per_T, steps_per_T, force_recalc, show_figs)

