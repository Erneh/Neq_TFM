"""
This script relies on obtaining all its inputs via command line

"""

#%% NEQ imports
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
#%% Defining functions to check if the result has been already calculated

def check_if_calculated(modifier_id, N_pot, E, Temp, mu, gamma, N_random_vector,
                        n_periods, meas_per_T, steps_per_T):
    N = 2**N_pot
    # Names of file and info on graphs
    folder_name = f'{modifier_id}/N={N}_E={float(E)}_Temp={float(Temp)}_mu={mu:.2f}_G={gamma:.3f}/Nrand={N_random_vector}_nT={n_periods}_measT={meas_per_T}_stT={steps_per_T}'
    # Check if file exists
    try:
        os.makedirs(f'Out/{folder_name}', exist_ok=False)
        os.rmdir(f'Out/{folder_name}')
        flag = False
    except FileExistsError:
        flag = True

    # Check if the data is correctly saved in the files and computation is finished
    if flag:
        try:
            EF_list = np.load(f'Out/{folder_name}/E.npy')
            if np.sum(np.abs(EF_list)) == 0:
                flag = False
        except FileNotFoundError:
            flag = False
    return flag


#%% Reading command line input
# Type of Light
modifier_id = sys.argv[1]
# Power to which the number of atoms is 'powered'
N_pot = int(sys.argv[2])
# Energy in pulse                        
E = float(sys.argv[3])     
# Temperature                       
Temp = float(sys.argv[4])
# Chemical potential
mu = float(sys.argv[5])
# Intensity param     (no units)
gamma = float(sys.argv[6])
# Amount of moments used to calculate
M = int(sys.argv[7])
# Amount of random vectors used in calculation
N_random_vector = int(sys.argv[8])

# # periods included in sims
n_periods = int(sys.argv[9])
# Amount of measures per period
meas_per_T = int(sys.argv[10])
# steps/T
steps_per_T = int(sys.argv[11])
# Force recalculation (default is false)
try:
    force_recalc = bool(sys.argv[12])
except IndexError:
    force_recalc = False

# Call general function to get the job done
neq_sim(modifier_id, N_pot, E, Temp, mu, gamma, M, N_random_vector, 
            n_periods, meas_per_T, steps_per_T, force_recalc, False)