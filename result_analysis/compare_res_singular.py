"""
Script made to analyze the output of the calculations. Will fail if the given
calcs do not exist!
"""
#%% First part: imports and analysis of just one specific case
import os
N_cores = 1
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
from core import DOS_sparse, load_data, str_parameters
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
# Benchmarking
from time import time
from datetime import timedelta

from core import frequency_analysis

# ------------------------------------------------------------------------------
#### PARAMETERS OF THE MODEL
## PHYSICAL
# Type of light               
modifier_id = 'linear'
# Hamiltonian type
type_ham = 'hbn'
# Parameters of the ham (only read if hbn)
ham_params = 0.5
# Energy in pulse                        
E = 1.1
# Temperature                       
Temp = 1e-9
# Chemical potential
mu = 0.01
# Intensity param     (no units)
gamma = 0.005     

## SIMULATION
# Size of hamiltonian
N_pot = 17
N = 2**N_pot
# Amount of periods to be simulated
n_periods = 100
# Simulation steps per period
steps_per_T = 1000
# Amount of measures per period
meas_per_T = 16
N_measures = meas_per_T*n_periods
# Amount of random vectors used in calculation
N_random_vector = 1
# Momenta
M = int(np.sqrt(N))
#M = 362

## RESULT ANALYSIS
# Range of searching the maximim frequency (in period^-1 units)
range_search = 1

## CALCULATED PARAMS
# Parameters of the laser
w = E/jcl.hbar_fs 
T = 2*np.pi/w    
t_vec = np.linspace(0,n_periods*T , steps_per_T*n_periods)   
# Amount of half multiples of E where the occupation is obtained
hE_reps = 2
# Broadening in the energies
t_vec_measures = np.linspace(0, n_periods*T, N_measures)


hE_list = [hE*E/2 for hE in range(-hE_reps, hE_reps+1)]



fig_title_info = f'N={2**N_pot}, T={Temp} K, $\\mu$={mu} eV, $\\Gamma$={gamma}, M={M}'
# ------------------------------------------------------------------------------
# First element
str1 = 'circle  16      0.5     1e-09   0.01    0.020   256     1   4       32      1000    basic   0.0     1'
modifier_id1, N_pot1, E1, Temp1, mu1, gamma1, M1, N_random_vector1, n_periods1, meas_per_T1, steps_per_T1, type_ham1, ham_param1 = str_parameters(str1)
# Loading the info in the .npy files
EF_list1, n_E_list1, dos_list1, dosn_list1 = load_data(modifier_id1, N_pot1, E1, 
    Temp1, mu1, gamma1, M1, N_random_vector1, n_periods1, meas_per_T1, steps_per_T1, 
    type_ham1, ham_param1, R=None)

# Readying identifier on graph
label1 = f"$\\hbar \\omega = {E1}$ eV"
color1 = 'blue'

# Parameters of the laser
N1 = 2**N_pot1
N_measures1 = meas_per_T1*n_periods1
w1 = E1/jcl.hbar_fs 
T1 = 2*np.pi/w1
t_vec1 = np.linspace(0,n_periods1*T1 , steps_per_T1*n_periods1)   
# Amount of half multiples of E where the occupation is obtained
hE_reps1 = 2
# Broadening in the energies
t_vec_measures1 = np.linspace(0, n_periods1*T1, N_measures1)


hE_list1 = [hE*E1/2 for hE in range(-hE_reps1, hE_reps1+1)]


# ------------------------------------------------------------------------------
# Second element
str2 = 'circle  16      1.0     1e-09   0.01    0.020   256     1   4       32      1000    basic   0.0     1'
modifier_id2, N_pot2, E2, Temp2, mu2, gamma2, M2, N_random_vector2, n_periods2, meas_per_T2, steps_per_T2, type_ham2, ham_param2 = str_parameters(str2)

# Loading the info in the .npy files
EF_list2, n_E_list2, dos_list2, dosn_list2 = load_data(modifier_id2, N_pot2, E2,
     Temp2, mu2, gamma2, M2, N_random_vector2, n_periods2, meas_per_T2, steps_per_T2, 
     type_ham2, ham_param2)

# Parameters of the laser
w2 = E/jcl.hbar_fs 
T2 = 2*np.pi/w2  
t_vec2 = np.linspace(0,n_periods2*T2 , steps_per_T2*n_periods2)  
N_measures2 = meas_per_T2*n_periods2 
# Broadening in the energies
t_vec_measures2 = np.linspace(0, n_periods2*T2, N_measures2)
hE_reps2 = 2

hE_list2 = [hE*E2/2 for hE in range(-hE_reps2, hE_reps2+1)]


label2 = f"$\\hbar \\omega = {E2}$ eV"
color2 = 'red'
# ------------------------------------------------------------------------------
# Third element
"""
N_pot = 17
mu = 0.50
# Names of file and info on graphs
folder_name = f'{modifier_id}/G={gamma:.3f}_E={float(E)}_Temp={Temp}_mu={mu:.2f}/N={N_pot}_M={M}_R={N_random_vector}_nT={n_periods}_measT={meas_per_T}_stT={steps_per_T}'

# Loading the info in the .npy files
EF_list3 = np.load(f'Out/{folder_name}/E.npy')
n_E_list3 = np.load(f'Out/{folder_name}/n_E.npy')

label3 = f"$\\mu={mu}$"
"""
# ------------------------------------------------------------------------------

# Energy window and lines
cmap = plt.cm.plasma 
norm = Normalize(vmin=t_vec.min(), vmax=t_vec.max())

# Extra graphs for the n for different times
hE_list = [hE*E/2 for hE in range(-hE_reps, hE_reps+1)]
color_list = ['olivedrab', 'red', 'orange', 'blue', 'darkviolet']
total_nhE = 2*hE_reps + 1 


occ_drop_list1, fourier_occ1, freq1, char_freq1, max_freq_ind1 = frequency_analysis(EF_list1, n_E_list1, hE_list1, t_vec_measures1, T1, range_search)
occ_drop_list2, fourier_occ2, freq2, char_freq2, max_freq_ind2 = frequency_analysis(EF_list2, n_E_list2, hE_list2, t_vec_measures2, T2, range_search)
#occ_drop_list3, fourier_occ3, freq3, char_freq3, max_freq_ind3 = frequency_analysis(EF_list3, n_E_list3, hE_list, t_vec_measures, T, range_search)



# General figure to contain all important data
fig, ax = plt.subplots()
ax.plot(np.array(t_vec_measures1), occ_drop_list1[3], c=color_list[2], 
        marker='.', ls='--', label=label1)
ax.plot(np.array(t_vec_measures2), occ_drop_list2[3], c=color_list[3], 
        marker='.', ls='--', label=label2)
#ax.plot(np.array(t_vec_measures)/T, occ_drop_list3[3], c=color_list[4], 
#        marker='.', ls='--', label=label3)
ax.set_xlabel('Time (Periods)')
ax.set_ylabel('$n(t)$')
ax.set_title(f"$E=0.5\\hbar \\omega$ eV")
fig.suptitle(fig_title_info)
ax.legend()
#%%

# General figure to contain all important data
fig, ax = plt.subplots()
ax.plot(np.array(t_vec_measures)/T, occ_drop_list1[4], c=color_list[2], 
        marker='.', ls='--', label=label1)
ax.plot(np.array(t_vec_measures)/T, occ_drop_list2[4], c=color_list[3], 
        marker='.', ls='--', label=label2)
#ax.plot(np.array(t_vec_measures)/T, occ_drop_list3[4], c=color_list[4], 
#        marker='.', ls='--', label=label3)
ax.set_xlabel('Time (Periods)')
ax.set_ylabel('$n(t)$')
ax.set_title(f"$E=1\\hbar \\omega$ eV")
fig.suptitle(fig_title_info)
ax.legend()

# ------------------------------------------------------------------------------
# FREQUENCY ANALYSIS OF THE RESULTS
# General figure
fig, ax = plt.subplots()
ax.plot(freq1*T/(2*np.pi), fourier_occ1[3], c=color_list[2], marker='.', ls='--', 
        label=f'{label1}, $\\omega_c={char_freq1[3]:.6f}$')
ax.plot(freq2*T/(2*np.pi), fourier_occ2[3], c=color_list[3], marker='.', ls='--', 
        label=f'{label2}, $\\omega_c={char_freq2[3]:.6f}$')
#ax.plot(freq3*T/(2*np.pi), fourier_occ3[3], c=color_list[4], marker='.', ls='--', 
#        label=f'{label3}, $\\omega_c={char_freq3[3]:.6f}$')
ax.set_xlabel('Ordinary frequency (period$^{-1}$)')
ax.set_ylabel('Amplitude')
ax.set_title('FFT in $E=0.5\\hbar \\omega$')
fig.suptitle(fig_title_info)
ax.legend()
#fig.savefig(f'Out/{folder_name}/FREQ_N(T).png', bbox_inches='tight')


