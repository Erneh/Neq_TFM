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

from ham_creation import create_graphene_ham
from lat_creation import get_positions_graphene
from core import DOS_sparse
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
# Benchmarking
from time import time
from datetime import timedelta

from core import frequency_analysis

# ------------------------------------------------------------------------------
####  Defining a hamiltonian (own)
N_pot = 17
N = 2**N_pot
N1 = N2 = int(np.sqrt(N))//2
S = get_positions_graphene(N1, N2)

Ham = create_graphene_ham(S, N1, N2, out_format='ELL')

dE = (Ham.bounds[1] - Ham.bounds[0])/2
# ------------------------------------------------------------------------------
#### PARAMETERS   
# Energy in pulse                        
E = 1            
# Temperature                       
Temp = 1e-9
# Chemical potential
mu = 0.01
# Amount of periods to be simulated
n_periods = 100
# Amount of half multiples of E where the occupation is obtained
hE_reps = 2
# ??
tau = 0.0 
# Amount of measures per period
meas_per_T = 8
N_measures = meas_per_T*n_periods
# Parameters of the laser
w = E/jcl.hbar_fs 
T = 2*np.pi/w                            
# Time of  (fs)
steps_per_T = 1000
t_vec = np.linspace(0,n_periods*T , steps_per_T*n_periods)      
# Type of light               
modifier_id = 'circle'
# Intensity param     (no units)
gamma = 0.010            

# Amount of random vectors used in calculation
N_random_vector = 1

# Momenta
#M = int(np.sqrt(N))
M = 362
# Broadening in the energies
broad = dE*np.pi/M
range_search = 1
t_vec_measures = np.linspace(0, n_periods*T, N_measures)
# ------------------------------------------------------------------------------
# First element
N_pot = 17
mu = 0.01
modifier_id = 'linear'
# Names of file and info on graphs
folder_name = f'{modifier_id}/G={gamma:.3f}_E={float(E)}_Temp={Temp}_mu={mu:.2f}/N={N_pot}_M={M}_R={N_random_vector}_nT={n_periods}_measT={meas_per_T}_stT={steps_per_T}'
fig_title_info = f'N={2**N_pot}, $\\hbar \\omega$={E} eV, T={Temp} K, $\\mu$={mu} eV, $\\Gamma$={gamma}, M={M}'
pulse_suptitle = fr'$E={E}, \omega = {w:.3f}$ fs$^{{-1}}, T={T:.3f}$ fs'
extra_text = f'Light: {modifier_id}\n$\\delta E={broad:.3f}$\n# Rand Vecs: {N_random_vector}\nMeasures/T={meas_per_T}\nSteps/T={steps_per_T}'

# Loading the info in the .npy files
EF_list1 = np.load(f'Out/{folder_name}/E.npy')
n_E_list1 = np.load(f'Out/{folder_name}/n_E.npy')

# Readying identifier on graph
label1 = f"linear"
# ------------------------------------------------------------------------------
# Second element
N_pot = 17
mu = 0.01
modifier_id = 'circle'
# Names of file and info on graphs
folder_name = f'{modifier_id}/G={gamma:.3f}_E={float(E)}_Temp={Temp}_mu={mu:.2f}/N={N_pot}_M={M}_R={N_random_vector}_nT={n_periods}_measT={meas_per_T}_stT={steps_per_T}'

# Loading the info in the .npy files
EF_list2 = np.load(f'Out/{folder_name}/E.npy')
n_E_list2 = np.load(f'Out/{folder_name}/n_E.npy')


label2 = f"circle"
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


occ_drop_list1, fourier_occ1, freq1, char_freq1, max_freq_ind1 = frequency_analysis(EF_list1, n_E_list1, hE_list, t_vec_measures, T, range_search)
occ_drop_list2, fourier_occ2, freq2, char_freq2, max_freq_ind2 = frequency_analysis(EF_list2, n_E_list2, hE_list, t_vec_measures, T, range_search)
#occ_drop_list3, fourier_occ3, freq3, char_freq3, max_freq_ind3 = frequency_analysis(EF_list3, n_E_list3, hE_list, t_vec_measures, T, range_search)


char_freq1*jcl.hbar_fs

# General figure to contain all important data
fig, ax = plt.subplots()
ax.plot(np.array(t_vec_measures)/T, occ_drop_list1[3], c=color_list[2], 
        marker='.', ls='--', label=label1)
ax.plot(np.array(t_vec_measures)/T, occ_drop_list2[3], c=color_list[3], 
        marker='.', ls='--', label=label2)
#ax.plot(np.array(t_vec_measures)/T, occ_drop_list3[3], c=color_list[4], 
#        marker='.', ls='--', label=label3)
ax.set_xlabel('Time (Periods)')
ax.set_ylabel('$n(t)$')
ax.set_title(f"$E=0.5\\hbar \\omega$ eV")
fig.suptitle(fig_title_info)
ax.legend()


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

#%% Second part of the analysis
# Plotting the characteristic frequency depending on the intensity
# ------------------------------------------------------------------------------
## SETTING COMMON CONSTANTS
N_pot = 17
N = 2**N_pot
E = 1            
# Temperature                       
Temp = 1e-9
# Chemical potential
mu = 0.01
# Amount of periods to be simulated
n_periods = 100
# Amount of half multiples of E where the occupation is obtained
hE_reps = 2
# ??
tau = 0.0 
# Amount of measures per period
meas_per_T = 8
N_measures = meas_per_T*n_periods
# Parameters of the laser
w = E/jcl.hbar_fs 
T = 2*np.pi/w                            
# Time of  (fs)
steps_per_T = 1000
t_vec = np.linspace(0,n_periods*T , steps_per_T*n_periods)      
# Type of light               
modifier_id = 'circle'
# Amount of random vectors used in calculation
N_random_vector = 1
range_search = 1
# Intensity param     (no units)
gamma_list = np.linspace(0.005, 0.025, 5)                   
hE_list = [hE*E/2 for hE in range(-hE_reps, hE_reps+1)]
color_list = ['olivedrab', 'red', 'orange', 'blue', 'darkviolet']
t_vec_measures = np.linspace(0, n_periods*T, N_measures)
M = int(N**0.5)

# Print of parameters to check results
print('PARAMETERS OF ANALYSIS')
print()
print(f'Type of light is {modifier_id}')
print(f'# of atoms: {N}')
print(f'Energy: {E} eV')
print(f'Temperature: {Temp} K')
print(f'Chem potential: {mu} eV')
print(f'# of periods: {n_periods}')
print(f'steps/period: {steps_per_T}')
print(f'# measures/T: {meas_per_T}')

char_freq1 = np.zeros((len(gamma_list), len(hE_list)))
char_freq2 = np.zeros((len(gamma_list), len(hE_list)))
for (g, gamma) in enumerate(gamma_list):
    if M == 0:
        M = int(np.sqrt(N))
    # ------------------------------------------------------------------------------
    # First element
    N_pot = 17
    mu = 0.01
    modifier_id = 'linear'
    # Names of file and info on graphs
    folder_name = f'{modifier_id}/G={gamma:.3f}_E={float(E)}_Temp={Temp}_mu={mu:.2f}/N={N_pot}_M={M}_R={N_random_vector}_nT={n_periods}_measT={meas_per_T}_stT={steps_per_T}'
    fig_title_info = f'N={2**N_pot}, $\\hbar \\omega$={E} eV, T={Temp} K, $\\mu$={mu} eV, $\\Gamma$={gamma}, M={M}'
    pulse_suptitle = fr'$E={E}, \omega = {w:.3f}$ fs$^{{-1}}, T={T:.3f}$ fs'
    extra_text = f'Light: {modifier_id}\n$\\delta E={broad:.3f}$\n# Rand Vecs: {N_random_vector}\nMeasures/T={meas_per_T}\nSteps/T={steps_per_T}'

    # Loading the info in the .npy files
    EF_list1 = np.load(f'Out/{folder_name}/E.npy')
    n_E_list1 = np.load(f'Out/{folder_name}/n_E.npy')

    # Readying identifier on graph
    label1 = f"linear"
    # ------------------------------------------------------------------------------
    # Second element
    N_pot = 17
    mu = 0.01
    modifier_id = 'circle'
    # Names of file and info on graphs
    folder_name = f'{modifier_id}/G={gamma:.3f}_E={float(E)}_Temp={Temp}_mu={mu:.2f}/N={N_pot}_M={M}_R={N_random_vector}_nT={n_periods}_measT={meas_per_T}_stT={steps_per_T}'

    # Loading the info in the .npy files
    EF_list2 = np.load(f'Out/{folder_name}/E.npy')
    n_E_list2 = np.load(f'Out/{folder_name}/n_E.npy')


    label2 = f"circle"
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
    
    occ_drop_list1, fourier_occ1, freq1, char_freq1[g], max_freq_ind1 = frequency_analysis(EF_list1, n_E_list1, hE_list, t_vec_measures, T, range_search)
    occ_drop_list2, fourier_occ2, freq2, char_freq2[g], max_freq_ind2 = frequency_analysis(EF_list2, n_E_list2, hE_list, t_vec_measures, T, range_search)
    # --------------------------------------------------------------------------
    # OCCUPATION(time) GRAPH
    # General figure to contain all important data
    fig, ax = plt.subplots()
    ax.plot(np.array(t_vec_measures)/T, occ_drop_list1[3], c='orange', 
            marker='.', ls='--', label=f'{label1}')
    ax.plot(np.array(t_vec_measures)/T, occ_drop_list2[3], c='blue', 
            marker='.', ls='--', label=label2)
    ax.set_xlabel('Time (Periods)')
    ax.set_ylabel('$n(t)$')
    ax.set_title(fr' in E=0.5$\hbar\omega$')
    fig.suptitle(fig_title_info)
    ax.legend()
    # --------------------------------------------------------------------------
    # FREQUENCY ANALYSIS OF THE RESULTS
    # General figure

    fig, ax = plt.subplots()
    ax.plot(freq1*T/(2*np.pi), fourier_occ1[3], c='orange', marker='.', ls='--', 
        label=f'{label1}, $\\omega_c={char_freq1[g, 3]:.6f}$')
    ax.plot(freq2*T/(2*np.pi), fourier_occ2[3], c='blue', marker='.', ls='--',
        label=f'{label2}, $\\omega_c={char_freq2[g, 3]:.6f}$')
    # Markers of max frequencies
    ax.scatter(freq1[max_freq_ind1[3]]*T/(2*np.pi), fourier_occ1[3, max_freq_ind1[3]], color='red',
        marker='*', zorder=2)
    ax.scatter(freq2[max_freq_ind2[3]]*T/(2*np.pi), fourier_occ2[3, max_freq_ind2[3]], color='cyan',
        marker='*', zorder=2)
    ax.set_xlabel('Normal Frequency (period$^{-1}$)')
    ax.set_ylabel('Amplitude')
    ax.set_title('FFT')
    fig.suptitle(fig_title_info)
    ax.legend()

    fig, ax = plt.subplots()
    ax.plot(freq1*T/(2*np.pi), fourier_occ1[4], c='orange', marker='.', ls='--', 
        label=f'{label1}, $\\omega_c={char_freq1[g, 4]:.6f}$')
    ax.plot(freq2*T/(2*np.pi), fourier_occ2[4], c='blue', marker='.', ls='--',
        label=f'{label2}, $\\omega_c={char_freq2[g, 4]:.6f}$')
    # Markers of max frequencies
    ax.scatter(freq1[max_freq_ind1[4]]*T/(2*np.pi), fourier_occ1[4, max_freq_ind1[4]], color='red',
        marker='*', zorder=2)
    ax.scatter(freq2[max_freq_ind2[4]]*T/(2*np.pi), fourier_occ2[4, max_freq_ind2[4]], color='cyan',
        marker='*', zorder=2)
    ax.set_xlabel('Normal Frequency (period$^{-1}$)')
    ax.set_ylabel('Amplitude')
    ax.set_title('FFT')
    fig.suptitle(fig_title_info)
    ax.legend()
    

# Final plot of the gammas
fig, ax = plt.subplots()
#ax.plot(gamma_list, char_freq[:,2], ls='--', c=color_list[2], marker='.',
#        label=f'$E = \\mu$ eV')
ax.plot(gamma_list, char_freq1[:,3], ls='--', c='orange', marker='.',
        label=label1)
ax.plot(gamma_list, char_freq2[:,3], ls='--', c='blue', marker='.',
        label=label2, alpha=0.8)
ax.set_xlabel(f'Intensity $\\Gamma$')
ax.set_ylabel('Angular frequency (las. period$^{-1}$)')
ax.set_title(f'Characteristic frequencies for different intensities ($E=0.5\\hbar \\omega$)')
ax.legend()
ax.set_xticks(gamma_list)
fig.suptitle(f'$N=2^{{{N_pot}}}$, $\\hbar \\omega$={E} eV, Temp={Temp} K, $\\mu$={mu} eV')

# Final plot of the gammas
fig, ax = plt.subplots()
#ax.plot(gamma_list, char_freq[:,2], ls='--', c=color_list[2], marker='.',
#        label=f'$E = \\mu$ eV')
ax.plot(gamma_list, char_freq1[4], ls='--', c='orange', marker='.',
        label=label1)
ax.plot(gamma_list, char_freq2[:,4], ls='--', c='blue', marker='.',
        label=label2, alpha=0.8)
ax.set_xlabel(f'Intensity $\\Gamma$')
ax.set_ylabel('Angular frequency (las. period$^{-1}$)')
ax.set_title(f'Characteristic frequencies for different intensities ($E=1\\hbar \\omega$)')
ax.legend()
ax.set_xticks(gamma_list)
fig.suptitle(f'$N=2^{{{N_pot}}}$, $\\hbar \\omega$={E} eV, Temp={Temp} K, $\\mu$={mu} eV')
