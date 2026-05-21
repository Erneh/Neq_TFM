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

#from ham_creation import create_graphene_ham
from lat_creation import get_positions_graphene
from core import DOS_sparse, load_data
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
gamma_list = np.linspace(0.000, 0.025, 6)                   

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
color_list = ['olivedrab', 'red', 'orange', 'blue', 'darkviolet']
t_vec_measures = np.linspace(0, n_periods*T, N_measures)


# Plotting the characteristic frequency depending on the intensity
# ------------------------------------------------------------------------------
## SETTING COMMON CONSTANTS

# Print of parameters to check results


char_freq1 = np.zeros((len(gamma_list), len(hE_list)))
char_freq2 = np.zeros((len(gamma_list), len(hE_list)))
for (g, gamma) in enumerate(gamma_list):
    fig_title_info = f'N={2**N_pot}, $\\hbar \\omega$={E} eV, T={Temp} K, $\\mu$={mu} eV, $\\Gamma$={gamma}, M={M}'
    if M == 0:
        M = int(np.sqrt(N))
    # ------------------------------------------------------------------------------
    # First element01
    modifier_id = 'linear'
    # Names of file and info on graphs
    # Loading the info in the .npy files
    EF_list1, n_E_list1, dos_list1, dosn_list1 = load_data(modifier_id, N_pot, E, Temp, mu, gamma, 
                        M, N_random_vector, n_periods, meas_per_T, steps_per_T, type_ham, ham_params, R=None)

    # Readying identifier on graph
    label1 = f"linear"
    # ------------------------------------------------------------------------------
    # Second element
    modifier_id = 'circle'
    # Loading the info in the .npy files
    EF_list2, n_E_list2, dos_list2, dosn_list2 = load_data(modifier_id, N_pot, E, Temp, mu, gamma, 
                        M, N_random_vector, n_periods, meas_per_T, steps_per_T, type_ham, ham_params, R=None)

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
    
    occ_drop_list1, fourier_occ1, freq1, char_freq1[g], max_freq_ind1 = frequency_analysis(EF_list1, dosn_list1, hE_list, t_vec_measures, T, range_search)
    occ_drop_list2, fourier_occ2, freq2, char_freq2[g], max_freq_ind2 = frequency_analysis(EF_list2, dosn_list2, hE_list, t_vec_measures, T, range_search)
    # --------------------------------------------------------------------------
    # FREQUENCY ANALYSIS OF THE RESULTS
    # General figure
    """
    # --------------------------------------------------------------------------
    # OCCUPATION(time) 0.5HW
    # General figure to contain all important data
    fig, ax = plt.subplots()
    ax.plot(np.array(t_vec_measures)/T, occ_drop_list1[3], c='blue', 
            marker='.', ls='--', label=f"$E=0.5\\hbar\\omega$ eV")
    ax.plot(np.array(t_vec_measures)/T, occ_drop_list2[3], c='red', 
            marker='.', ls='--', label=f"$E=0.5\\hbar\\omega$ eV")
    ax.set_xlabel('Time (Periods)')
    ax.set_ylabel('$n(t)$')
    ax.set_title(fr'Occupation')
    fig.suptitle(fig_title_info)
    ax.legend()
    
    # --------------------------------------------------------------------------
    # OCCUPATION(time)  1.0HW
    # General figure to contain all important data
    fig, ax = plt.subplots()
    ax.plot(np.array(t_vec_measures)/T, occ_drop_list1[4], c='darkviolet', 
            marker='.', ls='--', label=f"$E=1\\hbar\\omega$")
    ax.plot(np.array(t_vec_measures)/T, occ_drop_list2[4], c='darkviolet', 
            marker='.', ls='--', label=f"$E=1\\hbar\\omega$")
    ax.set_xlabel('Time (Periods)')
    ax.set_ylabel('$n(t)$')
    ax.set_title(fr'Occupation')
    fig.suptitle(fig_title_info)
    ax.legend()

    # --------------------------------------------------------------------------
    # OCCUPATION(time)  0.0HW
    # General figure to contain all important data
    fig, ax = plt.subplots()
    ax.plot(np.array(t_vec_measures)/T, occ_drop_list1[2], c='orange', 
            marker='.', ls='--', label=f"$E=0.0\\hbar\\omega$ eV")
    ax.plot(np.array(t_vec_measures)/T, occ_drop_list2[2], c='orange', 
            marker='.', ls='--', label=f"$E=0.0\\hbar\\omega$ eV")
    ax.set_xlabel('Time (Periods)')
    ax.set_ylabel('$n(t)$')
    ax.set_title(fr'Occupation')
    fig.suptitle(fig_title_info)
    ax.legend()
    """
    """
    # --------------------------------------------------------------------------
    # FREQUENCY ANALYSIS OF THE RESULTS
    # General figure 0.5HW
    max_freq = min([3, freq1[-1]*T/(2*np.pi), freq2[-1]*T/(2*np.pi)])
    props = dict(boxstyle='round', facecolor='white', edgecolor='grey', alpha=0.8)

    fig, ax = plt.subplots()
    ax.plot(freq1*T/(2*np.pi), fourier_occ1[3], c='blue', marker='.', ls='--', 
            label=f'$E = 0.5\\hbar\\omega$')
    ax.scatter(freq1[max_freq_ind1[3]]*T/(2*np.pi), fourier_occ1[3, max_freq_ind1[3]], color='cyan',
            marker='*', zorder=2)
    
    ax.plot(freq1*T/(2*np.pi), fourier_occ2[3], c='red', marker='.', ls='--', 
            label=f'$E = 0.5\\hbar\\omega$')
    # Markers of max frequencies
    ax.scatter(freq2[max_freq_ind2[3]]*T/(2*np.pi), fourier_occ2[3, max_freq_ind2[3]], color='darkred',
            marker='*', zorder=2)
    ax.vlines([range_search], 0, np.max(fourier_occ1), ls='-.', color='gray')
    ax.set_xlim(0, max_freq)
    wc_text = f'$\\omega_c = {char_freq1[g, 3]:.6f}$ fs$^{{-1}}$\n$\\omega_c = {char_freq2[g, 3]:.6f}$ fs$^{{-1}}$'
    ax.text(0.72, 0.98, wc_text, transform=ax.transAxes,
            verticalalignment='top', bbox=props)
    ax.legend(loc=(0.47, 0.7815), labelspacing=0.8, edgecolor='grey')
    ax.set_xlabel('Normal Frequency (period$^{-1}$)')
    ax.set_ylabel('Amplitude')
    ax.set_title('FFT')
    fig.suptitle(fig_title_info)
    """
    """THESE DONT WORK YET
    
    # General figure 1.0HW
    max_freq = min([3, freq1[-1]*T/(2*np.pi), freq2[-1]*T/(2*np.pi)])
    props = dict(boxstyle='round', facecolor='white', edgecolor='grey', alpha=0.8)

    fig, ax = plt.subplots()
    ax.plot(freq1*T/(2*np.pi), fourier_occ1[3], c='blue', marker='.', ls='--', 
            label=f'$E = 0.5\\hbar\\omega$')
    # Markers of max frequencies2)
    ax.scatter(freq1[max_freq_ind1[3]]*T/(2*np.pi), fourier_occ1[3, max_freq_ind1[3]], color='cyan',
            marker='*', zorder=2)
    ax.plot(freq2*T/(2*np.pi), fourier_occ2[3], c='blue', marker='.', ls='--', 
            label=f'$E = 0.5\\hbar\\omega$')
    # Markers of max frequencies
    ax.scatter(freq2[max_freq_ind2[3]]*T/(2*np.pi), fourier_occ2[3, max_freq_ind2[3]], color='cyan',
            marker='*', zorder=2)
    ax.vlines([range_search], 0, np.max(fourier_occ1), ls='-.', color='gray')
    ax.set_xlim(0, max_freq)
    wc_text = f'$\\omega_c = {char_freq1[g, 3]:.6f}$ fs$^{{-1}}$\n$\\omega_c = {char_freq2[g, 3]:.6f}$ fs$^{{-1}}$'
    ax.text(0.72, 0.98, wc_text, transform=ax.transAxes,
            verticalalignment='top', bbox=props)
    ax.legend(loc=(0.47, 0.7815), labelspacing=0.8, edgecolor='grey')
    ax.set_xlabel('Normal Frequency (period$^{-1}$)')
    ax.set_ylabel('Amplitude')
    ax.set_title('FFT')
    fig.suptitle(fig_title_info)

    # General figure 0.0HW
    max_freq = min([3, freq1[-1]*T/(2*np.pi), freq2[-1]*T/(2*np.pi)])
    props = dict(boxstyle='round', facecolor='white', edgecolor='grey', alpha=0.8)

    fig, ax = plt.subplots()
    ax.plot(freq1*T/(2*np.pi), fourier_occ1[3], c='blue', marker='.', ls='--', 
            label=f'$E = 0.5\\hbar\\omega$')
    # Markers of max frequencies2)
    ax.scatter(freq1[max_freq_ind1[3]]*T/(2*np.pi), fourier_occ1[3, max_freq_ind1[3]], color='cyan',
            marker='*', zorder=2)
    ax.plot(freq2*T/(2*np.pi), fourier_occ2[3], c='blue', marker='.', ls='--', 
            label=f'$E = 0.5\\hbar\\omega$')
    # Markers of max frequencies
    ax.scatter(freq2[max_freq_ind2[3]]*T/(2*np.pi), fourier_occ2[3, max_freq_ind2[3]], color='cyan',
            marker='*', zorder=2)
    ax.vlines([range_search], 0, np.max(fourier_occ1), ls='-.', color='gray')
    ax.set_xlim(0, max_freq)
    wc_text = f'$\\omega_c = {char_freq1[g, 3]:.6f}$ fs$^{{-1}}$\n$\\omega_c = {char_freq2[g, 3]:.6f}$ fs$^{{-1}}$'
    ax.text(0.72, 0.98, wc_text, transform=ax.transAxes,
            verticalalignment='top', bbox=props)
    ax.legend(loc=(0.47, 0.7815), labelspacing=0.8, edgecolor='grey')
    ax.set_xlabel('Normal Frequency (period$^{-1}$)')
    ax.set_ylabel('Amplitude')
    ax.set_title('FFT')
    fig.suptitle(fig_title_info)
    """

    # --------------------------------------------------------------------------
    # FREQUENCY ANALYSIS OF THE RESULTS
    # General figure 0.5HW
    omega = np.sqrt(3)*2.7*gamma*np.pi/jcl.hbar_fs
    fig, ax = plt.subplots()
    ax.plot(freq1/omega, fourier_occ1[3], c='blue', marker='.', ls='--', 
            label=f'{label1}')
    ax.scatter(freq1[max_freq_ind1[3]]/omega, fourier_occ1[3, max_freq_ind1[3]], color='cyan',
            marker='*', zorder=2)
    
    ax.plot(freq2/omega, fourier_occ2[3], c='red', marker='.', ls='--', 
            label=f'{label2}')
    # Markers of max frequencies
    ax.scatter(freq2[max_freq_ind2[3]]/omega, fourier_occ2[3, max_freq_ind2[3]], color='darkred',
            marker='*', zorder=2)
    ax.vlines(1, 0, np.max(fourier_occ1), ls='-.', color='gray')
    ax.set_xlabel('Angular freq. $/ \\Omega$')
    ax.set_ylabel('Amplitude')
    ax.set_title('FFT on $E=0.5\\hbar\\omega$')
    ax.legend()
    fig.suptitle(fig_title_info)
# Final plot of the gammas
fig, ax = plt.subplots()
#ax.plot(gamma_list, char_freq[:,2], ls='--', c=color_list[2], marker='.',
#        label=f'$E = \\mu$ eV')
ax.plot(gamma_list, char_freq1[:,3], ls='--', c=color_list[3], marker='.',
        label=f'{label1}')
ax.plot(gamma_list, char_freq2[:,3], ls='--', c=color_list[4], marker='.',
        label=f'{label2}')
ax.set_xlabel(f'Intensity $\\Gamma$')
ax.set_ylabel('Angular frequency (las. period$^{-1}$)')
ax.set_title(f'$\\omega_c$ in $E= 0.5\\hbar\\omega$')
ax.legend()
ax.set_xticks(gamma_list)
fig.suptitle(f'$N={{{2**N_pot}}}$, $\\hbar\\omega$={E} eV, Temp={Temp} K, $\\mu$={mu} eV')
