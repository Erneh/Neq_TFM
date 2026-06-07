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
from core import DOS_sparse, load_data, load_data_dict, frequency_analysis_dict
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
gamma_list = np.linspace(0.000, 0.050, 11)                   

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


hE_list = [hE*E/2 for hE in range(-hE_reps, hE_reps+1)]
color_list = ['olivedrab', 'red', 'orange', 'blue', 'darkviolet']



# ------------------------------------------------------------------------------
## CREATING A DICTIONARY FOR EACH COMPARATION TO MAKE
res1 = {'modifier_id': 'linear',
        'type_ham': 'basic',
        'ham_params': 0.50,  
        'hw' : 1.0,                       
        'Temp' :  1e-9,
        'mu' : 0.01,                 
        'N_pot' : 19,
        'n_periods' : 500,
        'steps_per_T' : 1000,
        'meas_per_T' : 16,
        'N_random_vector' : 5,
        'M' : int(np.sqrt(2**19)),
        'label': 'linear graphene',
        'color1': 'blue',
        'color2': 'cyan',
        'ls' : 'dashed'}


res2 = {'modifier_id': 'circle',
        'type_ham': 'basic',
        'ham_params': 0.50,  
        'hw' : 1.0,                       
        'Temp' :  1e-9,
        'mu' : 0.01,                 
        'N_pot' : 19,
        'n_periods' : 500,
        'steps_per_T' : 1000,
        'meas_per_T' : 16,
        'N_random_vector' : 5,
        'M' : int(np.sqrt(2**19)),
        'label': 'circle graphene',
        'color1': 'red',
        'color2': 'magenta',
        'ls' : 'dotted'}

res3 = {'modifier_id': 'linear',
        'type_ham': 'hbn',
        'ham_params': 0.50,  
        'hw' : 1.1,                       
        'Temp' :  1e-9,
        'mu' : 0.01,                 
        'N_pot' : 19,
        'n_periods' : 500,
        'steps_per_T' : 1000,
        'meas_per_T' : 16,
        'N_random_vector' : 5,
        'M' : int(np.sqrt(2**19)),
        'label': 'circle graphene',
        'color1': 'blue',
        'color2': 'cyan',
        'ls' : 'dashed'}

res4 = {'modifier_id': 'circle',
        'type_ham': 'hbn',
        'ham_params': 0.50,  
        'hw' : 1.1,                       
        'Temp' :  1e-9,
        'mu' : 0.01,                 
        'N_pot' : 19,
        'n_periods' : 500,
        'steps_per_T' : 1000,
        'meas_per_T' : 16,
        'N_random_vector' : 5,
        'M' : int(np.sqrt(2**19)),
        'label': 'circle graphene',
        'color1': 'red',
        'color2': 'magenta',
        'ls' : 'dotted'}

list_res = [res1, res2, res3, res4]
# Initialize arrays to save results and add auxiliary stuff
for res in list_res:
    res.update({'char_freq' : np.zeros((len(gamma_list), len(hE_list)))})
    res['w'] = res['hw']/jcl.hbar_fs 
    res['T'] = 2*np.pi/res['w']    
    res['t_vec'] = np.linspace(0,n_periods*res['T'] , res['steps_per_T']*res['n_periods'])   
    res['t_vec_measures'] = np.linspace(0, n_periods*res['T'], res['meas_per_T']*res['n_periods'])

#%% True Loop Start
for (g, gamma) in enumerate(gamma_list):
    fig_title_info = f'N={2**N_pot}, $\\hbar \\omega$={E} eV, T={Temp} K, $\\mu$={mu} eV, $\\Gamma$={gamma}, M={M}'
    for res in list_res:
        res['gamma'] = gamma
        load_data_dict(res)
        frequency_analysis_dict(res)
        res['char_freq'][g] = res['char_freq_s']
    
    
    # --------------------------------------------------------------------------
    # OCCUPATION(time) 0.5HW
    # General figure to contain all important data
    fig, ax = plt.subplots()
    for res in list_res:
        ax.plot(res['t_vec_measures']/res['T'], res['occ_drop_list'][3], c=res['color1'], 
                marker='.', ls=res['ls'], label=res['label'])
    ax.set_xlabel('Time (Periods)')
    ax.set_ylabel('$n(t)$')
    ax.set_title(f'n(t) in E=$0.5\\hbar\\omega$')
    fig.suptitle(fig_title_info)
    ax.legend()
    
    # --------------------------------------------------------------------------
    # OCCUPATION(time)  1.0HW
    # General figure to contain all important data
    fig, ax = plt.subplots()
    for res in list_res:
        ax.plot(res['t_vec_measures']/res['T'], res['occ_drop_list'][4], c=res['color1'], 
                marker='.', ls=res['ls'], label=res['label'])
    ax.set_xlabel('Time (Periods)')
    ax.set_ylabel('$n(t)$')
    ax.set_title(f'DOS*n(t) in E=$1.0\\hbar\\omega$')
    fig.suptitle(fig_title_info)
    ax.legend()

    # --------------------------------------------------------------------------
    # OCCUPATION(time)  0.0HW
    # General figure to contain all important data
    fig, ax = plt.subplots()
    for res in list_res:
        ax.plot(res['t_vec_measures']/res['T'], res['occ_drop_list'][2], c=res['color1'], 
                marker='.', ls=res['ls'], label=res['label'])
    ax.set_xlabel('Time (Periods)')
    ax.set_ylabel('$n(t)$')
    ax.set_title(f'DOS*n(t) in E=$0.0\\hbar\\omega$')
    fig.suptitle(fig_title_info)
    ax.legend()
    
    
    # --------------------------------------------------------------------------
    # FREQUENCY ANALYSIS OF THE RESULTS
    # General figure 0.5HW
    if gamma != 0:
        props = dict(boxstyle='round', facecolor='white', edgecolor='grey', alpha=0.8)
        omega = np.sqrt(3)*2.7*gamma*np.pi/jcl.hbar_fs
        
 
        fig, ax = plt.subplots()
        for res in list_res:
            ax.plot(res['freq']/omega, res['fourier_occ'][3], c=res['color1'], marker='.', ls=res['ls'], label=res['label'])
            ax.scatter(res['freq'][res['max_freq_ind'][3]]/omega, res['fourier_occ'][3, res['max_freq_ind'][3]], 
                    color=res['color2'], marker='*', zorder=2)
            
        x_min, x_max = ax.get_xlim()
        y_min, y_max = ax.get_ylim()
        ax.set_ylim(y_min, y_max)
        ax.vlines(1, 0, y_max, ls='-', color='gray', zorder=1)
    
        #ax.set_xlim(0, max_freq/omega)
        #ax.vlines(np.arange(w, 4*w, w)/omega, 0, np.max(fourier_occ1), ls=(0, (1, 2)), color='darkslategrey')
        ax.set_xlabel('Angular freq. $/ \\Omega$')
        ax.set_ylabel('Amplitude')
        ax.set_title('FFT on $E=0.5\\hbar\\omega$')
        #ax.text(w/omega, np.max(fourier_occ1)*6/8, '$\\omega = \\omega_L$', verticalalignment='center', horizontalalignment='center', bbox=props)
        ax.legend()
        fig.suptitle(fig_title_info)

# Final plot of the gammas
fig, ax = plt.subplots()
for res in list_res:
    ax.plot(gamma_list, res['char_freq'][:,3],ls=res['ls'], c=res['color1'], marker='.',
            label=res['label'])

ax.set_xlabel(f'Intensity $\\Gamma$')
ax.set_ylabel('Angular frequency (las. period$^{-1}$)')
ax.set_title(f'$\\omega_c$ in $E= 0.5\\hbar\\omega$')
ax.legend()
ax.set_xticks(gamma_list)
fig.suptitle(f'$N={{{2**N_pot}}}$, $\\hbar\\omega$={E} eV, Temp={Temp} K, $\\mu$={mu} eV')


#%%
# Comparison for the presentation
char_freqs_og = np.load('Out/char_freqs.npy')[:len(gamma_list)]
omega_list = np.sqrt(3)*2.7*gamma_list*np.pi/jcl.hbar_fs
fig, ax = plt.subplots()
ax.plot(gamma_list, char_freq1[:,3], ls='--', c=color1, marker='.',
        label=f'hBN {label1}')
ax.plot(gamma_list, char_freq2[:,3], ls='--', c=color2, marker='.',
        label=f'hBN {label2}')
ax.plot(gamma_list, char_freqs_og, ls='--', c='green', marker='.',
        label=f'graphene circle')
ax.plot(gamma_list, lin_graphene, ls='--', c='magenta', marker='.',
        label=f'graphene linear')
ax.plot(gamma_list, omega_list, c='gray', label='$\\Omega(\\Gamma)$')
ax.set_xlabel(f'Intensity $\\Gamma$')
ax.set_ylabel('Angular frequency (las. period$^{-1}$)')
ax.set_title(f'$\\omega_c$ in $E= 0.5\\hbar\\omega$')
ax.legend()
ax.set_xticks(gamma_list)
fig.suptitle(f'$N={{{2**N_pot}}}$, $\\hbar\\omega$={E} eV, Temp={Temp} K, $\\mu$={mu} eV')
