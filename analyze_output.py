<<<<<<< HEAD
"""
Script made to analyze the output of the calculations. Will fail if the given
calcs do not exist!
"""
#%% First part: imports and analysis of just one specific case
=======
>>>>>>> Massive amount of changes and reorganizing of code. New results
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
<<<<<<< HEAD
from core import DOS_sparse, frequency_analysis
=======
from core import DOS_sparse
>>>>>>> Massive amount of changes and reorganizing of code. New results
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
# Benchmarking
from time import time
from datetime import timedelta

<<<<<<< HEAD

# ------------------------------------------------------------------------------
####  Defining a hamiltonian (own)
N_pot = 17
N = 2**N_pot
=======
#%% NEQ Calcs
# ------------------------------------------------------------------------------
####  Defining a hamiltonian (own)
N = 2**17
>>>>>>> Massive amount of changes and reorganizing of code. New results
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
<<<<<<< HEAD
mu = 0.40
# Amount of periods to be simulated
n_periods = 100
=======
mu = 0.01
# Amount of periods to be simulated
n_periods = 50
>>>>>>> Massive amount of changes and reorganizing of code. New results
# Amount of half multiples of E where the occupation is obtained
hE_reps = 2
# ??
tau = 0.0 
# Amount of measures per period
<<<<<<< HEAD
meas_per_T = 8
=======
meas_per_T = 32
>>>>>>> Massive amount of changes and reorganizing of code. New results
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
<<<<<<< HEAD
gamma = 0.010        
=======
gamma = 0.010                   
>>>>>>> Massive amount of changes and reorganizing of code. New results

# Intensity of the laser
Phi0 = jcl.hbar_fs*2*np.pi
A0 = gamma*Phi0/(2*3**0.5*jcl.a_cc)

# Amount of random vectors used in calculation
<<<<<<< HEAD
N_random_vector = 5

# Momenta
M = int(np.sqrt(N))
#M = 362
=======
N_random_vector = 10

# Momenta
M = int(np.sqrt(N))
>>>>>>> Massive amount of changes and reorganizing of code. New results
# Broadening in the energies
broad = dE*np.pi/M

t_vec_measures = np.linspace(0, n_periods*T, N_measures)
<<<<<<< HEAD

# Range of searching the maximim frequency (in period^-1 units)
range_search = 1


=======
>>>>>>> Massive amount of changes and reorganizing of code. New results
# Print of parameters to check results
print('PARAMETERS OF ANALYSIS')
print()
print(f'Type of light is {modifier_id}')
print(f'# of atoms: {N}')
print(f'Energy: {E} eV')
print(f'Intensity param: {gamma}')
print(f'Temperature: {Temp} K')
print(f'Chem potential: {mu} eV')
print(f'# of periods: {n_periods}')
print(f'steps/period: {steps_per_T}')
print(f'# measures/T: {meas_per_T}')


# Names of file and info on graphs
<<<<<<< HEAD
folder_name = f'{modifier_id}/G={gamma:.3f}_E={float(E)}_Temp={Temp}_mu={mu:.2f}/N={N_pot}_M={M}_R={N_random_vector}_nT={n_periods}_measT={meas_per_T}_stT={steps_per_T}'
fig_title_info = f'$N={{{2**N_pot}}}$, E={E} eV, T={Temp} K, $\\mu$={mu} eV, $\\Gamma$={gamma}, M={M}'
pulse_suptitle = fr'$E={E}, \omega = {w:.3f}$ fs$^{{-1}}, T={T:.3f}$ fs'
extra_text = f'Light: {modifier_id}\n$\\delta E={broad:.3f}$\n# Rand Vecs: {N_random_vector}\nMeasures/T={meas_per_T}\nSteps/T={steps_per_T}'

# Loading the info in the .npy files available at the moment
if N_random_vector == 1:
    EF_list = np.load(f'Out/{folder_name}/E.npy')
    n_E_list = np.load(f'Out/{folder_name}/n_E.npy')
    dos_list = np.load(f'Out/{folder_name}/dos_E.npy')
else:
    # See the amount calculated in the corresponding folder
    R_calc = len(os.listdir(f'Out/{folder_name}/Ene_R'))
    print(f'{R_calc}/{N_random_vector} calculations finished! Showing results...')
    EF_list = np.load(f'Out/{folder_name}/Ene_R/1.npy')
    n_E_list = np.zeros_like(EF_list)
    dos_list = np.zeros_like(EF_list)
    for i in range(1, R_calc+1):
        n_E_list += np.load(f'Out/{folder_name}/noc_R/{i}.npy')
        dos_list += np.load(f'Out/{folder_name}/dos_R/{i}.npy')
    n_E_list /= R_calc
    dos_list /= R_calc
            
# Remaking the graphs of oscillations in time
hE_list = [hE*E/2 for hE in range(-hE_reps, hE_reps+1)]


# Energy window and lines
cmap = plt.cm.plasma 
norm = Normalize(vmin=t_vec.min(), vmax=t_vec.max())
min_e, max_e = -2.5, 2.5
hw_lines_step = 0.5
hw_hlines = [i*E for i in [hw_lines_step*i for i in range(1, int(max_e/hw_lines_step+1))]]
hw_hlines += [i*E for i in [-hw_lines_step*i for i in range(1, int(max_e/hw_lines_step+1))]]

# Text box
props = dict(boxstyle='round', facecolor='white', alpha=1.0)

fig, ax = plt.subplots()
for i in range(N_measures):
    ax.plot(EF_list[i,:], n_E_list[i,:], color=cmap(norm(t_vec_measures[i])), label=f't={round(t_vec_measures[i]/T, 3)}T')
#ax.legend()
cbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax, orientation='vertical', label='Time (periods)')
ax.set_xlabel('Energy (eV)')
ax.set_ylabel('$n(\\varepsilon)$')
ax.set_xlim(min_e, max_e)
ax.vlines(hw_hlines, 0, 1, color='grey', ls='--', alpha=0.5, zorder=1)
fig.suptitle(fig_title_info)
ax.set_title('Occupation Number')
ax.text(min_e + 0.5, 0.2, extra_text, bbox=props)
fig.savefig(f'Out/{folder_name}/N(E).png', bbox_inches='tight')


# Density of states graph
fig, ax = plt.subplots()
for i in range(N_measures):
    ax.plot(EF_list[i,:], dos_list[i,:], color=cmap(norm(t_vec_measures[i])), label=f't={round(t_vec_measures[i]/T, 3)}T')
#ax.legend()
cbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax, orientation='vertical', label='Time (periods)')
ax.set_xlabel('Energy (eV)')
ax.set_ylabel('$n(\\varepsilon)$')
#ax.vlines(hw_hlines, 0, np.max(dos_list), color='grey', ls='--', alpha=0.5, zorder=1)
ax.vlines([0], 0, np.max(dos_list), color='grey', ls='--', alpha=0.8, zorder=1)
fig.suptitle(fig_title_info)
ax.set_title('Density of States')
fig.savefig(f'Out/{folder_name}/dos(E).png', bbox_inches='tight')

# Extra graphs for the n for different times
hE_list = [hE*E/2 for hE in range(-hE_reps, hE_reps+1)]
color_list = ['olivedrab', 'red', 'orange', 'blue', 'darkviolet']

occ_drop_list, fourier_occ, freq, char_freq, max_freq_ind = frequency_analysis(EF_list, n_E_list, hE_list, t_vec_measures, T, range_search)
# General figure to contain all important data
fig, ax = plt.subplots()
reescale = np.max(occ_drop_list[3]) / np.max(occ_drop_list[4]) / 2
ax.plot(np.array(t_vec_measures)/T, occ_drop_list[4]*reescale, c='darkviolet', 
        marker='.', ls='--', label=f"$E=1E+\\mu$ eV $\\cdot$ {reescale:.3f}")
ax.plot(np.array(t_vec_measures)/T, occ_drop_list[3], c='blue', 
        marker='.', ls='--', label=f"$E=0.5E+\\mu$ eV")
ax.plot(np.array(t_vec_measures)/T, occ_drop_list[2], c='orange', 
        marker='.', ls='--', label=f"$E=\\mu$ eV")
ax.set_xlabel('Time (Periods)')
ax.set_ylabel('$n(t)$')
ax.set_title(fr'Occupation')
fig.suptitle(fig_title_info)
ax.legend()
fig.savefig(f'Out/{folder_name}/N(T).png', bbox_inches='tight')

# ------------------------------------------------------------------------------
# FREQUENCY ANALYSIS OF THE RESULTS
# General figure
max_freq = min(20, freq[-1]*T/(2*np.pi))

fig, ax = plt.subplots()
ax.plot(freq*T/(2*np.pi), fourier_occ[2], c='orange', marker='.', ls='--', 
        label=f'$E = \\mu$ eV, $f_c={char_freq[2]:.6f}$')
ax.plot(freq*T/(2*np.pi), fourier_occ[3], c='blue', marker='.', ls='--', 
        label=f'$E = 0.5E + \\mu$ eV, $f_c={char_freq[3]:.6f}$')
ax.plot(freq*T/(2*np.pi), fourier_occ[4], c='darkviolet', marker='.', ls='--', 
        label=f'$E = 1E+\\mu$ eV, $f_c={char_freq[4]:.6f}$')
# Markers of max frequencies
ax.scatter(freq[max_freq_ind[2]]*T/(2*np.pi), fourier_occ[2, max_freq_ind[2]], color='red',
           marker='*', zorder=2)
ax.scatter(freq[max_freq_ind[3]]*T/(2*np.pi), fourier_occ[3, max_freq_ind[3]], color='cyan',
           marker='*', zorder=2)
ax.scatter(freq[max_freq_ind[4]]*T/(2*np.pi), fourier_occ[4, max_freq_ind[4]], color='magenta',
           marker='*', zorder=2)
ax.set_xlabel('Normal Frequency (period$^{-1}$)')
ax.set_ylabel('Amplitude')
ax.set_title('FFT')
fig.suptitle(fig_title_info)
ax.legend()
ax.set_xlim(0, max_freq)
ax.vlines([range_search], 0, np.max(fourier_occ), ls='-.', color='gray')
fig.savefig(f'Out/{folder_name}/FREQ_N(T).png', bbox_inches='tight')

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
modifier_id = 'linear'
# Amount of random vectors used in calculation
N_random_vector = 1
range_search = 1
# Intensity param     (no units)
gamma_list = np.linspace(0.005, 0.050, 10)                   
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

char_freq = np.zeros((len(gamma_list), len(hE_list)))
for (g, gamma) in enumerate(gamma_list):
    if M == 0:
        M = int(np.sqrt(N))
    

    # Names of file and info on graphs
    folder_name = f'{modifier_id}/G={gamma:.3f}_E={float(E)}_Temp={Temp}_mu={mu:.2f}/N={N_pot}_M={M}_R={N_random_vector}_nT={n_periods}_measT={meas_per_T}_stT={steps_per_T}'
    fig_title_info = f'$N={{{2**N_pot}}}$, E={E} eV, T={Temp} K, $\\mu$={mu} eV, $\\Gamma$={gamma}, M={M}'
    pulse_suptitle = fr'$E_0={E}, \omega = {w:.3f}$ fs$^{{-1}}, T={T:.3f}$ fs'

    # Loading the info in the .npy files
    EF_list = np.load(f'Out/{folder_name}/E.npy')
    n_E_list = np.load(f'Out/{folder_name}/n_E.npy')

    dE = EF_list[0, 1] - EF_list[0, 0]
    broad = dE*np.pi/M
    extra_text = f'Light: {modifier_id}\n$\\delta E={broad:.3f}$\n# Rand Vecs: {N_random_vector}\nMeasures/T={meas_per_T}\nSteps/T={steps_per_T}'

    # Getting the characteristic period of each system
    occ_drop_list, fourier_occ, freq, char_freq[g], max_freq_ind = frequency_analysis(EF_list, n_E_list, hE_list, t_vec_measures, T, range_search)
    # --------------------------------------------------------------------------
    # OCCUPATION(time) GRAPH
    # General figure to contain all important data
    fig, ax = plt.subplots()
    reescale = np.max(occ_drop_list[3]) / np.max(occ_drop_list[4]) / 2
    ax.plot(np.array(t_vec_measures)/T, occ_drop_list[4]*reescale, c='darkviolet', 
            marker='.', ls='--', label=f"$E=1E+\\mu$ eV $\\cdot$ {reescale:.3f}")
    ax.plot(np.array(t_vec_measures)/T, occ_drop_list[3], c='blue', 
            marker='.', ls='--', label=f"$E=0.5E+\\mu$ eV")
    ax.plot(np.array(t_vec_measures)/T, occ_drop_list[2], c='orange', 
            marker='.', ls='--', label=f"$E=\\mu$ eV")
    ax.set_xlabel('Time (Periods)')
    ax.set_ylabel('$n(t)$')
    ax.set_title(fr'Occupation')
    fig.suptitle(fig_title_info)
    ax.legend()
    fig.savefig(f'Out/{folder_name}/N(T).png', bbox_inches='tight')
    # --------------------------------------------------------------------------
    # FREQUENCY ANALYSIS OF THE RESULTS
    # General figure
    max_freq = min(20, freq[-1]*T/(2*np.pi))

    fig, ax = plt.subplots()
    ax.plot(freq*T/(2*np.pi), fourier_occ[2], c='orange', marker='.', ls='--', 
        label=f'$E = \\mu$ eV, $f_c={char_freq[g, 2]:.6f}$')
    ax.plot(freq*T/(2*np.pi), fourier_occ[3], c='blue', marker='.', ls='--',
        label=f'$E = 0.5E + \\mu$ eV, $f_c={char_freq[g, 3]:.6f}$')
    ax.plot(freq*T/(2*np.pi), fourier_occ[4], c='darkviolet', marker='.', ls='--', 
        label=f'$E = 1E+\\mu$ eV, $f_c={char_freq[g, 4]:.6f}$')
    # Markers of max frequencies
    ax.scatter(freq[max_freq_ind[2]]*T/(2*np.pi), fourier_occ[2, max_freq_ind[2]], color='red',
        marker='*', zorder=2)
    ax.scatter(freq[max_freq_ind[3]]*T/(2*np.pi), fourier_occ[3, max_freq_ind[3]], color='cyan',
        marker='*', zorder=2)
    ax.scatter(freq[max_freq_ind[4]]*T/(2*np.pi), fourier_occ[4, max_freq_ind[4]], color='magenta',
        marker='*', zorder=2)
    ax.set_xlabel('Normal Frequency (period$^{-1}$)')
    ax.set_ylabel('Amplitude')
    ax.set_title('FFT')
    fig.suptitle(fig_title_info)
    ax.legend()
    ax.set_xlim(0, max_freq)
    ax.vlines([range_search], 0, np.max(fourier_occ), ls='-.', color='gray')
    fig.savefig(f'Out/{folder_name}/FREQ_N(T).png', bbox_inches='tight')
    
    

# Final plot of the gammas
fig, ax = plt.subplots()
#ax.plot(gamma_list, char_freq[:,2], ls='--', c=color_list[2], marker='.',
#        label=f'$E = \\mu$ eV')
ax.plot(gamma_list, char_freq[:,3], ls='--', c=color_list[3], marker='.',
        label=f'$E = 0.5E_0 + \\mu$ eV')
ax.plot(gamma_list, char_freq[:,4], ls='--', c=color_list[4], marker='.',
        label=f'$E = 1E_0 + \\mu$ eV')
ax.set_xlabel(f'Intensity $\\Gamma$')
ax.set_ylabel('Angular frequency (las. period$^{-1}$)')
ax.set_title(f'Characteristic frequencies for different intensities')
ax.legend()
ax.set_xticks(gamma_list)
fig.suptitle(f'$N=2^{{{N_pot}}}$, E_0={E} eV, Temp={Temp} K, $\\mu$={mu} eV')


#%% Fitting to see if a relation is even possible
# Importing the results from Floquet
imported_inds = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
model_select = 'real'
gamma_fl = np.load(f'/home/eperez/Code/Floquet_tfm/Outr/GAP/gamma_list.npy')[imported_inds]
gap0_fl = np.load(f'/home/eperez/Code/Floquet_tfm/Outr/GAP/gap0_{model_select}.npy')[imported_inds]
gap1_fl = np.load(f'/home/eperez/Code/Floquet_tfm/Outr/GAP/gap1_{model_select}.npy')[imported_inds]
gap2_fl = np.load(f'/home/eperez/Code/Floquet_tfm/Outr/GAP/gap2_{model_select}.npy')[imported_inds]

# Fitting of the curve
from scipy.stats import linregress
# GAP1 REGRESSION
lin_reg1 = linregress(gap1_fl, char_freq[:,3]*jcl.hbar_fs)
slope1 = lin_reg1.slope
intercept1 = lin_reg1.intercept
R2_1 = lin_reg1.rvalue**2
test_gaps1 = np.linspace(np.min(gap1_fl), np.max(gap1_fl), 100)
lin_predict1 = slope1*test_gaps1 + intercept1

# GAP2 REGRESSION
lin_reg2 = linregress(gap2_fl, char_freq[:,4]*jcl.hbar_fs)
slope2 = lin_reg2.slope
intercept2 = lin_reg2.intercept
R2_2 = lin_reg2.rvalue**2
test_gaps2 = np.linspace(np.min(gap2_fl), np.max(gap2_fl), 100)
lin_predict2 = slope2*test_gaps2 + intercept2

# Graph of the results
fig, ax = plt.subplots(dpi=200)
ax.plot(gap1_fl, char_freq[:,3]*jcl.hbar_fs, ls='--', marker='.', color='blue', label='$\\Delta_1$ points')
ax.plot(test_gaps1, lin_predict1, color='cyan', label=f'$y={slope1:.3f}x+{intercept1:.3f}, R^2={R2_1:.3f} $')

ax.plot(gap2_fl, char_freq[:,4]*jcl.hbar_fs, ls='--', marker='.', color='orange', label='$\\Delta_2$ points')
ax.plot(test_gaps2, lin_predict2, color='red', label=f'$y={slope2:.3f}x+{intercept2:.3f}, R^2={R2_2:.3f} $')
ax.set_xlabel('Gap size (eV)')
ax.set_ylabel('$\\hbar \\omega_{char}$ (eV)')
ax.set_title('Comparison Floquet vs NEQ')
ax.legend()

#%% Trying Rabi stuff to see if theoretical results make sense using gap
t = -2.7
w_laser = E/jcl.hbar_fs
w1 = np.sqrt(3)*t*gamma_fl*np.pi/jcl.hbar_fs/2
#w1 = 3*t**2*gamma_list**2*np.pi**2
rep = 1
delta_w = (1-rep)*w_laser
#delta_w = w_laser - gap1_fl/jcl.hbar_fs
rabi_freq = np.sqrt(delta_w**2 + w1**2)


fig, ax = plt.subplots(dpi=200)
ax.plot(gamma_list, char_freq[:,3], ls='--', c=color_list[3], marker='.',
        label=f'Numerical ang. freq')
ax.plot(gamma_fl, rabi_freq, ls='--', c='red', marker='.',
        label=f'Rabi ang. freq')
ax.set_title('Characteristic Frequencies in $E=0.5E_0$')
ax.set_xlabel('Intensity parameter $\\Gamma$')
ax.set_ylabel('Ang. frequency')
ax.legend()

=======
folder_name = f'circle/N={N}_E={float(E)}_Temp={Temp}_mu={mu:.2f}_G={gamma:.3f}/Nrand={N_random_vector}_nT={n_periods}_measT={meas_per_T}_stT={steps_per_T}'
fig_title_info = f'N={N}, E={E} eV, T={Temp} K, $\\mu$={mu} eV, $\\Gamma$={gamma}'
pulse_suptitle = fr'$E={E}, \omega = {w:.3f}$ fs$^{{-1}}, T={T:.3f}$ fs'
extra_text = f'Light: circle\n$\\delta E={broad:.3f}$\nMeasures/T={meas_per_T}\nSteps/T={steps_per_T}'

# Loading the info in the .npy files
EF_list = np.load(f'Out/{folder_name}/E.npy')
n_E_list = np.load(f'Out/{folder_name}/n_E.npy')

# Remaking the graphs of oscillations in time
hE_list = [hE*E/2 for hE in range(-hE_reps, hE_reps+1)]
total_nhE = 2*hE_reps + 1 
occ_drop_list = np.zeros((total_nhE, N_measures))
for (i, hE) in enumerate(hE_list):
    ind_time = np.where(EF_list[i,:] > hE + mu)[0][0]
    occ_drop_list[i] = [n_E_list[i,ind_time] for i in range(N_measures)]

    fig, ax = plt.subplots()
    ax.plot(np.array(t_vec_measures)/T, occ_drop_list[i], c='blue', marker='.', ls='--')
    ax.set_xlabel('Time (Periods)')
    ax.set_ylabel('$n(\\varepsilon)$')
    ax.set_title(fr'Occupation in $E={hE/E:.1f}E$ + $\mu$')
    fig.suptitle(fig_title_info)
    #fig.savefig(f'Out/{folder_name}/N(T)_{hE:.1f}.png', bbox_inches='tight')

# ------------------------------------------------------------------------------
# FREQUENCY ANALYSIS OF THE RESULTS
fourier_occ = np.abs(np.fft.rfft(occ_drop_list))
dt = t_vec_measures[1] - t_vec_measures[0]
df = 1/dt/N_measures
freq = np.arange(0, fourier_occ.shape[1], 1)*df
for (i, hE) in enumerate(hE_list):
    max_freq_ind = np.where(fourier_occ[i] == max(fourier_occ[i]))[0]

    fig, ax = plt.subplots()
    ax.set_ylim(0, 10)
    ax.plot(freq*T, fourier_occ[i], c='blue', marker='.', ls='--')
    ax.set_xlabel('Frequency (period$^{-1}$)')
    ax.set_ylabel('Amplitude')
    ax.set_title(fr'FFT in $E={hE/E:.1f}E$ + $\mu$')
    fig.suptitle(fig_title_info)
    #fig.savefig(f'Out/{folder_name}/N(T)_{hE:.1f}.png', bbox_inches='tight')
>>>>>>> Massive amount of changes and reorganizing of code. New results
