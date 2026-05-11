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
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
# Benchmarking
from time import time
from datetime import timedelta
#%% Other important functions

<<<<<<< HEAD
def check_if_calculated(modifier_id, N_pot, E, Temp, mu, gamma, M, N_random_vector,
                        n_periods, meas_per_T, steps_per_T):
    N = 2**N_pot
    # Names of file and info on graphs
    folder_name = f'{modifier_id}/G={gamma:.3f}_E={float(E)}_Temp={Temp}_mu={mu:.2f}/N={N_pot}_M={M}_R={N_random_vector}_nT={n_periods}_measT={meas_per_T}_stT={steps_per_T}'
=======
def check_if_calculated(modifier_id, N_pot, E, Temp, mu, gamma, N_random_vector,
                        n_periods, meas_per_T, steps_per_T):
    N = 2**N_pot
    # Names of file and info on graphs
    folder_name = f'{modifier_id}/N={N}_E={float(E)}_Temp={float(Temp)}_mu={mu:.2f}_G={gamma:.3f}/Nrand={N_random_vector}_nT={n_periods}_measT={meas_per_T}_stT={steps_per_T}'
>>>>>>> Massive amount of changes and reorganizing of code. New results
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
<<<<<<< HEAD
            np.load(f'Out/{folder_name}/n_E.npy')
            np.load(f'Out/{folder_name}/dosn_E.npy')
            np.load(f'Out/{folder_name}/dos_E.npy')
=======
>>>>>>> Massive amount of changes and reorganizing of code. New results
            if np.sum(np.abs(EF_list)) == 0:
                flag = False
        except FileNotFoundError:
            flag = False
    return flag

<<<<<<< HEAD
#%% NEQ Calcs in Full
# ------------------------------------------------------------------------------
####  Defining a hamiltonian (own)
def neq_sim(modifier_id, N_pot, E, Temp, mu, gamma, M, N_random_vector, 
=======
#%% NEQ Calcs
# ------------------------------------------------------------------------------
####  Defining a hamiltonian (own)
def neq_sim(modifier_id, N_pot, E, Temp, mu, gamma, N_random_vector, 
>>>>>>> Massive amount of changes and reorganizing of code. New results
            n_periods, meas_per_T, steps_per_T, force_recalc=False, show_figs=True):
    N = 2**N_pot
    N1 = N2 = int(np.sqrt(N))//2
    S = get_positions_graphene(N1, N2)

    Ham = create_graphene_ham(S, N1, N2, out_format='ELL')

    dE = (Ham.bounds[1] - Ham.bounds[0])/2
    # ------------------------------------------------------------------------------
    #### PARAMETERS   
    hE_reps = 2
    # ??
    tau = 0.0 

    N_measures = meas_per_T*n_periods
    # Parameters of the laser
    w = E/jcl.hbar_fs 
    T = 2*np.pi/w                            
    # Time of  (fs)

    t_vec = np.linspace(0,n_periods*T , steps_per_T*n_periods)      
                
    # Intensity of the laser
<<<<<<< HEAD
    #Phi0 = jcl.hbar_fs*2*np.pi
    #A0 = gamma*Phi0/(2*3**0.5*jcl.a_cc)
    A0 = np.pi*gamma/(3**0.5*jcl.a_cc)
=======
    Phi0 = jcl.hbar_fs*2*np.pi
    A0 = gamma*Phi0/(2*3**0.5*jcl.a_cc)


    # Momenta
    M = int(np.sqrt(N))
>>>>>>> Massive amount of changes and reorganizing of code. New results
    # Broadening in the energies
    broad = dE*np.pi/M

    # Different possibilities depending on type of light
    if modifier_id == 'circle':
        # Polarization (right or left)
        pol = 'r'
        modifier_params = (A0, w, pol)

    elif modifier_id == 'linear':
        modifier_params = (A0, w)

    elif modifier_id == 'linear_packed':
        Tp = T
        modifier_params = (A0, w, Tp)

    obs_list = [['n', N_measures, M]]


    # Print of parameters to check results
    print('PARAMETERS OF CALCULATION')
    print()
    print(f'Type of light is {modifier_id}')
    print(f'# of atoms: {N}')
<<<<<<< HEAD
    print(f'# of moments: {M}')
=======
>>>>>>> Massive amount of changes and reorganizing of code. New results
    print(f'Energy: {E} eV')
    print(f'Intensity param: {gamma}')
    print(f'Temperature: {Temp} K')
    print(f'Chem potential: {mu} eV')
    print(f'# of Random Vectors: {N_random_vector}')
    print(f'# of periods: {n_periods}')
    print(f'steps/period: {steps_per_T}')
    print(f'# measures/T: {meas_per_T}')

<<<<<<< HEAD
    # Say approximate calculation time
    aprox_time = 265/(2**17*362*20*1000*4)*(N*M*n_periods*steps_per_T*meas_per_T)
    print()
    print(f'Approximate time: {aprox_time:.5f}')
=======
>>>>>>> Massive amount of changes and reorganizing of code. New results
    # Set up cmap and norm for the graphs in the future
    cmap = plt.cm.plasma 
    norm = Normalize(vmin=t_vec.min(), vmax=t_vec.max())

    # Times where the measurements take place
    t_vec_measures = np.linspace(0, n_periods*T, N_measures)
    
<<<<<<< HEAD
    # Names of file and info on 
    folder_name = f'{modifier_id}/G={gamma:.3f}_E={float(E)}_Temp={Temp}_mu={mu:.2f}/N={N_pot}_M={M}_R={N_random_vector}_nT={n_periods}_measT={meas_per_T}_stT={steps_per_T}'
    fig_title_info = f'$N={{{2**N_pot}}}$, E={E} eV, T={Temp} K, $\\mu$={mu} eV, $\\Gamma$={gamma}, M={M}, R={N_random_vector}'
    pulse_suptitle = fr'$E={E}, \omega = {w:.3f}$ fs$^{{-1}}, T={T:.3f}$ fs'
    extra_text = f'Light: {modifier_id}\n$\\delta E={broad:.3f}$\n# Rand Vecs: {N_random_vector}\nMeasures/T={meas_per_T}\nSteps/T={steps_per_T}'

    already_calc = check_if_calculated(modifier_id, N_pot, E, Temp, mu, gamma, M, N_random_vector,
=======
    # Names of file and info on graphs
    folder_name = f'{modifier_id}/N={N}_E={float(E)}_Temp={Temp}_mu={mu:.2f}_G={gamma:.3f}/Nrand={N_random_vector}_nT={n_periods}_measT={meas_per_T}_stT={steps_per_T}'
    fig_title_info = f'N={N}, E={E} eV, T={Temp} K, $\\mu$={mu} eV, $\\Gamma$={gamma}'
    pulse_suptitle = fr'$E={E}, \omega = {w:.3f}$ fs$^{{-1}}, T={T:.3f}$ fs'
    extra_text = f'Light: {modifier_id}\n$\\delta E={broad:.3f}$\n# Rand Vecs: {N_random_vector}\nMeasures/T={meas_per_T}\nSteps/T={steps_per_T}'

    already_calc = check_if_calculated(modifier_id, N_pot, E, Temp, mu, gamma, N_random_vector,
>>>>>>> Massive amount of changes and reorganizing of code. New results
                            n_periods, meas_per_T, steps_per_T)
    perform_calc = force_recalc or (not already_calc)

    if perform_calc:
<<<<<<< HEAD
        print('Results not found. Calculation starts!')
        # Create said folder
        os.makedirs(f'Out/{folder_name}', exist_ok=True)
        if N_random_vector > 1:
            os.makedirs(f'Out/{folder_name}/Ene_R', exist_ok=True)
            os.makedirs(f'Out/{folder_name}/noc_R', exist_ok=True)
            os.makedirs(f'Out/{folder_name}/dos_R', exist_ok=True)
            os.makedirs(f'Out/{folder_name}/dosn_R', exist_ok=True)

        # -----------------------------------------------------------------------------
        ## PRE-CALL: calculating the pulse in the given timeframe
=======
        print('Calculation starts!')
        # Create said folder
        os.makedirs(f'Out/{folder_name}', exist_ok=True)

        # -----------------------------------------------------------------------------
        ## PRE-CALL: calculating the pulse in the given timeframe

        

        
>>>>>>> Massive amount of changes and reorganizing of code. New results
        if modifier_id == 'circle':
            Pulse_x = np.cos(t_vec*w)
            Pulse_y = np.sin(t_vec*w)

            # Pulse on X
            fig, ax = plt.subplots()
            ax.plot(t_vec/T, Pulse_x, color='blue')
            ax.set_xlabel('Time (Periods)')
            ax.set_ylabel('Pulse amplitude')
            fig.suptitle(pulse_suptitle)
            ax.set_title('Pulse shape in x axis')
            fig.savefig(f'Out/{folder_name}/PULSE_X.png', bbox_inches='tight')
            for i in range(len(t_vec_measures)):
                ax.vlines(t_vec_measures[i]/T, -1, 1, color=cmap(norm(t_vec_measures[i])))
            fig.savefig(f'Out/{folder_name}/PULSE_X_meas.png', bbox_inches='tight')
<<<<<<< HEAD
            plt.close(fig)
=======
            
>>>>>>> Massive amount of changes and reorganizing of code. New results

            # Pulse on Y
            fig, ax = plt.subplots()
            ax.plot(t_vec/T, Pulse_y, color='blue')
            ax.set_xlabel('Time (Periods)')
            ax.set_ylabel('Pulse amplitude')
            fig.suptitle(pulse_suptitle)
            ax.set_title('Pulse shape in y axis')
            fig.savefig(f'Out/{folder_name}/PULSE_Y.png', bbox_inches='tight')
            for i in range(len(t_vec_measures)):
                ax.vlines(t_vec_measures[i]/T, -1, 1, color=cmap(norm(t_vec_measures[i])))
            fig.savefig(f'Out/{folder_name}/PULSE_Y_meas.png', bbox_inches='tight')
<<<<<<< HEAD
            plt.close(fig)
=======
            
>>>>>>> Massive amount of changes and reorganizing of code. New results

        elif modifier_id == 'linear':
            Pulse = np.sin(t_vec*w)

            fig, ax = plt.subplots()
            ax.plot(t_vec/T, Pulse, color='blue')
            ax.set_xlabel('Time (Periods)')
            ax.set_ylabel('Pulse amplitude')
            fig.suptitle(pulse_suptitle)
            ax.set_title('Pulse shape')
            fig.savefig(f'Out/{folder_name}/PULSE.png', bbox_inches='tight')
            for i in range(len(t_vec_measures)):
                ax.vlines(t_vec_measures[i]/T, -1, 1, color=cmap(norm(t_vec_measures[i])))
            fig.savefig(f'Out/{folder_name}/PULSE_meas.png', bbox_inches='tight')
<<<<<<< HEAD
            plt.close(fig)
=======

>>>>>>> Massive amount of changes and reorganizing of code. New results
        elif modifier_id == 'linear_packed':
            Pulse = np.sin(t_vec*w)

            fig, ax = plt.subplots()
            ax.plot(t_vec/T, Pulse, color='blue')
            ax.set_xlabel('Time (Periods)')
            ax.set_ylabel('Pulse amplitude')
            fig.suptitle(pulse_suptitle)
            ax.set_title('Pulse shape')
            fig.savefig(f'Out/{folder_name}/PULSE.png', bbox_inches='tight')
            for i in range(len(t_vec_measures)):
                ax.vlines(t_vec_measures[i]/T, -1, 1, color=cmap(norm(t_vec_measures[i])))
            fig.savefig(f'Out/{folder_name}/PULSE_meas.png', bbox_inches='tight')
<<<<<<< HEAD
            plt.close(fig)
=======

>>>>>>> Massive amount of changes and reorganizing of code. New results

        taux = time()
        print('Main calculation starts now.')
        trep = time()
<<<<<<< HEAD
        n_mat, dos_n_mat, dos_mat, t_vec_measures = jcl.kpm_rho_neq(Ham,t_vec,tau,modifier_id,modifier_params,Temp,mu,obs_list,M)
        n_mat_total = np.zeros_like(n_mat)
        dos_n_mat_total = np.zeros_like(dos_n_mat)
        dos_mat_total = np.zeros_like(dos_mat)

        n_mat_total += n_mat 
        dos_n_mat_total += dos_n_mat
        dos_mat_total += dos_mat
        if N_random_vector > 1:
            np.save(f'Out/{folder_name}/Ene_R/1.npy', n_mat[:,:,0])
            np.save(f'Out/{folder_name}/noc_R/1.npy', n_mat[:,:,1])
            np.save(f'Out/{folder_name}/dosn_R/1.npy', dos_n_mat[:,:,1])
            np.save(f'Out/{folder_name}/dos_R/1.npy', dos_mat[:,:,1])

=======
        n_mat, dos_n_mat, t_vec_measures = jcl.kpm_rho_neq(Ham,t_vec,tau,modifier_id,modifier_params,Temp,mu,obs_list,M)
        n_mat_total, dos_n_mat_total = np.zeros_like(n_mat), np.zeros_like(dos_n_mat)
        n_mat_total += n_mat 
        dos_n_mat_total += dos_n_mat
>>>>>>> Massive amount of changes and reorganizing of code. New results
        print('Calc #1 done!')
        print(f'Time elapsed: {timedelta(seconds=time() - trep)}')
        for i in range(2, N_random_vector+1):
            trep = time()
<<<<<<< HEAD
            n_mat, dos_n_mat, dos_mat, t_vec_measures = jcl.kpm_rho_neq(Ham,t_vec,tau,modifier_id,modifier_params,Temp,mu,obs_list,M)

            n_mat_total += n_mat 
            dos_n_mat_total += dos_n_mat
            dos_mat_total += dos_mat

            # Saving specific random vector data
            np.save(f'Out/{folder_name}/Ene_R/{i}.npy', n_mat[:,:,0])
            np.save(f'Out/{folder_name}/noc_R/{i}.npy', n_mat[:,:,1])
            np.save(f'Out/{folder_name}/dosn_R/{i}.npy', dos_n_mat[:,:,1])
            np.save(f'Out/{folder_name}/dos_R/{i}.npy', dos_mat[:,:,1])

=======
            n_mat, dos_n_mat, t_vec_measures = jcl.kpm_rho_neq(Ham,t_vec,tau,modifier_id,modifier_params,Temp,mu,obs_list,M)
            n_mat_total += n_mat 
            dos_n_mat_total += dos_n_mat
>>>>>>> Massive amount of changes and reorganizing of code. New results
            print(f'Calc #{i} done!')
            print(f'Time elapsed: {timedelta(seconds=time() - trep)}')
        n_mat_total /= N_random_vector
        dos_n_mat_total /= N_random_vector
        print('Calculation finished!')
        # Saving results
        EF_list = n_mat_total[:,:,0]
        n_list = n_mat_total[:,:,1]
        ndos_list = dos_n_mat_total[:,:,1]
<<<<<<< HEAD
        dos_list = dos_mat_total[:,:,1]
        np.save(f'Out/{folder_name}/E.npy', EF_list)
        np.save(f'Out/{folder_name}/n_E.npy', n_list)
        np.save(f'Out/{folder_name}/dosn_E.npy', ndos_list)
        np.save(f'Out/{folder_name}/dos_E.npy', dos_list)
=======
        np.save(f'Out/{folder_name}/E.npy', EF_list)
        np.save(f'Out/{folder_name}/n_E.npy', n_list)
        np.save(f'Out/{folder_name}/dosn_E.npy', ndos_list)
>>>>>>> Massive amount of changes and reorganizing of code. New results
        
        print(f'Time elapsed: {timedelta(seconds=time() - taux)}')
    else:
        print('Calculation already performed. Using saved results...')
        EF_list = np.load(f'Out/{folder_name}/E.npy')
        n_list = np.load(f'Out/{folder_name}/n_E.npy')
        ndos_list = np.load(f'Out/{folder_name}/dosn_E.npy')
<<<<<<< HEAD
        dos_list = np.load(f'Out/{folder_name}/dos_E.npy')
=======
>>>>>>> Massive amount of changes and reorganizing of code. New results
    # ------------------------------------------------------------------------------
    #%% Graph results
    # Energy window and lines
    min_e, max_e = -2.5, 2.5
    hw_lines_step = 0.5
    hw_hlines = [i*E for i in [hw_lines_step*i for i in range(1, int(max_e/hw_lines_step+1))]]
    hw_hlines += [i*E for i in [-hw_lines_step*i for i in range(1, int(max_e/hw_lines_step+1))]]

    # Text box
    props = dict(boxstyle='round', facecolor='white', alpha=1.0)

<<<<<<< HEAD
    # Occupation graph
=======
>>>>>>> Massive amount of changes and reorganizing of code. New results
    fig, ax = plt.subplots()
    for i in range(N_measures):
        ax.plot(EF_list[i,:], n_list[i,:], color=cmap(norm(t_vec_measures[i])), label=f't={round(t_vec_measures[i]/T, 3)}T')
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

<<<<<<< HEAD
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

    #%% Extra graphs for the n for different times
    hE_list = [hE*E/2 for hE in range(-hE_reps, hE_reps+1)]
    color_list = ['olivedrab', 'red', 'orange', 'blue', 'darkviolet']
    total_nhE = 2*hE_reps + 1 
    occ_drop_list = np.zeros((total_nhE, N_measures))
    
    for (i, hE) in enumerate(hE_list):
        ind_time = np.where(EF_list[i,:] > hE)[0][0]
        occ_drop_list[i] = np.array([n_list[i,ind_time] for i in range(N_measures)])
        # Specific figure for each calculation
=======

    #%% Extra graphs for the n for different times
    hE_list = [hE*E/2 for hE in range(-hE_reps, hE_reps+1)]
    total_nhE = 2*hE_reps + 1 
    occ_drop_list = np.zeros((total_nhE, N_measures))
    for (i, hE) in enumerate(hE_list):
        ind_time = np.where(EF_list[i,:] > hE)[0][0]
        occ_drop_list[i] = np.array([n_list[i,ind_time] for i in range(N_measures)])

>>>>>>> Massive amount of changes and reorganizing of code. New results
        fig, ax = plt.subplots()
        ax.plot(np.array(t_vec_measures)/T, occ_drop_list[i], c='blue', marker='.', ls='--')
        ax.set_xlabel('Time (Periods)')
        ax.set_ylabel('$n(\\varepsilon)$')
        ax.set_title(fr'Occupation in $E={hE/E:.1f}E$')
        fig.suptitle(fig_title_info)
        fig.savefig(f'Out/{folder_name}/N(T)_{hE:.1f}.png', bbox_inches='tight')
<<<<<<< HEAD
        plt.close(fig)
    # General figure to contain all important data
    fig, ax = plt.subplots()
    reescale = np.max(occ_drop_list[3]) / np.max(occ_drop_list[4]) / 2
    ax.plot(np.array(t_vec_measures)/T, occ_drop_list[4]*reescale, c=color_list[4], 
            marker='.', ls='--', label=f"$E=1E+\\mu$ eV $\\cdot$ {reescale:.3f}")
    ax.plot(np.array(t_vec_measures)/T, occ_drop_list[3], c=color_list[3], 
            marker='.', ls='--', label=f"$E=0.5E+\\mu$ eV")
    ax.plot(np.array(t_vec_measures)/T, occ_drop_list[2], c=color_list[2], 
            marker='.', ls='--', label=f"$E=\\mu$ eV")
    ax.set_xlabel('Time (Periods)')
    ax.set_ylabel('$n(t)$')
    ax.set_title(fr'Occupation')
    fig.suptitle(fig_title_info)
    ax.legend()
    fig.savefig(f'Out/{folder_name}/N(T).png', bbox_inches='tight')

    # ------------------------------------------------------------------------------
    # FREQUENCY ANALYSIS OF THE RESULTS
    ## Performing the Fourier transform
    fourier_occ = np.abs(np.fft.rfft(occ_drop_list))
    ## Obtaining the spacing between frequencies
    dt = t_vec_measures[1] - t_vec_measures[0]
    df = 1/dt/N_measures
    freq = np.arange(0, fourier_occ.shape[1], 1)*df

    # Eliminating the first element of the arrays of freq and occ, as it always
    # explodes
    freq = freq[1:]
    fourier_occ = fourier_occ[:,1:]

    char_freq = np.zeros(len(hE_list))
    for (i, hE) in enumerate(hE_list):
        # Calculating characteristic frequency
        max_freq_ind = (np.where(fourier_occ[i] == max(fourier_occ[i][1:])))[0][0]
        char_freq[i] = freq[max_freq_ind]
        fig, ax = plt.subplots()
=======

    # ------------------------------------------------------------------------------
    # FREQUENCY ANALYSIS OF THE RESULTS
    fourier_occ = np.abs(np.fft.rfft(occ_drop_list))
    dt = t_vec_measures[1] - t_vec_measures[0]
    df = 1/dt/N_measures
    freq = np.arange(0, fourier_occ.shape[1], 1)*df
    for (i, hE) in enumerate(hE_list):
        fig, ax = plt.subplots()
        ax.set_ylim(0, 10)
>>>>>>> Massive amount of changes and reorganizing of code. New results
        ax.plot(freq*T, fourier_occ[i], c='blue', marker='.', ls='--')
        ax.set_xlabel('Frequency (period$^{-1}$)')
        ax.set_ylabel('Amplitude')
        ax.set_title(fr'FFT in $E={hE/E:.1f}E$ + $\mu$')
        fig.suptitle(fig_title_info)
        fig.savefig(f'Out/{folder_name}/FREQ_N(T)_{hE:.1f}.png', bbox_inches='tight')
<<<<<<< HEAD
        plt.close(fig)
    # General figure
    fig, ax = plt.subplots()
    ax.plot(freq*T, fourier_occ[2], c=color_list[2], marker='.', ls='--', 
            label=f'$E = \\mu$ eV, $f_c={char_freq[2]:.6f}$')
    ax.plot(freq*T, fourier_occ[3], c=color_list[3], marker='.', ls='--', 
            label=f'$E = 0.5E + \\mu$ eV, $f_c={char_freq[3]:.6f}$')
    ax.plot(freq*T, fourier_occ[4], c=color_list[4], marker='.', ls='--', 
            label=f'$E = 1E+\\mu$ eV, $f_c={char_freq[4]:.6f}$')
    ax.set_xlabel('Frequency (period$^{-1}$)')
    ax.set_ylabel('Amplitude')
    ax.set_title('FFT')
    fig.suptitle(fig_title_info)
    ax.legend()
    fig.savefig(f'Out/{folder_name}/FREQ_N(T).png', bbox_inches='tight')
    if not show_figs:
        plt.close()
    return EF_list, n_list, ndos_list, dos_list


if __name__ == '__main__':
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
=======
    if not show_figs:
        plt.close()
    return EF_list, n_list, ndos_list
>>>>>>> Massive amount of changes and reorganizing of code. New results
