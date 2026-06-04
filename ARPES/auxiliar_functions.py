import numpy as np
import matplotlib.pyplot as plt
import os

import jclsquant as jcl
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

def plot_pulse(modifier_params, t_vec, t_vec_measures):
    modifier_id, w = modifier_params[0:2]
    T = 2*np.pi/w
    E = jcl.hbar_fs*w
    pulse_suptitle = f'$\\hbar\\omega={round(E, 2)}$ eV, $\\omega$ = {w:.3f} fs$^{{-1}}$, $T={T:.3f}$ fs'
    cmap = plt.cm.plasma 
    norm = Normalize(vmin=t_vec.min(), vmax=t_vec.max())
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
        for i in range(len(t_vec_measures)):
            ax.vlines(t_vec_measures[i]/T, -1, 1, color=cmap(norm(t_vec_measures[i])))

        # Pulse on Y
        fig, ax = plt.subplots()
        ax.plot(t_vec/T, Pulse_y, color='blue')
        ax.set_xlabel('Time (Periods)')
        ax.set_ylabel('Pulse amplitude')
        fig.suptitle(pulse_suptitle)
        ax.set_title('Pulse shape in y axis')
        for i in range(len(t_vec_measures)):
            ax.vlines(t_vec_measures[i]/T, -1, 1, color=cmap(norm(t_vec_measures[i])))

    elif modifier_id == 'linear':
        Pulse = np.sin(t_vec*w)

        fig, ax = plt.subplots()
        ax.plot(t_vec/T, Pulse, color='blue')
        ax.set_xlabel('Time (Periods)')
        ax.set_ylabel('Pulse amplitude')
        fig.suptitle(pulse_suptitle)
        ax.set_title('Pulse shape')
        for i in range(len(t_vec_measures)):
            ax.vlines(t_vec_measures[i]/T, -1, 1, color=cmap(norm(t_vec_measures[i])))

    elif modifier_id == 'linear_packed':
        Tp = modifier_params[2]
        Pulse = np.sin(t_vec*w)/np.cosh((t_vec - 2*Tp)/0.5673/Tp)

        fig, ax = plt.subplots()
        ax.plot(t_vec/T, Pulse, color='blue')
        ax.set_xlabel('Time (Periods)')
        ax.set_ylabel('Pulse amplitude')
        fig.suptitle(pulse_suptitle)
        ax.set_title('Pulse shape')
        for i in range(len(t_vec_measures)):
            ax.vlines(t_vec_measures[i]/T, -1, 1, color=cmap(norm(t_vec_measures[i])))
    plt.show()
    return fig, ax

def load_data(modifier_id, N_pot, E, Temp, mu, gamma, M, N_random_vector, 
              path_type, nk, n_periods, meas_per_T, steps_per_T, type_ham, mass):
    if type_ham == 'basic':
        folder_name = f'{modifier_id}/{path_type}/G={gamma:.3f}_E={float(E)}_Temp={Temp}_mu={mu:.2f}/N={N_pot}_M={M}_nT={n_periods}_measT={meas_per_T}_stT={steps_per_T}_nk={nk}'

    elif type_ham == 'jcl':
        folder_name = f'{modifier_id}{type_ham}/G={gamma:.3f}_E={float(E)}_Temp={Temp}_mu={mu:.2f}/N={N_pot}_M={M}_nT={n_periods}_measT={meas_per_T}_stT={steps_per_T}_nk={nk}'
    
    elif type_ham in ['hbn', 'jclhbn']:
        if mass is None:
            mass = 0.5
        folder_name = f'{modifier_id}{type_ham}/{path_type}/G={gamma:.3f}_E={float(E)}_Temp={Temp}_mu={mu:.2f}_m={mass:.2f}/N={N_pot}_M={M}_nT={n_periods}_measT={meas_per_T}_stT={steps_per_T}_nk={nk}'
   
    save_path = f'ARPES/Out/{folder_name}'
    try:
        cR = len(os.listdir(f'{save_path}/dosn_f'))
        
        if N_random_vector > cR:
            print(f'There are {cR} calculations ready! R not reached, using {cR}/{cR}')
            lR = cR
        else:
            print(f'There are {cR} calculations ready! Enough for R, using {N_random_vector}/{cR}')
            lR = N_random_vector
        EF_list = np.load(f'{save_path}/E.npy')
        t_vec_meas = np.load(f'{save_path}/t_meas.npy')
        Ham = jcl.load_ell(f'{save_path}/Ham.npz')
        n_mat = np.load(f'{save_path}/n/1.npy')
        dosn_mat = np.load(f'{save_path}/dosn/1.npy')
        dosn_f = np.load(f'{save_path}/dosn_f/1.npy')
        for r in range(2, 1+lR):
            n_f += np.load(f'{save_path}/n_f/{r}.npy')
            dosn_f += np.load(f'{save_path}/dosn_f/{r}.npy')
        return EF_list, t_vec_meas, dosn_f/lR, n_mat/lR, dosn_mat/lR
    
    except FileNotFoundError:
        cR = 0
        print('There is no data to load here! Run the calculations or check if a parameter has been misspelled')
        print(f'python3 ARPES/neq_f.py {modifier_id} {N_pot} {E} {Temp} {mu:.2f} {gamma:.3f} {M} {N_random_vector} {path_type} {nk} {n_periods} {meas_per_T} {steps_per_T} {type_ham} {mass:.2f}')
        return 1
    


def string_to_parameters(str_parameters):
    params = str_parameters.split(' ')
    # Type of Light
    modifier_id = params[0]
    # Power to which the number of atoms is 'powered'
    N_pot = int(params[1])
    # Energy in pulse                        
    E = float(params[2])     
    # Temperature                       
    Temp = float(params[3])
    # Chemical potential
    mu = float(params[4])
    # Intensity param     (no units)
    gamma = float(params[5])
    # Amount of moments used to calculate
    M = int(params[6])
    # Amount of random vectors used in calculation
    N_random_vector = int(params[7])
    # Type of path for the system to take
    path_type = params[8]
    # Number of k-points in the first segment
    nk = int(params[9])
    # # periods included in sims
    n_periods = int(params[10])
    # Amount of measures per period
    meas_per_T = int(params[11])
    # steps/T
    steps_per_T = int(params[12])
    # Type of hamiltonian used in the calculations
    type_ham = params[13]
    # Parameter of the given hamiltonian
    mass = float(params[14])
    return modifier_id, N_pot, E, Temp, mu, gamma, M, N_random_vector, path_type, nk, n_periods, meas_per_T, steps_per_T, type_ham, mass