"""
Creates the file containing the parameters to be added to the queue
"""

import numpy as np
import os
import subprocess

def check_if_calculated(modifier_id, N_pot, E, Temp, mu, gamma, M, N_random_vector,
                        n_periods, meas_per_T, steps_per_T):
    N = 2**N_pot
    # Names of file and info on graphs
    folder_name = f'{modifier_id}/G={gamma:.3f}_E={float(E)}_Temp={Temp}_mu={mu:.2f}/N={N_pot}_M={M}_R={N_random_vector}_nT={n_periods}_measT={meas_per_T}_stT={steps_per_T}'
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
            np.load(f'Out/{folder_name}/n_E.npy')
            np.load(f'Out/{folder_name}/dosn_E.npy')
            np.load(f'Out/{folder_name}/dos_E.npy')
            if np.sum(np.abs(EF_list)) == 0:
                flag = False
        except FileNotFoundError:
            flag = False
    return flag

# Name of the file
file_name = 'aa_calcs'


# Type of Light
modifier_id_list = ['circle']

# Power to which the number of atoms is 'powered'
N_pot_list = [19]
# Energy in pulse                        
E_list = [0.5]     
# Temperature                       
Temp_list = [1e-9]
# Chemical potential
mu_list = [0.01]
# Intensity param     (no units)
gamma_list = np.linspace(0.000, 0.050, 11)
# Amount of random vectors used in calculation
N_random_vector_list = [5]
M_list = [0]
# # periods included in sims
n_periods_list = [200]
# Amount of measures per period
meas_per_T_list = [32]
# steps/T
steps_per_T_list = [1000]

# Force the calcs to be redone even if found
# 0: Recalculations are not forced
# 1: Recalculations ARE forced
force_recalc = 1

file = open(f'Calcs_files/{file_name}.txt', 'w')

for modifier_id in modifier_id_list:
    for N_pot in N_pot_list:
        for E in E_list:
            for Temp in Temp_list:
                for mu in mu_list:
                    for gamma in gamma_list:
                        for M in M_list:
                            for N_random_vector in N_random_vector_list:
                                for n_periods in n_periods_list:
                                    for meas_per_T in meas_per_T_list:
                                        for steps_per_T in steps_per_T_list:
                                            flag = check_if_calculated(modifier_id, N_pot, E, Temp, mu, gamma, M, N_random_vector,
                            n_periods, meas_per_T, steps_per_T)
                                            if not flag or force_recalc==1:
                                                if M == 0:
                                                    M = int((2**N_pot)**0.5)
                                                file.write(f'{modifier_id} {N_pot} {E} {Temp} {mu:.2f} {gamma:.3f} {M} {N_random_vector} {n_periods} {meas_per_T} {steps_per_T} {force_recalc}\n')
file.close()


file = open(f'Calcs_files/{file_name}.txt', 'r')
n_lines = sum([1 for line in file])
file.close()
print(f'There are {n_lines} calculations queued!')

### After checking the file, this is to launch the actual calculations!


command = f'cat {file_name}.txt | xargs -L 1 -P 6 ./launcher.sh'
print('screen -S mass_neq')
print(command)