"""
Creates the file containing the parameters to be added to the queue
"""

import numpy as np
import os
import subprocess

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

# Name of the file
file_name = 'checking'


# Type of Light
modifier_id_list = ['linear', 'circle']

# Power to which the number of atoms is 'powered'
N_pot_list = [17]
# Energy in pulse                        
E_list = [1.0]     
# Temperature                       
Temp_list = [1e-9]
# Chemical potential
mu_list = [0.01]
# Intensity param     (no units)
gamma_list = [0.010]
# Amount of random vectors used in calculation
N_random_vector_list = [1, 3, 5, 10]

# # periods included in sims
n_periods_list = [20, 50, 100]
# Amount of measures per period
meas_per_T_list = [4, 8, 16, 32]
# steps/T
steps_per_T_list = [200, 500, 1000]

# Force the calcs to be redone even if found
# 0: Recalculations are not forced
# 1: Recalculations ARE forced
force_recalc = 0

file = open(f'{file_name}.txt', 'w')

for modifier_id in modifier_id_list:
    for N_pot in N_pot_list:
        for E in E_list:
            for Temp in Temp_list:
                for mu in mu_list:
                    for gamma in gamma_list:
                        for N_random_vector in N_random_vector_list:
                            for n_periods in n_periods_list:
                                for meas_per_T in meas_per_T_list:
                                    for steps_per_T in steps_per_T_list:
                                        flag = check_if_calculated(modifier_id, N_pot, E, Temp, mu, gamma, N_random_vector,
                        n_periods, meas_per_T, steps_per_T)
                                        if not flag or force_recalc==1:
                                            file.write(f'{modifier_id} {N_pot} {E} {Temp} {mu} {gamma} {N_random_vector} {n_periods} {meas_per_T} {steps_per_T} {force_recalc}\n')
file.close()


file = open(f'{file_name}.txt', 'r')
n_lines = sum([1 for line in file])
file.close()
print(f'There are {n_lines} calculations queued!')

### After checking the file, this is to launch the actual calculations!


command = f'cat {file_name}.txt | xargs -L 1 -P 6 ./launcher.sh'
print('screen -S mass_neq')
print(command)