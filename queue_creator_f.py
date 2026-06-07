"""
Creates the file containing the parameters to be added to the queue
"""

import numpy as np
import os
import subprocess

from core import check_if_calculated

# Name of the file
file_name = 'test'


# Type of Light
modifier_id_list = ['linear', 'circle']
ham_type_list = ['basic', 'hbn']
mass_list = [0.50]
# Power to which the number of atoms is 'powered'
N_pot_list = [18, 20]
# Energy in pulse                        
E_list = [1.0]     
# Temperature                       
Temp_list = [1e-9]
# Chemical potential
mu_list = [0.00]
# Intensity param     (no units)
gamma_list = [0.010]
# Amount of random vectors used in calculation
N_random_vector_list = [3]
M_list = [0]
path_type_list = ['full']
nk_list = [100]

# # periods included in sims
n_periods_list = [1]
# Amount of measures per period
meas_per_T_list = [16]
# steps/T
steps_per_T_list = [1000]


file = open(f'ARPES/{file_name}.txt', 'w')
for N_random_vector in N_random_vector_list:
    for modifier_id in modifier_id_list:
        for N_pot in N_pot_list:
            for E in E_list:
                for Temp in Temp_list:
                    for mu in mu_list:
                        for gamma in gamma_list:
                            for M in M_list:
                                for path_type in path_type_list:
                                    for nk in nk_list:
                                        for n_periods in n_periods_list:
                                            for meas_per_T in meas_per_T_list:
                                                for steps_per_T in steps_per_T_list:
                                                    for ham_type in ham_type_list:
                                                        for mass in mass_list:
                                                            if M == 0:
                                                                M = int((2**N_pot)**0.5)
                                                            file.write(f'{modifier_id} {N_pot} {E} {Temp} {mu:.2f} {gamma:.3f} {M} {N_random_vector} {path_type} {nk} {n_periods} {meas_per_T} {steps_per_T} {ham_type} {mass:.2f}\n')
file.close()

file_name = 'test'
file = open(f'ARPES/{file_name}.txt', 'r')
n_lines = sum([1 for line in file])
file.close()
print(f'There are {n_lines} calculations queued!')

### After checking the file, this is to launch the actual calculations!


command = f"grep -v '^#'  ARPES/{file_name}.txt | xargs -L 1 -P 6 ./launcher_f.sh"
print('screen -S mass_arpes')
print(command)