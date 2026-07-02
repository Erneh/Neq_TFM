"""
Creates the file containing the parameters to be added to the queue
"""

import numpy as np
import os
import subprocess
from core import check_if_calculated

# Name of the file
file_name = 'absolute_calcs'


# Type of Light
modifier_id_list = ['circln', 'linear']
ham_type_list = ['basic', 'hbn']
ham_param_list = [0.50]
# Power to which the number of atoms is 'powered'
N_pot_list = [20]   
# Temperature                       
Temp_list = [1e-9]
# Chemical potential
mu_list = [0.01]
# (gamma, nT, meas/T, E)
timeframe_intensity_list = [(0.000, 100, 4, 1.0),
                            (0.005, 500, 4, 1.0),
                            (0.010, 200, 8, 1.0),
                            (0.015, 150, 8, 1.0),
                            (0.020, 100, 16, 1.0),
                            (0.025, 100, 16, 1.0),
                            (0.030, 80, 24, 1.0),
                            (0.035, 80, 24, 1.0),
                            (0.040, 50, 32, 1.0),
                            (0.045, 50, 32, 1.0),
                            (0.050, 50, 32, 1.0),
                            (0.000, 100, 4, 0.5),
                            (0.005, 200, 8, 0.5),
                            (0.010, 200, 8, 0.5),
                            (0.015, 150, 8, 0.5),
                            (0.020, 100, 16, 0.5),
                            (0.025, 100, 16, 0.5),
                            (0.030, 80, 24, 0.5),
                            (0.035, 80, 24, 0.5),
                            (0.040, 50, 32, 0.5),
                            (0.045, 50, 32, 0.5),
                            (0.050, 50, 32, 0.5)]
# Amount of random vectors used in calculation
N_random_vector_list = [5]
M_list = [0]
# # periods included in sims
# steps/T
steps_per_T_list = [1000]

# Force the calcs to be redone even if found
# 0: Recalculations are not forced
# 1: Recalculations ARE forced
force_recalc = 1

file = open(f'Calcs_files/{file_name}.txt', 'w')
for N_random_vector in N_random_vector_list:
    for modifier_id in modifier_id_list:
        for N_pot in N_pot_list:
            for Temp in Temp_list:
                for mu in mu_list:
                    for M in M_list:
                        for steps_per_T in steps_per_T_list:
                            for ham_type in ham_type_list:
                                for ham_param in ham_param_list:
                                    for (gamma, n_periods, meas_per_T, E) in timeframe_intensity_list:
                                        flag = check_if_calculated(modifier_id, N_pot, E, Temp, mu, gamma, M, N_random_vector, n_periods, meas_per_T, steps_per_T, ham_type)
                                        if M == 0:
                                            M = int((2**N_pot)**0.5)
                                        if not flag or force_recalc==1:
                                            file.write(f'{modifier_id} {N_pot} {E} {Temp} {mu:.2f} {gamma:.3f} {M} {N_random_vector} {n_periods} {meas_per_T} {steps_per_T} {ham_type} {ham_param:.2f} {force_recalc}\n')
file.close()


file = open(f'Calcs_files/{file_name}.txt', 'r')
n_lines = sum([1 for line in file])
file.close()
print(f'There are {n_lines} calculations queued!')

### After checking the file, this is to launch the actual calculations!


command = f'cat Calcs_files/{file_name}.txt | xargs -L 1 -P 6 ./launcher.sh'
print('screen -S mass_neq')
print(command)