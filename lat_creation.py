import numpy as np
import matplotlib.pyplot as plt

import scipy.sparse as sp



# General Lattice points, ordered in the same way as the original hamiltonian
def get_positions(N1, N2, rLat, Rat):
    # Setting up the displacement by the lattice
    N1_list = np.arange(0, N1, 1, dtype=int)
    N2_list = np.arange(0, N2, 1, dtype=int)

    N1_exp = np.kron(N1_list, np.ones(N2, dtype=int))
    N2_exp = np.kron(np.ones(N1, dtype=int), N2_list)
    N_grid = np.vstack([N1_exp, N2_exp])

    # Expanding by # of atoms in original unit cell
    N_grid = np.kron(N_grid, np.ones(len(Rat))).T
    R_grid = N_grid @ rLat

    # Preparing the positions of the atoms
    Rat_exp = np.kron(np.ones(N1*N2), Rat.T).T

    # Final atom positions
    return Rat_exp + R_grid


def get_positions_graphene(N1, N2):
    a_l = 0.24595
    a_cc =a_l/3**0.5
    a1=np.array([np.sqrt(3)*a_cc/2,a_cc/2])
    a2=np.array([0,a_cc])

    Rat = np.zeros((4, 2))
    Rat[1] = a2
    Rat[2] = a2 + a1
    Rat[3] = 2*a2 + a1

    A1 = np.array([a_l, 0])
    A2 = np.array([0.0, 3**0.5*a_l])
    rLat = np.array([A1, A2])

    return get_positions(N1, N2, rLat, Rat)