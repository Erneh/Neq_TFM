import numpy as np
import matplotlib.pyplot as plt

import jclsquant as jcl

from ham_creation import create_hex_ham
from lat_creation import get_positions_graphene


N_pot = 17
N = 2**N_pot
N1 = N2 = int(np.sqrt(N))//2
M = 0.5
t = -2.7
a = 0.25
S = get_positions_graphene(N1, N2, a_l=a)
Ham_ell = create_hex_ham(S, N1, N2, t=t, M=M, a_l=a, out_format='ELL')
Ham, Vx, Vy, conjugate_mat = create_hex_ham(S, N1, N2, t=t, M=M, a_l=a, out_format='CSR')

# DOS calculations
res = jcl.kpm_dos(Ham_ell)

fig, ax = plt.subplots(dpi=200)
ax.plot(res[:,0], res[:,1], c='blue')
ax.set_xlabel('Energy (eV)')
ax.set_ylabel('DOS')
ax.set_xlim(-2, 2)
print(res)
