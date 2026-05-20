import numpy as np
import scipy.sparse as sp

from scipy.spatial import cKDTree

import jclsquant as jcl

# Base general hamiltonian with 1NN neighbour definition
def create_ham(S, N1, N2, rLat, thresh_1NN, t=-2.7, M=0.0, per=True, out_format = 'CSR'):
    N = S.shape[0]
    # Only if rLat is square!!! In any other case, things change
    Lx, Ly = N1*rLat[0][0], N2*rLat[1][1]
    L = np.array([Lx, Ly])
    if per:
        ckd = cKDTree(S, boxsize=(Lx, Ly))
        inds = ckd.query_pairs(thresh_1NN, output_type='ndarray')
    else:
        ckd = cKDTree(S)
        inds = ckd.query_pairs(thresh_1NN, output_type='ndarray')
        
    inds_conj = inds[:,::-1]
    conj_list = np.concatenate((np.zeros(len(inds), dtype=int), np.ones(len(inds_conj), dtype=int)))
    inds = np.concatenate((inds, inds_conj))
    
    N_1NN = len(inds)
    Ham = sp.csr_matrix((np.ones(N_1NN)*t, (inds[:,0], inds[:,1])), shape=(N, N), 
                        dtype=np.complex128)
    Mass = np.kron(np.ones(N//2), np.array([M, -M]))
    Ham += sp.diags(Mass) 

    # Position matrices
    Pos_dif = S[inds[:,1]] - S[inds[:,0]] 
    # Make periodic part
    Pos_dif = Pos_dif - np.floor(Pos_dif/L + 0.5)*L
    # Relative dist matrices
    Vx = sp.csr_matrix((Pos_dif[:,0], (inds[:,0], inds[:,1])), shape=Ham.shape) + 1e-16*Ham
    Vy = sp.csr_matrix((Pos_dif[:,1], (inds[:,0], inds[:,1])), shape=Ham.shape) + 1e-16*Ham

    # Conjugate matrix
    conjugate_mat = sp.csr_matrix((conj_list, (inds[:,0], inds[:,1])), shape=(N, N))
    if np.abs(M) != 0:
        conjugate_mat += sp.diags(np.ones(N)*2, shape=(N, N))
    # Different outs depending on format specified
    if out_format == 'CSR':
        return Ham, Vx, Vy, conjugate_mat
    
    elif out_format == 'ELL':
        H_ell = jcl.ell_matrix(Ham, Vx, Vy, conjugate_mat)
        H_ell.Omega = Lx*Ly
        return H_ell


def create_hex_ham(S, N1, N2, t=-2.7, M=0.0, a_l = 0.24595, per=True, out_format='CSR'):
    """
    Creates a hamiltonian pertaining to graphene with less inputs than normal and in
    the ELL format to be compatible with the rest of the code

    By default, it does graphene, but it can make anythinw with the same lattice
    vectors
    """
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

    if per:
        per_info = (N1, N2)
    else:
        per_info = None

    return create_ham(S, N1, N2, rLat, a_cc*1.1, t, M, per_info, out_format)