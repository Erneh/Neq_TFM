import numpy as np
import scipy.sparse as sp
import os

from scipy.spatial import cKDTree
from scipy.fft import fft, fftfreq, rfft, rfftfreq   

# Random vector generation
def random_vector(N, n_rand):
    return np.exp(np.random.random((N, n_rand))*2j*np.pi)

# Kerner generation to avoid excessive noise
def jackson_kernel(m, M):
    c1 = (np.pi)/(M+1)
    return ((M+1-m)*np.cos(m*c1) + np.sin(c1*m)/np.tan(c1))/(M+1)


# Getting the density of states!
def DOS_sparse(Ham, E_tol, E_list, M, n_random = 10):
    # Prepare hamiltonian
    N = Ham.shape[0]
    dE = (E_tol[1] - E_tol[0])/2
    mE = (E_tol[1] + E_tol[0])/2
    Ham_n = (Ham - mE*sp.eye(N))/dE
    # Normalizing energies
    En_list = (E_list - mE)/dE
    # Initializing Cheb_H coefficients
    Cheb_H = np.zeros(M+1)
    Cheb_E = np.zeros((M+1, len(En_list)))
    rvec = random_vector(N, n_random)

    # Recursive vectors for H
    v1 = np.copy(rvec)
    v2 = Ham_n.dot(v1)

    # Recursive vectors for E
    E_1 = np.ones(len(En_list))
    E_2 = En_list

    for i in range(M+1):
        # Updating newest vector
        v3 = 2*Ham_n.dot(v2) - v1
        E_3 = 2*E_2*En_list - E_1
        
        # Calculating the terms
        Cheb_H[i] = np.sum(np.diag(np.conj(rvec.T) @ v1).real)/n_random
        Cheb_E[i] = E_1
        # Renewing saved terms
        v1, v2 = v2, v3
        E_1, E_2 = E_2, E_3

        
    # Energy calculation
    delta_E = 2*Cheb_E/(np.pi*np.sqrt(1-En_list[None,:]**2))
    delta_E[0] = 1/(np.pi*np.sqrt(1-En_list[None,:]**2))

    g_Jackson = jackson_kernel(np.arange(0, M+1, 1), M)
    return np.sum(delta_E*g_Jackson[:,None]*Cheb_H[:,None], axis=0)/dE/N


def check_if_calculated(modifier_id, N_pot, E, Temp, mu, gamma, M, N_random_vector,
                        n_periods, meas_per_T, steps_per_T, type_ham, ham_params=None):
    N = 2**N_pot
    # Names of file and info on graphs
    if type_ham == 'hbn':
        if ham_params is None:
            mass = 0.5
        else:
            mass = ham_params
        folder_name = f'{modifier_id}{type_ham}/G={gamma:.3f}_E={float(E)}_Temp={Temp}_mu={mu:.2f}_m={mass:.2f}/N={N_pot}_M={M}_R={N_random_vector}_nT={n_periods}_measT={meas_per_T}_stT={steps_per_T}'
    else:
        folder_name = f'{modifier_id}{type_ham}/G={gamma:.3f}_E={float(E)}_Temp={Temp}_mu={mu:.2f}/N={N_pot}_M={M}_R={N_random_vector}_nT={n_periods}_measT={meas_per_T}_stT={steps_per_T}'
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


def load_data(modifier_id, N_pot, E, Temp, mu, gamma, M, N_random_vector,
                        n_periods, meas_per_T, steps_per_T, type_ham, ham_params, R=None):
    '''
    Script to load the data of an already performed calculation
    '''
    if type_ham == 'hbn':
        if ham_params is None:
            mass = 0.5
        else:
            mass = ham_params
        folder_name = f'{modifier_id}{type_ham}/G={gamma:.3f}_E={float(E)}_Temp={Temp}_mu={mu:.2f}_m={mass:.2f}/N={N_pot}_M={M}_R={N_random_vector}_nT={n_periods}_measT={meas_per_T}_stT={steps_per_T}'
    else:
        folder_name = f'{modifier_id}{type_ham}/G={gamma:.3f}_E={float(E)}_Temp={Temp}_mu={mu:.2f}/N={N_pot}_M={M}_R={N_random_vector}_nT={n_periods}_measT={meas_per_T}_stT={steps_per_T}'

    # Loading the info in the .npy files available at the moment
    if N_random_vector == 1:
        print(f'1/1 calculations finished! Loading results...')
        EF_list = np.load(f'Out/{folder_name}/E.npy')
        n_E_list = np.load(f'Out/{folder_name}/n_E.npy')
        dos_list = np.load(f'Out/{folder_name}/dos_E.npy')
        dosn_list = np.load(f'Out/{folder_name}/dosn_E.npy')
    else:
        # Deciding the amount of random vectors being used in the sim
        R_total = len(os.listdir(f'Out/{folder_name}/Ene_R'))
        if R is None:
            R = R_total

        print(f'{R_total}/{N_random_vector} calculations finished! Loading {R} results...')
        EF_list = np.load(f'Out/{folder_name}/Ene_R/1.npy')
        n_E_list = np.zeros_like(EF_list)
        dos_list = np.zeros_like(EF_list)
        dosn_list = np.zeros_like(EF_list)
        for i in range(1, R+1):
            n_E_list += np.load(f'Out/{folder_name}/noc_R/{i}.npy')
            dos_list += np.load(f'Out/{folder_name}/dos_R/{i}.npy')
            dosn_list += np.load(f'Out/{folder_name}/dosn_R/{i}.npy')
        n_E_list /= R
        dos_list /= R
        dosn_list /= R
    return EF_list, n_E_list, dos_list, dosn_list


# Completely experimental!!! Do not use, it WONT WORK
def create_ham_rel(S, rLat, thresh_1NN, t=-2.7, M=0.0, per=None):
    N = S.shape[0]
    # Creating inverse to 
    rLat_inv = np.linalg.inv(rLat)
    # Adding the hopping terms
    S_rel = S @ rLat_inv

    if per is None:
        ckd = cKDTree(S_rel)
        inds = ckd.query_pairs(thresh_1NN, output_type='ndarray')
    else:
        ckd = cKDTree(S_rel, boxsize=(per[0], per[1]))
        inds = ckd.query_pairs(thresh_1NN, output_type='ndarray')

    N_1NN = len(inds)
    Ham = sp.csr_matrix((np.ones(N_1NN)*t, (inds[:,0], inds[:,1])), shape=(N, N), 
                        dtype=np.complex128)
    Mass = np.kron(np.ones(N//2), np.array([M, -M]))
    Ham += Mass 
    return Ham
    

def extract_occ_time(t_vec_measures, EF_list, n_E_list, hE_list):
    total_nhE = len(hE_list)
    N_measures = len(t_vec_measures)
    occ_drop_list = np.zeros((total_nhE, N_measures))
    
    for (i, hE) in enumerate(hE_list):
        ind_time = np.where(EF_list[i,:] > hE)[0][0]
        #occ_drop_list[i] = np.array([n_E_list[i,ind_time] for i in range(N_measures)])
        occ_drop_list[i] = n_E_list[:,ind_time]
    return occ_drop_list, N_measures


def frequency_analysis(EF_list, n_E_list, hE_list, t_vec_measures, T, range_search=1):
    """
    Returns ANGULAR frequency and chararcteristic freq of a given occupation
    """
    occ_drop_list, N_measures = extract_occ_time(t_vec_measures, EF_list, n_E_list, hE_list)

    fourier_occ = np.abs(fft(occ_drop_list))[:,:N_measures//2]
    dt = t_vec_measures[1] - t_vec_measures[0]
    #df = 1/dt/N_measures
    #freq = np.arange(0, fourier_occ.shape[1], 1)*df*(2*np.pi)
    freq = fftfreq(N_measures, dt)[:N_measures//2]*(2*np.pi)
    df = freq[1] - freq[0]
    char_freq = np.zeros(occ_drop_list.shape[0])
    max_freq_ind = np.zeros(occ_drop_list.shape[0], dtype=np.int64)
    # Eliminating the first element of the arrays of freq and occ, as it always explodes
    freq = freq[1:]
    fourier_occ = fourier_occ[:,1:]
    # Checking only the frequencies between 0 and 1 (in laser period units)
    inf_freqs = freq/(2*np.pi)*T < range_search - df/(2*np.pi)
    for (i, hE) in enumerate(hE_list):
        #max_freq_ind[i] = (np.where(fourier_occ[i, inf_freqs] == max(fourier_occ[i][inf_freqs])))[0][0]
        max_freq_ind[i] = fourier_occ[i, inf_freqs].argmax()
        # Rule out noisy situation where no real mode is detected
        if fourier_occ[i, inf_freqs][max_freq_ind[i]] < 1e-10:
            max_freq_ind[i] = 0
        char_freq[i] = freq[max_freq_ind[i]]
    return occ_drop_list, fourier_occ, freq, char_freq, max_freq_ind
