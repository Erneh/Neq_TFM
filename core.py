import numpy as np
import scipy.sparse as sp
import os

from scipy.spatial import cKDTree
from scipy.fft import fft, fftfreq, rfft, rfftfreq   


hbar = 0.6582
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


def str_parameters(str_params):
    """
    Assign parameters by using the same string that is used to calculate them in 
    the first place
    """
    # Getting the parameters from the strings
    params = [i for i in str_params.split(' ') if len(i) > 0]
    # Assigning to specific variables
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

    # # periods included in sims
    n_periods = int(params[8])
    # Amount of measures per period
    meas_per_T = int(params[9])
    # steps/T
    steps_per_T = int(params[10])
    # Type of hamiltonian used in the calculations
    type_ham = params[11]
    # Parameter of the given hamiltonian
    ham_param = float(params[12])

    return modifier_id, N_pot, E, Temp, mu, gamma, M, N_random_vector, n_periods, meas_per_T, steps_per_T, type_ham, ham_param

def str_to_dict(str_params):
    # Getting the parameters from the strings
    params = [i for i in str_params.split(' ') if len(i) > 0]
    # Assigning to specific variables
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

    # # periods included in sims
    n_periods = int(params[8])
    # Amount of measures per period
    meas_per_T = int(params[9])
    # steps/T
    steps_per_T = int(params[10])
    # Type of hamiltonian used in the calculations
    type_ham = params[11]
    # Parameter of the given hamiltonian
    ham_param = float(params[12])
    res = {'modifier_id': modifier_id,
        'type_ham': type_ham,
        'ham_params': ham_param,  
        'hw' : E,                       
        'Temp' :  Temp,
        'mu' : mu,                 
        'N_pot' : N_pot,
        'n_periods' : n_periods,
        'steps_per_T' : steps_per_T,
        'meas_per_T' : meas_per_T,
        'N_random_vector' : N_random_vector,
        'M' : M,
        'gamma' : gamma
        }
    return res

def extract_from_dict(res_d):
    # Assigning to specific variables
    modifier_id = res_d['modifier_id']
    # Power to which the number of atoms is 'powered'
    N_pot = res_d['N_pot']
    # Energy in pulse                        
    E = res_d['hw']
    # Temperature                       
    Temp = res_d['Temp']
    # Chemical potential
    mu = res_d['mu']
    # Intensity param     (no units)
    gamma = res_d['gamma']
    # Amount of moments used to calculate
    M = res_d['M']
    # Amount of random vectors used in calculation
    N_random_vector = res_d['N_random_vector']
    # # periods included in sims
    n_periods = res_d['n_periods']
    # Amount of measures per period
    meas_per_T = res_d['meas_per_T']
    # steps/T
    steps_per_T = res_d['steps_per_T']
    # Type of hamiltonian used in the calculations
    type_ham = res_d['type_ham']
    # Parameter of the given hamiltonian
    ham_param = res_d['ham_params']
    return modifier_id, N_pot, E, Temp, mu, gamma, M, N_random_vector, n_periods, meas_per_T, steps_per_T, type_ham, ham_param


def load_data_dict(res_d, R=None, out_file_loc=''):
    modifier_id, N_pot, E, Temp, mu, gamma, M, N_random_vector, n_periods, meas_per_T, steps_per_T, type_ham, ham_param = extract_from_dict(res_d)
    EF_list, n_E_list, dos_list, dosn_list = load_data(modifier_id, N_pot, E, Temp, mu, gamma, M, N_random_vector, n_periods, meas_per_T, steps_per_T, type_ham, ham_param, R, out_file_loc)
    res_d['EF_list'] = EF_list
    res_d['n_E_list'] = n_E_list
    res_d['dos_list'] = dos_list
    res_d['dosn_list'] = dosn_list
    return res_d


def load_data(modifier_id, N_pot, E, Temp, mu, gamma, M, N_random_vector,
                        n_periods, meas_per_T, steps_per_T, type_ham, ham_params, R=None, out_file_loc=''):
    '''
    Script to load the data of an already performed calculation
    '''
    if type_ham == 'hbn':
        if ham_params is None:
            mass = 0.5
        else:
            mass = ham_params
        folder_name = f'{modifier_id}{type_ham}/G={gamma:.3f}_E={float(E)}_Temp={Temp}_mu={mu:.2f}_m={mass:.2f}/N={N_pot}_M={M}_R={N_random_vector}_nT={n_periods}_measT={meas_per_T}_stT={steps_per_T}'
    elif type_ham == 'basic':
        folder_name = f'{modifier_id}/G={gamma:.3f}_E={float(E)}_Temp={Temp}_mu={mu:.2f}/N={N_pot}_M={M}_R={N_random_vector}_nT={n_periods}_measT={meas_per_T}_stT={steps_per_T}'
    else:
        folder_name = f'{modifier_id}{type_ham}/G={gamma:.3f}_E={float(E)}_Temp={Temp}_mu={mu:.2f}/N={N_pot}_M={M}_R={N_random_vector}_nT={n_periods}_measT={meas_per_T}_stT={steps_per_T}'

    # Loading the info in the .npy files available at the moment
    if N_random_vector == 1:
        print(f'1/1 calculations finished! Loading results...')
        EF_list = np.load(f'{out_file_loc}/Out/{folder_name}/E.npy')
        n_E_list = np.load(f'{out_file_loc}/Out/{folder_name}/n_E.npy')
        dos_list = np.load(f'{out_file_loc}/Out/{folder_name}/dos_E.npy')
        dosn_list = np.load(f'{out_file_loc}/Out/{folder_name}/dosn_E.npy')
    else:
        # Deciding the amount of random vectors being used in the sim
        R_total = len(os.listdir(f'{out_file_loc}/Out/{folder_name}/Ene_R'))
        if R is None:
            R = R_total

        print(f'{R_total}/{N_random_vector} calculations finished! Loading {R} results...')
        EF_list = np.load(f'{out_file_loc}/Out/{folder_name}/Ene_R/1.npy')
        n_E_list = np.zeros_like(EF_list)
        dos_list = np.zeros_like(EF_list)
        dosn_list = np.zeros_like(EF_list)
        for i in range(1, R+1):
            n_E_list += np.load(f'{out_file_loc}/Out/{folder_name}/noc_R/{i}.npy')
            dos_list += np.load(f'{out_file_loc}/Out/{folder_name}/dos_R/{i}.npy')
            dosn_list += np.load(f'{out_file_loc}/Out/{folder_name}/dosn_R/{i}.npy')
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


def frequency_analysis(EF_list, n_E_list, hE_list, t_vec_measures, T):
    """
    Returns ANGULAR frequency and chararcteristic freq of a given occupation
    """
    occ_drop_list, N_measures = extract_occ_time(t_vec_measures, EF_list, n_E_list, hE_list)

    fourier_occ = np.abs(np.fft.fft(occ_drop_list))[:,:N_measures//2]
    dt = t_vec_measures[1] - t_vec_measures[0]
    #df = 1/dt/N_measures
    #freq = np.arange(0, fourier_occ.shape[1], 1)*df*(2*np.pi)
    freq = np.fft.fftfreq(N_measures, dt)[:N_measures//2]*(2*np.pi)
    df = freq[1] - freq[0]
    char_freq = np.zeros(occ_drop_list.shape[0])
    max_freq_ind = np.zeros(occ_drop_list.shape[0], dtype=np.int64)

    # Eliminating the first element of the arrays of freq and occ, as it always explodes
    freq = freq[1:]
    fourier_occ = fourier_occ[:,1:]

    # Preparing a gaussian filder
    w = (2*np.pi/T)
    max_las_freq = round(freq[-1] / w)
    std = 2*df
    gaussians = (std*np.sqrt(2*np.pi))*np.exp(-0.5*(freq[None,:]-w*np.arange(1, max_las_freq+0.1, 1)[:,None])**2/std**2)
    f_gaussian = 1 - np.sum(gaussians, axis=0) / np.sum(gaussians, axis=0).max()
    f_gaussian[f_gaussian < 0.6] = 0.0
    # Checking only the frequencies between 0 and 1 (in laser period units)
    max_occ = occ_drop_list.max()
    for (i, hE) in enumerate(hE_list):
        #max_freq_ind[i] = (np.where(fourier_occ[i, inf_freqs] == max(fourier_occ[i][inf_freqs])))[0][0]
        max_freq_ind[i] = (fourier_occ[i]*f_gaussian).argmax()
        char_freq[i] = freq[max_freq_ind[i]]
        # Rule out noisy situation where no real mode is detected
        if max(np.abs(fourier_occ[i][max_freq_ind[i]] - fourier_occ[i])) < 1:
            max_freq_ind[i] = 0
            char_freq[i] = .0
        
    return occ_drop_list, fourier_occ, freq, char_freq, max_freq_ind


def frequency_analysis_dict(res, hE_list = None, use_dosn=False, hE_reps=2):
    hw = res['hw']
    w = hw/hbar
    T = 2*np.pi/w
    t_vec_measures = np.linspace(0, res['n_periods']*T, res['n_periods']*res['meas_per_T'])
    if hE_list is None:
        hE_list = [hE*hw/2 for hE in range(-hE_reps, hE_reps+1)]
    if use_dosn:
        occ_drop_list, fourier_occ, freq, char_freq, max_freq_ind = frequency_analysis(res['EF_list'], res['dosn_list'], hE_list, t_vec_measures, T)
    else:
        occ_drop_list, fourier_occ, freq, char_freq, max_freq_ind = frequency_analysis(res['EF_list'], res['n_E_list'], hE_list, t_vec_measures, T)
    res['occ_drop_list'] = occ_drop_list
    res['fourier_occ'] = fourier_occ
    res['freq'] = freq
    res['char_freq_s'] = char_freq
    res['max_freq_ind'] = max_freq_ind
