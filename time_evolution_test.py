#%%
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

from scipy.linalg import expm
import jclsquant as jcl
from tqdm import tqdm

# Constants of the system
hw = 1.0
hbar = jcl.hbar_fs 
omega_l = hw / hbar
T = 2*np.pi/omega_l

gamma = 0.025
omega1 = np.sqrt(3)*2.7*gamma*np.pi/jcl.hbar_fs
#rho0 = 0.5*np.ones((2, 2))
#rho0 = np.array([[0, 0], [0, 1]])
delta_w =hw

def Ham(t):
    Ham = np.zeros((2, 2), dtype=np.complex128)
    Ham[0, 0] = delta_w
    Ham[0, 1] = omega1*np.exp(1j*omega_l*t)
    Ham[1, 0] = omega1*np.exp(-1j*omega_l*t)
    Ham[1, 1] = -delta_w
    Ham *= jcl.hbar_fs/2
    return Ham


autV, autE = np.linalg.eigh(Ham(0))
rho0m = np.outer(autE[:,0], autE[:,0])
rho0p =np.outer(autE[:,1], autE[:,1])
psi_0=autE[:,0]
broad = 0.01

dE = broad / 1.5
EF_list = np.arange(autV[0] - 20*dE, autV[1] + 20*dE, dE)


# Actual evolution
N_periods = 100
steps_per_T = 10
t_vec = np.linspace(0, N_periods*T, steps_per_T*N_periods)
dt = t_vec[1]
rho = np.copy(rho0m)
Tracem = np.zeros((len(t_vec)))
Tracep = np.zeros((len(t_vec)))
purity = np.zeros((len(t_vec)),dtype=complex)
DOS_list=np.zeros((len(t_vec), len(EF_list)))
psi=psi_0
for (i, t) in tqdm(enumerate(t_vec)):
 
    autV, autE = np.linalg.eigh(Ham(t))
    rho0m = np.outer(autE[:,0], autE[:,0])
    rho0p =np.outer(autE[:,1], autE[:,1])

        
    Tracem[i] = np.trace(rho@rho0m)
    Tracep[i] = np.trace(rho@rho0p)

    U  = expm(-1j*Ham(t)*dt/jcl.hbar_fs)
    rho = U @ rho 
    psi=U@psi
    purity[i]=np.vdot(autE[:,0],psi)

    for (j, E_step) in enumerate(EF_list):
        DOS_list[i, j] = broad/np.pi*np.trace(np.linalg.inv(np.conjugate(U.T)@(E_step*np.eye(2) - Ham(t))@(E_step*np.eye(2) - Ham(t)) + np.eye(2)*broad**2)@rho)
       
    #rho = np.conj(autE.T) @ np.conj(evol_eigen) @ (autE @ rho @ np.conj(autE.T))@ evol_eigen@ autE

print(purity)
plt.plot(t_vec,np.abs(purity))


plt.plot(EF_list,DOS_list[1000,:])

for i in range(len(t_vec)):

    plt.plot(EF_list,DOS_list[i,:])

print(EF_list[20])
plt.plot(t_vec/T,DOS_list[:,20])
plt.xlim([0,10])



fig, ax = plt.subplots()
cmap = plt.cm.plasma 
norm = Normalize(vmin=t_vec.min()/T, vmax=t_vec.max()/T)
ax.plot(t_vec/T, Tracem[:], c='blue', label='$E_-$')
#ax.plot(t_vec/T, Tracep[:,len(EF_list)//2], c='red', label='$E_+$')
ax.set_xlabel('Time (laser periods)')
ax.set_title('Tr($\\rho \\rho_0$)')
ax.legend()
a=np.array([1,0])
np.outer(a.T,a)

amplitude = np.fft.rfft(Tracem[:])
freqs = np.fft.rfftfreq(len(t_vec), dt)

np.shape(amplitude)
fig, ax = plt.subplots()
ax.plot(jcl.hbar_fs*2*np.pi*freqs, np.abs(amplitude))
ax.set_xlim(0, 3)

print(2*np.pi/(150*T))
amplitude = np.fft.rfft(DOS_list[:,20])
freqs = np.fft.rfftfreq(len(t_vec), dt)

np.shape(amplitude)
fig, ax = plt.subplots()
ax.plot(jcl.hbar_fs*2*np.pi*freqs, np.abs(amplitude))
ax.set_xlim(0, 3)


amplitude = np.fft.rfft(np.abs(purity)**2)
freqs = np.fft.rfftfreq(len(t_vec), dt)

np.shape(amplitude)
fig, ax = plt.subplots()
ax.plot(jcl.hbar_fs*2*np.pi*freqs, np.abs(amplitude))
ax.set_xlim(0, 3)


#for i in range(len(t_vec)):
    #ax.plot(EF_list, DOS_list[i], label=f't={round(t_vec[i]/T, 3)}T')
#cbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax, orientation='vertical', label='Time (periods)')

#%% Time evolution for an specific state

# Constants of the system
hw = 1.0
hbar = jcl.hbar_fs 
omega_l = hw / hbar
T = 2*np.pi/omega_l

gamma = 0.005
omega1 = np.sqrt(3)*2.7*gamma*np.pi/jcl.hbar_fs
#rho0 = 0.5*np.ones((2, 2))
#rho0 = np.array([[0, 0], [0, 1]])
omega0 = omega_l

def Ham(t):
    Ham = np.zeros((2, 2), dtype=np.complex128)
    Ham[0, 0] = omega0
    Ham[0, 1] = omega1*np.exp(1j*omega_l*t)
    Ham[1, 0] = omega1*np.exp(-1j*omega_l*t)
    Ham[1, 1] = -omega0
    Ham *= jcl.hbar_fs/2
    return Ham


psi0 = np.array([1, 0], dtype=complex)


# Actual evolution
N_periods = 10
steps_per_T = 20000
t_vec = np.linspace(0, N_periods*T, steps_per_T*N_periods)
dt = t_vec[1]
psi = np.zeros((len(t_vec), 2))
psi[0] = np.copy(psi0)
for i in tqdm(range(1, len(t_vec))):
    U  = expm(-1j*Ham(t_vec[i])*dt/jcl.hbar_fs)
    psi[i]=U@psi[i-1]

fig, ax = plt.subplots()
ax.plot(t_vec/T, np.abs(psi[:,1]))
ax.plot(t_vec/T, np.abs(psi[:,0]))
ax.plot(t_vec/T, np.sum(psi**2, axis=1)**0.5)
