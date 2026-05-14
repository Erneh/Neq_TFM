import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import jclsquant as jcl
from scipy.linalg import expm
# Constants of the sim
hw = 0.5
w = hw/jcl.hbar_fs
T = 2*np.pi/w      
e = 1
nu = 1
gamma = 0.005
Phi0 = 2*np.pi*jcl.hbar_fs
A0 = gamma*Phi0/(2*jcl.a_cc*np.sqrt(3))
vf = 3*jcl.a_cc*2.7/(2*jcl.hbar_fs)
kx, ky = 0, 0
km = np.sqrt(kx**2 + ky**2)
alpha = np.arctan2(ky, kx)

# Try to evolve using exponential matrices
def H0(t):
    C = A0/(1j*nu*w*jcl.hbar_fs*t)
    return jcl.hbar_fs*vf*nu*np.array([[0, km*np.exp(-1j*nu*alpha) + C*(np.exp(1j*nu*w*t)-1)],
                                [km*np.exp(1j*nu*alpha) + C*(np.exp(-1j*nu*w*t)-1), 0]])
State0 =np.array([1, 0])
N_periods = 20
steps_per_period = 1000
t0, tf = 0, 100*T
t_vec = np.linspace(t0, tf, N_periods*steps_per_period)
dt = t_vec[1] - t_vec[0]
State = np.zeros((len(t_vec), 2))
State[0] = State0
U = expm(-1j*dt*H0(dt)/jcl.hbar_fs)
for i in range(1, len(t_vec)):
    State[i] = U @ State[i-1]



# Equations from the perturbation theory
def fun(t, z):
    cp, cm = z
    dcp = np.exp(2j*vf*km*t)*vf*e*A0/(2j*jcl.hbar_fs)*(-np.exp(-1j*nu*(2*alpha+w*t)+nu*np.exp(1j*nu*w*t)))*cm
    dcm = np.exp(-2j*vf*km*t)*vf*e*A0/(2j*jcl.hbar_fs)*(-np.exp(1j*nu*(2*alpha+w*t)+nu*np.exp(-1j*nu*w*t)))*cp
    return [dcp, dcm]

y0 = [1, 0]


sol = solve_ivp(fun, (t0, tf), y0)
sol.t
