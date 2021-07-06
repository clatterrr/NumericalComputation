import numpy as np
'''
https://github.com/nickabattista/Holy_Grail

dv/dt = D * laplacian(v) + v*(v-a)*(v-1) - w - I(t)
dw/dt = eps*(v - gamma*w)

v(x,t): membrane potential
w(x,t): blocking mechanism
D: diffusion rate of potential
a: threshold potential
gamma: resetting rate
eps: strength of blocking
I(t):initial condition for applied activation
'''

D = 1 # Diffusion coefficient
a = 0.3 # Threshold potential
gamma = 1 # Resetting rate
eps = 0.001 # Blocking strength
I_mag = 0.05 # Activation strength

N = 800
L = 2000
dx = L / N

timeFinal = 10
Npulse = 10
pulse = timeFinal / Npulse
dt = 0.001
i1 = 0.475 # fraction of total length where current starts
i2 = 0.525 # fraction of total length where current ends

v = np.zeros((N+1))
time = 0
while(time < timeFinal):
    
    DDv = np.zeros((N+1))
    DDv[0] = (v[1] - 2*v[0] + v[N])/(dx**2)
    DDv[N] = (v[0] - 2*v[N] + v[N-1])/(dx**2)
    for i in range(1,N):
        DDv[i] = (v[i-1] + v[i+1] - 2*v[i])/(dx**2)
        
    # activation time
    if t > pulse_time:
        for i in range(N):
            app[j] = I_mag
    
    time = time + 10
