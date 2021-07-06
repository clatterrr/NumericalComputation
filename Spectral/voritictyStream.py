import numpy as np
import random
# https://github.com/nickabattista/Holy_Grail
Nx = 256
a1 = np.zeros((Nx,Nx))
for i in range(Nx):
    for j in range(Nx):
        a1[i,j] = Nx - Nx * 3 // 2 + 1 / 2 + j
        
aD = 25
aR = 38
vort = np.zeros((Nx,Nx))
b1 = np.zeros((Nx,Nx))
b2 = np.zeros((Nx,Nx))
b3 = np.zeros((Nx,Nx))
radius1 = 64

for i in range(Nx):
    for j in range(Nx):
        vort[i,j] = random.random()*2 - 1
        if (a1[i,j]**2 + a1[i,j]**2) < radius1**2:
            vort[i,j] = 0.4
        if (a1[i,j]**2 + (a1[i,j] - aD)**2) < radius1**2:
            vort[i,j] = -0.5
        if ((a1[i,j] + aR)**2 + (a1[i,j] - 2*aD)**2) < radius1**2:
            vort[i,j] = 0.5
            
vort_hat = np.fft.fft2(vort)

kMatx = np.zeros((Nx,Nx),dtype = complex)
kMaty = np.zeros((Nx,Nx),dtype = complex)
for i in range(Nx):
    i0 = i
    if i > Nx//2:
        i0 = i - Nx
    kMatx[:,i] =  1j * i0 
    kMaty[i,:] =  1j * i0 

KLaplace = kMatx**2 + kMaty**2

time = 0
timeFinal = 5
while(time < timeFinal):
    
    psi_hat = np.zeros((Nx,Nx),dtype = complex)
    for i in range(Nx):
        for j in range(Nx):
            if i + j > 1:
                psi_hat[i,j] = -vort_hat[i,j] / KLaplace[i,j]
    
    u = np.fft.ifft2(kMaty * psi_hat).real
    v = np.fft.ifft2(-kMaty * psi_hat).real
    
    vort_x = np.fft.ifft2(kMatx * vort_hat)
    vort_y = np.fft.ifft2(kMaty * vort_hat)
    
    advect = u * vort_x + v * vort_y
    advect_hat = np.fft.fft2(advect)
    
    dt = 0.01
    nu = 0.001
    vort_hat = ( (1 + dt/2*nu*KLaplace )*vort_hat - dt*advect_hat ) / (1 - dt/2*nu*KLaplace)
    
    time = time + 10