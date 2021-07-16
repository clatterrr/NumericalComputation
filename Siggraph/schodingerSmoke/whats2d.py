import numpy as np
import math
volsize = [4,2]
volres = [64,32]
hbar = 0.1 # 普朗克常量
dt = 1 / 48
tmax = 0.05

dx = volsize[0] / volres[0]
dy = volsize[1] / volres[1]

norm = np.sqrt(1*1 + 0.01*0.01)
psi1 = np.ones((64,32),dtype = complex)/norm
psi2 = np.ones((64,32),dtype = complex)/100/norm

mask = np.zeros((64,32),dtype = complex)
phase = np.zeros((64,32))
for i in range(64):
    phase[i,:] = 0.6250 * i
    for j in range(32):
        kx = (i - volres[0] / 2 + 1) / volsize[0]
        ky = (j - volres[1] / 2 + 1) / volsize[1]
        fac = -4 * np.pi**2 * hbar
        lam = fac *(kx**2 + ky**2)
        mask[i,j] = np.exp(1j * lam * dt / 2)
        
    
nozzle_cen = [2-1.7,1-0.034,1+0.066]
nozzle_len = 0.5
nozzle_rad = 0.5

n_particles = 50

for ite in range(10):
    vx = np.zeros((64,32))
    vy = np.zeros((64,32))
    for i in range(64):
        for j in range(32):
            if ((i*dx - nozzle_cen[0])**2 + (j*dy - nozzle_cen[1])**2) < nozzle_rad**2:
               amp1 = abs(psi1[i,j])
               amp2 = abs(psi2[i,j])
               psi1[i,j] = amp1 * np.exp(1j * phase[i,j])
               psi2[i,j] = amp2 * np.exp(1j * phase[i,j])
               
               ixp = int((i+1)%64)
               jxp = int((j+1)%32)
               term = np.conjugate(psi1[i,j])*psi1[ixp,j]
               vx[i,j] = math.atan2(term.imag,term.real) * hbar
               term = np.conjugate(psi2[i,j])*psi2[i,jxp]
               vy[i,j] = math.atan2(term.imag,term.real) * hbar
               
    div = np.zeros((64,32))
    for i in range(64):
        for j in range(32):
               ixm = int((i-1 + 64) % 64)
               jym = int((j-1 + 32) % 32)
               div[i,j] = (vx[i,j] - vx[ixm,j])/dx + (vy[i,j] - vy[i,jym])/dy
               
    f = np.fft.fft2(div)
    fac = np.zeros((64,32))
    for i in range(64):
        for j in range(32):
            sx = np.sin(np.pi * i) / volres[0]
            sy = np.sin(np.pi * j) / volres[1]
            denom = sx**2 + sy**2
            fac[i,j] = - 4 / denom
    fac[0,0] = 0
    f = f * fac
    p = np.fft.ifft(f)
    
    eip = np.exp(1j * p)
    psi1 = psi1 * eip
    psi2 = psi2 * eip
               
iter_max = 1
for ite in range(iter_max):
    psi1shift = np.fft.fft(psi1)
    