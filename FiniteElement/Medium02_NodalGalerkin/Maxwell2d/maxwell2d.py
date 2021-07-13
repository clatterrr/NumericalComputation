import numpy as np
import scipy.io as scio
from scipy.special import gamma
import math
import matplotlib.pyplot as plt
"""
D:\FluidSim\FluidSim\CompressibeNewgood\nodal-dg-master\nodal-dg-master
"""
Ez = scio.loadmat('Ez.mat')['Ez']
Np = 66
K = 146
Hx = np.zeros((Np,K))
Hy = np.zeros((Np,K))
resHx = np.zeros((Np,K))
resHy = np.zeros((Np,K))
resEz = np.zeros((Np,K))
N = 2


def JacobiGQ(alpha,beta,N):
    J3 = np.zeros((N))
    J = np.zeros((N + 1,N + 1))
    h1 = np.linspace(0,N,N+1)*2 + alpha + beta
    for i in range(N):
        J3[i] = np.sqrt((i + 1)*(i + 1 + beta)*(i + 1 + alpha
                )*(i + 1 + beta + alpha)/(h1[i] + 1)/(h1[i] + 3))
        J[i,i+1] = J[i+1,i] = 2 / (h1[i] + 2) * J3[i]
    D,V = np.linalg.eig(J)
    w0 = np.transpose(V[0,:])**2 * 2**(alpha + beta + 1)
    w1 = gamma(alpha + 1) * gamma(beta + 1)
    w2 = gamma(alpha + beta + 1)*(alpha + beta + 1)
    w = w0 * w1 / w2
    return w,D

def JacobiGL(alpha,beta,N):
    w,xint = JacobiGQ(alpha + 1, beta + 1, N - 2)
    localx = np.zeros((xint.shape[0] + 2))
    localx[0] = -1
    localx[-1] = 1
    localx[1:-1] = xint[:]
    return localx
    
def Vandermonde1D(N,r):
    VID = np.zeros((len(r),N+1))
    for j in range(N+1):
        VID[:,j] = JacobiP(r[:], 0,0,j)

def JacobiP(x,alpha,beta,N):
    PL = np.zeros((N+1,len(x)))
    gamma0 = 2**(alpha + beta + 1)/(alpha + beta + 1)/gamma(alpha + beta + 1)
    gamma0 = gamma0 * gamma(alpha + 1) * gamma(beta + 1)
    PL[0,:] = 1 / np.sqrt(gamma0)
    gamma1 = (alpha + 1)*(beta + 1)/(alpha + beta + 3)*gamma0
    PL[1,:] = ((alpha + beta + 2)*x/2 + (alpha - beta)/2)/np.sqrt(gamma1)
    aold = 2 / (2 + alpha + beta)*np.sqrt((alpha + 1)*(beta + 1)/(alpha + beta + 3))
    for i in range(N-1):
        hi = 2 * i + alpha + beta + 2
        anew = 2 / (hi + 2) * np.sqrt((i+2)*(i+2+alpha+beta)
                *(i+2+alpha)*(i+2+beta)/(hi+1)/(hi+3))
        bnew = -(alpha**2 - beta**2)/ hi /(hi + 2)
        PL[i+2,:] = 1 / anew*(-aold*PL[i,:] + (x - bnew)*PL[i+1,:])
        aold = anew
    

def Warpfactor(N,rout):
    LGLr = JacobiGL(0, 0, N)
    req = np.linspace(-1,1,N+1)
    
    return req

waste,rLGL = JacobiGQ(0,0,2)
rmin = abs(rLGL[0] - rLGL[1])

N = 10
alpopt = np.array([0,0,1.4152,0.1001,0.2751,0.9800,1.0999,1.2832,
                   1.3648,1.4773,1.4959,1.5743,1.5700,1.6223,1.6253])
alpha = alpopt[N - 1]
Np = int((N + 1) * (N + 2) / 2)
L1 = np.zeros((Np))
L2 = np.zeros((Np))
L3 = np.zeros((Np))
sk = 0
for n in range(N+1):
    for m in range(N+1-n):
        L1[sk] = (n - 1) / N
        L3[sk] = (m - 1) / N
        sk = sk + 1
L2 = 1 - L1 - L2
x = -L2 + L3
y = (-L2 - L3 + 2*L1) / np.sqrt(3)
blend1 = 4 * L2 * L3
blend2 = 4 * L1 * L3
blend3 = 4 * L1 * L2

warpf1 = Warpfactor(N, L3 - L2)
        

whatx = scio.loadmat('whatx.mat')['x']
whaty = scio.loadmat('whaty.mat')['y']
vx1 = whatx[0,:]
vx2 = whatx[10,:]
vx3 = whatx[65,:]
vy1 = whaty[0,:]
vy2 = whaty[10,:]
vy3 = whaty[65,:]

# Semi-perimeter and area
len1 = np.sqrt((vx1 - vx2)**2 + (vy1 - vy2)**2)
len2 = np.sqrt((vx2 - vx3)**2 + (vy2 - vy3)**2)
len3 = np.sqrt((vx3 - vx1)**2 + (vy3 - vy1)**2)
sper = (len1 + len2 + len3) / 2
Area = np.sqrt(sper * (sper - len1) * (sper - len2) * (sper - len3))
# Compute scale using radius of inscribed circle
dtscale = Area / sper
dt = min(dtscale) * rmin * 2 / 3
time = 0
Finaltime = 1
rk4a = np.array([0,-0.4179,-1.1922,-1.6978,-1.5142])
rk4b = np.array([0.1497,0.3792,0.8230,0.6995,0.1531])
Nfp = 11
Nfaces = 3
K = 146
vmapM = scio.loadmat('vmapM.mat')['vmapM']
vmapP = scio.loadmat('vmapP.mat')['vmapP']

nx = np.zeros((33,146))
ny = np.zeros((33,146))

while(time < Finaltime):
    
    time = time + dt
    for intrk in range(5):
        dHx = np.zeros((Nfp * Nfaces,K))
        dHy = np.zeros((Nfp * Nfaces,K))
        dEz = np.zeros((Nfp * Nfaces,K))
        
        for i in range(66):
            for j in range(146):
                idx0 = int(vmapM[j * 66 + i] % 66)
                idy0 = int(vmapM[j * 66 + i] // 66)
                idx1 = int(vmapP[j * 66 + i] % 66)
                idy1 = int(vmapP[j * 66 + i] // 66)
                dHx[i,j] = Hx[idx0,idy0] - Hx[idx1,idy1]
                dHy[i,j] = Hy[idx0,idy0] - Hy[idx1,idy1]
                dEz[i,j] = Ez[idx0,idy0] - Ez[idx1,idy1]
                
                alpha = 1
                ndotdH = nx * dHx + ny * dHy
                fluxHx = ny * dEz + alpha * (ndotdH * nx - dHx)
                fluxHy = - nx * dEz + alpha * (ndotdH * ny - dHy)
                fluxEz = - nx * dHy + ny * dHx - alpha * dEz
        
        