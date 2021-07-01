import numpy as np
import scipy.io as scio
"""
Element Based Galerkin

Chapter08
"""
ngl = 17 # 高斯点的数量


def legendre_poly(p,x):
    L1 = 0
    L11 = 0
    L12 = 0
    L0 = 1
    L01 = 0
    L02 = 0
    for j in range(p):
        L2 = L1
        L21 = L11
        L22 = L12
        L1 = L0
        L11 = L01
        L12 = L02
        a = (2 * j + 1) / (j + 1)
        b = j / (j + 1)
        L0 = a * x * L1 - b * L2
        L01 = a * (L1 + x*L11) - b * L21
        L02 = a * (2 * L11 + x*L12) - b*L22
    return L0,L01,L02
    

def legendre_gauss_lobatto(ngl):
    p = ngl - 1
    ph = int(np.floor((p + 1)/2))
    localxgl = np.zeros((ngl))
    localwgl = np.zeros((ngl))
    for i in range(ph):
        x = np.cos((2 * i + 1) * np.pi / (2 * p + 1)) 
        for k in range(20):
            
            # Legrende Poly
            L0,L01,L02 = legendre_poly(p,x)
            
            dx =  -(1-x**2)*L01/(-2*x*L01 + (1-x**2)*L02)
            x = x + dx
            if abs(dx) < 1e-20:
                break
        localxgl[p - i] = x
        localwgl[p - i] = 2 / (p * (p + 1)*L0**2)
        
    # Check For Zero Root
    if p + 1 != 2 * ph:
        x = 0
        L0,L01,L02 = legendre_poly(p,x)
        localxgl[ph] = x
        localwgl[ph] = 2 / (p * (p + 1)*L0**2)
        
    for i in range(ph):
        localxgl[i] = -localxgl[p - i]
        localwgl[i] = localwgl[p - i]
    return localxgl,localwgl

def legendre_gauss(ngl):
    p = ngl - 1
    ph = int(np.floor((p + 1)/2))
    xgl = np.zeros((ngl))
    wgl = np.zeros((ngl))
    for i in range(ph):
        x = np.cos((2 * i + 1) * np.pi / (2 * p + 1)) 
        for k in range(20):
            
            # Legrende Poly
            L0,L01,L02 = legendre_poly(p + 1,x)
            
            dx = 2 / (p * (p + 1)*L0**2)
            x = x + dx
            if abs(dx) < 1e-20:
                break
        xgl[p - i] = x
        wgl[p - i] = 2 / ((1 - x**2)*L01**2)
        
    # Check For Zero Root
    if p + 1 != 2 * ph:
        x = 0
        L0,L01,L02 = legendre_poly(p + 1,x)
        xgl[ph] = x
        wgl[ph] = 2 / ((1 - x**2)*L01**2)
        
    for i in range(ph):
        xgl[i] = -xgl[p - i]
        wgl[i] = wgl[p - i]
    return xgl,wgl    
    

# legendre_gauss_lobatto
xgl,wgl = legendre_gauss_lobatto(ngl)
intergration_points = 1
if intergration_points == 1:
    xnq,wnq = legendre_gauss_lobatto(ngl)
elif intergration_points == 2:
    xnq,wnq = legendre_gauss(ngl)
nq = ngl
psi = np.zeros((ngl,ngl))
dpsi = np.zeros((ngl,ngl))
for m in range(nq):
    x1 = xnq[m]
    for i in range(ngl):
        xi = xgl[i]
        psi[i,m] = 1
        dpsi[i,m] = 0
        for j in range(ngl):
            xj = xgl[j]
            if i != j:
                psi[i,m] = psi[i,m] * (x1 - xj) / (xi - xj)
            ddpsi = 1
            if i != j:
                for k in range(ngl):
                    xk = xgl[k]
                    if(i != k)&(k != j):
                        ddpsi = ddpsi * (x1 - xk) / (xi - xk)
                dpsi[i,m] = dpsi[i,m] + ddpsi / (xi - xj)

# Create Grid
coord = np.zeros((17,4))
nelem = 4
xmin = -1
xmax = 1
dx = (xmax - xmin) / nelem
ip = 1
coord[0,0] = xmin
for i in range(nelem):
    x0 = xmin + i * dx
    coord[0,i] = x0
    for j in range(1,ngl):
        ip = ip + 1
        coord[j,i] = (xgl[j] + 1) * dx / 2 + x0
        
        
# Exact Solution
npoin = 68
qe = np.zeros((npoin))
qex = np.zeros((npoin))
fe = np.zeros((npoin))
icase = 1 # 1 for Homogeneous Dirichlet BCs
for e in range(nelem):
    for i in range(ngl):
        x = coord[i,e]
        ip = e * 17 + i
        if icase == 1:
            c = 2 * np.pi
            qe[ip] = np.sin(c * x)
            qex[ip] = c * np.cos(c * x)
            fe[ip] = - c**2 * np.sin(c * x)
        elif icase == 2:
            c = np.pi
            qe[ip] = np.cos(c * x) + 1
            qex[ip] = c * np.sin(c * x)
            fe[ip] = - c**2 * np.cos(c * x)
        elif icase == 3:
            c = 2 * np.pi
            qe[ip] = np.cos(c * x)
            qex[ip] = c * np.sin(c * x)
            fe[ip] = - c**2 * np.cos(c * x)
            
            
        


Mmat = np.zeros((npoin,npoin))
for e in range(nelem):
    x = np.zeros((ngl))
    for i in range(ngl):
        x[i] = coord[i,e]
    dx = x[ngl - 1] - x[0]
    jac = dx / 2
    for l in range(nq):
        wq = wnq[l] * jac
        for i in range(ngl):
            ip = e * 17 + i
            hi = psi[i,l]
            for j in range(ngl):
                jp = e * 17 + j
                hj = psi[j,l]
                Mmat[ip,jp] += wq * hi * hj
                

Lmat = np.zeros((npoin,npoin))
DhatmatQ = np.zeros((npoin,npoin))
Dhatmatq = np.zeros((npoin,npoin))
RHSvector = np.dot(Mmat,fe)

elliptic_method_type = 1
if elliptic_method_type == 1:
    for i in range(npoin):
        q = np.zeros((npoin))
        q[i] = 1
        
        