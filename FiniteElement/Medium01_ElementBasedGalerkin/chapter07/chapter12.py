import numpy as np
import scipy.io as scio
"""
Element Based Galerkin

Chapter12
"""
ngl = 5 # 高斯点的数量


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
    
Q = 6
P = 5
# legendre_gauss_lobatto
xgl,wgl = legendre_gauss_lobatto(P)
xnq,wnq = legendre_gauss_lobatto(Q)
psi = np.zeros((P,Q))
dpsi = np.zeros((P,Q))
for l in range(Q):
    xl = xnq[l]
    for i in range(P):
        xi = xgl[i]
        psi[i,l] = 1
        dpsi[i,l] = 0
        for j in range(P):
            xj = xgl[j]
            if i != j:
                psi[i,l] = psi[i,l] * (xl - xj) / (xi - xj)
            ddpsi = 1
            if i != j:
                for k in range(P):
                    xk = xgl[k]
                    if (k != i) & (k != j):
                        ddpsi = ddpsi * (xl - xk) / (xi - xk)
                dpsi[i,l] = dpsi[i,l] + ddpsi / (xi - xj)
                
coord = np.zeros((289,2))
    
ip = 0
dx = dy = 0.5
xmin = -1
ymin = -1
for k in range(0,4):
    l1 = 1
    if k == 0:
        l1 = 0
    y0 = ymin + k * dy
    for l in range(l1,ngl):
        yy = (xgl[l] + 1) * dy / 2 + y0
        for i in range(0,4):
            j1 = 1
            if i == 0:
                j1 = 0
            x0 = xmin + i * dx
            for j in range(j1,ngl):
                xx = (xgl[j] + 1) * dx / 2 + x0
                coord[ip,0] = xx
                coord[ip,1] = yy
                ip += 1
# metrics
nelem = 16
nq = 6
ksi_x = np.zeros((nelem,nq,nq))
ksi_y = np.zeros((nelem,nq,nq))
eta_x = np.zeros((nelem,nq,nq))
eta_y = np.zeros((nelem,nq,nq))
jac = np.zeros((nelem,nq,nq))

x_ksi = np.zeros((nq,nq))
x_eta = np.zeros((nq,nq))
y_ksi = np.zeros((nq,nq))
y_eta = np.zeros((nq,nq))
x = np.zeros((ngl,ngl))
y = np.zeros((ngl,ngl))
intma = scio.loadmat('intma.mat')['intma']

testx = np.ones((5,6))
summ = 0
for i in range(5):
    for j in range(6):
        testx[i,j] = i / 6
        summ += testx[i,j] * psi[i,j]
for ie in range(nelem):
    
    for j in range(ngl):
        for i in range(ngl):
            ip = intma[ie,i,j] - 1
            x[i,j] = coord[ip,0]
            y[i,j] = coord[ip,1]
            
    for l in range(nq):
        for k in range(nq):
            for j in range(ngl):
                for i in range(ngl):
                    x_ksi[l,k] += dpsi[i,k]*psi[j,l]*x[i,j]
                    x_eta[l,k] += psi[i,k]*dpsi[j,l]*x[i,j]
                    y_ksi[l,k] += dpsi[i,k]*psi[j,l]*y[i,j]
                    y_eta[l,k] += psi[i,k]*dpsi[j,l]*y[i,j]
    
    for j in range(nq):
        for i in range(nq):
            xjac = x_ksi[i,j] * y_eta[i,j] - x_eta[i,j] * y_ksi[i,j]
            ksi_x[ie,i,j] = 1 / xjac * y_eta[i,j]
            ksi_y[ie,i,j] = -1 / xjac * x_eta[i,j]
            eta_x[ie,i,j] = -1 / xjac * y_ksi[i,j]
            eta_y[ie,i,j] = 1 / xjac * x_ksi[i,j]
            jac[ie,i,j] = wnq[i] * wnq[j] * abs(xjac)
            
npoin = 289
rhs = np.zeros((npoin,2))
inode = np.zeros((ngl,ngl))
q = np.zeros((npoin))

for np in range(npoin):
    q[:] = 0
    q[np] = 1
    for e in range(nelem):
        
        for l in range(nq):
            for k in range(nq):
                wq = jac[e,k,l]
                e_x = ksi_x[e,k,l]
                e_y = ksi_y[e,k,l]
                n_x = eta_x[e,k,l]
                n_y = eta_y[e,k,l]
                
                q_k = 0
                for j in range(ngl):
                    for i in range(ngl):
                        q_k += psi[i,k] * psi[j,l] * q[intma[ie,i,j] - 1]
                for j in range(ngl):
                    for i in range(ngl):
                        ip = intma[ie,i,j] - 1
                        h_e = dpsi[i,k] * psi[j,l]
                        h_n = psi[i,k] * dpsi[j,l]
                        dhdx = h_e * e_x + h_n * n_x
                        dhdy = h_e * e_y + h_n * n_y
                        rhs[ip,0] -= wq * dhdx * q_k
                        rhs[ip,1] -= wq * dhdy * q_k
                        
    test = 1
            