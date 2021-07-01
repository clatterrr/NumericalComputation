import numpy as np
"""

Element Based Galerkin

Chapter07

"""
nelem = 20
coor = np.zeros((2,nelem))
for i in range(nelem):
    coor[0,i] = -1 + 0.1*i
    coor[1,i] = -0.9 + 0.1*i 

ngl = 2 # 高斯积分点的数量
nq = 3
wnq = np.array([1/3,4/3,1/3])

# Compute Lagrange Polynomial and Derivatives
psi = np.array([[1,0.5,0],[0,0.5,1]])
dpsi = np.array([[-0.5,-0.5,-0.5],[0.5,0.5,0.5]])

mass = np.zeros((ngl,ngl,nelem))
for ie in range(nelem):
    x = coor[:,ie]
    dx = x[ngl - 1] - x[0]
    jac = dx / 2
    for k in range(nq):
        wq = wnq[k] * jac
        for i in range(ngl):
            hi = psi[i,k]
            for j in range(ngl):
                hj = psi[j,k]
                mass[i,j,ie] += wq * hi * hj
                
diff = np.zeros((ngl,ngl))
for i in range(ngl):
    for j in range(ngl):
        for k in range(nq):
            wk = wnq[k]
            dhdx = dpsi[i,k]
            hjk = psi[j,k]
            diff[i,j] += wk * dhdx * hjk
            
laplacian = np.zeros((ngl,ngl,nelem))
for ie in range(nelem):
    x = coor[:,ie]
    dx = x[ngl - 1] - x[0]
    jac = dx / 2
    dksidx = 2 / dx
    
    # LGL Intergration
    for i in range(ngl):
        for j in range(ngl):
            for k in range(nq):
                wq = wnq[k] * jac
                dhdx_ik = dpsi[i,k]*dksidx
                dhdx_jk = dpsi[j,k]*dksidx
                laplacian[i,j,ie] += wq * dhdx_ik * dhdx_jk
                
npoin = 21
Mmat = np.zeros((npoin,npoin))
Dmat = np.zeros((npoin,npoin))
Fmat = np.zeros((npoin,npoin))
Lmat = np.zeros((npoin,npoin))
inode = np.zeros((2,20),dtype = int)
periodic = np.zeros((21),dtype = int)
for i in range(20):
    inode[0,i] = i
    inode[1,i] = i + 1
    periodic[i] = i

u = 2
visc = 0.01
for ie in range(nelem):
    for i in range(ngl):
        ip = periodic[inode[i,ie]]
        for j in range(ngl):
            jp = periodic[inode[j,ie]]
            Mmat[ip,jp] += mass[i,j,ie]
            Dmat[ip,jp] += u * diff[i,j]
            Lmat[ip,jp] += laplacian[i,j,ie]
Mmat[npoin-1,npoin-1] = 1 # Countinues Galerkin
Lmat_hat = np.dot(np.linalg.inv(Mmat),Lmat)
HVmat = Lmat_hat.copy()
Rmat = np.dot(np.linalg.inv(Mmat),Dmat - Fmat)
Dmat_hat = Rmat - visc * HVmat

q0 = np.zeros((21))
q0[8] = q0[9] = q0[10] = 1
tmax = 33
for t in range(tmax):
    
    dtt = 0.0076
    # RK4
    