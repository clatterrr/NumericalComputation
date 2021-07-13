import numpy as np
# http://www.netlib.org/templates/matlab/
nx = 8
A = np.zeros((nx*nx,nx*nx))

def inDomain(i,j):
    if i < 0 or j < 0 or i >= nx  or j >= nx:
        return 0
    return 1
div = np.zeros((nx,nx))
b = np.zeros((nx * nx))
div[3,1] = 1
div[4,1] = 1
for i in range(nx):
    for j in range(nx):
        idx = j * nx + i
        b[idx] = div[i,j]
for k in range(nx*nx):
    A[k,k] = 4
    i = int(k % nx)
    j = int(k / nx)
    if inDomain(i-1, j):
        A[k,k-1] = -1
    if inDomain(i+1, j):
        A[k,k+1] = -1
    if inDomain(i, j-1):
        A[k,k-nx] = -1
    if inDomain(i, j+1):
        A[k,k+nx] = -1
x = np.zeros((nx * nx))
p = np.zeros((nx * nx)) # 方向
res = b - np.dot(A,x) # 残差
bnrm2 = np.linalg.norm(b)
tmax = 100
restld = res.copy()
rho1 = 1
for t in range(0,tmax):
    z = res # 预处理方程组的残差
    rho = np.dot(np.transpose(res),z)
    if rho == 0:
        break
    beta = 0
    if t > 0:
        beta = rho / rho1
    p = z + beta*p
    Ap = np.dot(A,p)
    alpha = rho / np.dot(np.transpose(p),Ap)
    x = x + alpha * p
    res = res - alpha * Ap
    error = np.linalg.norm(res) / bnrm2
    if error < 1e-8:
        break
    rho1 = rho
    
    
pre = np.zeros((nx,nx))
for i in range(nx):
    for j in range(nx):
        idx = j * nx + i
        pre[i,j] = x[idx]
        
    
    
    