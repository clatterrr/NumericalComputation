import numpy as np
import random
"""
D:\FluidSim\FluidSim\MultiGri\Poisson_FDM_Multigrid-master

https://github.com/EnigmaHuang/Poisson_FDM_Multigrid
"""
p = 4
n = 2**p - 1
N = n * n
A = np.zeros((N,N))
for k in range(N):
    i = int(k % n)
    j = int(k // n)
    A[k,k] = 4
    if i < n-1:
        A[k,k+1] = -1
    if i > 0:
        A[k,k-1] = -1
    if j < n-1:
        A[k,k+n] = -1
    if j > 0:
        A[k,k-n] = -1
        
b = np.zeros((N))
directn = 3 * 3
    
def Gauss(Amat,bmat,xmat):
    eps = 1e-3 # 可允许误差范围很大
    err = 1
    size = Amat.shape[0]
    while(err > eps):
        err = 0
        for i in range(size):
            summ = 0
            for j in range(size):
                summ += Amat[i,j] * xmat[j]
            xmat[i] = xmat[i] + 0.1 * (bmat[i] - summ)
            err = max(err,bmat[i] - summ)
            test = 1
        # print('err = %f',err)
    return xmat
        
def MultiGri(size,Amat,bmat,xmat):
    if size < 9:
        xmat = np.dot(np.linalg.inv(Amat),bmat)
        return xmat
    xmat = Gauss(Amat,bmat,xmat)
    fine_N = size
    fine_n = int(np.sqrt(fine_N + 1))
    coarse_n = int(np.floor((fine_n-1)/2))
    coarse_N = coarse_n * coarse_n
    R = np.zeros((coarse_N,fine_N))
    k = 0
    for j in range(1,fine_n-1,2):
        for i in range(1,fine_n-1,2):
            
            idx = j * fine_n + i
            R[k,idx-fine_n-1] = 1 / 16
            R[k,idx-fine_n] = 1 / 8
            R[k,idx-fine_n+1] = 1 / 16
            
            R[k,idx-1] = 1 / 16
            R[k,idx] = 1 / 8
            R[k,idx+1] = 1 / 16
            
            R[k,idx+fine_n-1] = 1 / 16
            R[k,idx+fine_n] = 1 / 8
            R[k,idx+fine_n+1] = 1 / 16
            
            k = k + 1
            
            
            
    coarse_A = np.dot(np.dot(R,Amat),4*np.transpose(R))
    
    P = np.transpose(R) * 4
            
    r = bmat - np.dot(Amat,xmat)
    coarse_b = np.dot(R,r)
    x0 = np.zeros((coarse_N))
    d_x = MultiGri(coarse_N, coarse_A, coarse_b,x0)
    
    xmat = xmat + np.dot(P,d_x)
    xmat = Gauss(Amat, bmat,xmat)
    
    return xmat
    
    
            
    

for i in range(n):
    b[i] = 1
    
exactAinv = np.linalg.inv(A)
exactx = np.dot(exactAinv,b)    

eps = 0.0028
norm = 1
x = np.zeros((N))
while(norm > eps):
    x = MultiGri(N,A,b,x)
    r = b - np.dot(A,x)
    norm = np.linalg.norm(r)
    print('norm = %f',norm)

