import numpy as np
import random
"""
D:\FluidSim\FluidSim\MultiGri\Poisson_FDM_Multigrid-master

https://github.com/EnigmaHuang/Poisson_FDM_Multigrid
"""
p = 5
n = 2**p - 1
A = np.zeros((n,n))
A[0,0] = 2
A[0,1] = -1
for i in range(1,n-1):
    A[i,i-1] = -1
    A[i,i] = 2
    A[i,i+1] = -1
A[n-1,n-1] = 2
A[n-1,n-2] = -1
b = np.zeros((n))
directn = 16

def VCycle(Amat,directn):
    
    while(n > directn):
        coarse_n = int(np.floor((n-1)/2))
        R = np.zeros((coarse_n,n))
        for i in range(coarse_n):
            col = 2 * i
            R[i,col] = 1 / 4
            R[i,col+1] = 1 / 2
            R[i,col+2] = 1 / 4
    
def Gauss(Amat,bmat,xmat):
    eps = 1e-2 # 可允许误差范围很大
    err = 1
    size = Amat.shape[0]
    while(err > eps):
        err = 0
        for i in range(size):
            summ = 0
            for j in range(size):
                summ += Amat[i,j] * xmat[j]
            xmat[i] = xmat[i] + 0.5 * (bmat[i] - summ)
            err = max(err,bmat[i] - summ)
            test = 1
    return xmat
        
def MultiGri(size,Amat,bmat,xmat):
    if size < 16:
        xmat = np.dot(np.linalg.inv(Amat),bmat)
        return xmat
    xmat = Gauss(Amat, bmat,xmat)
    fine_n = size
    coarse_n = int(np.floor((fine_n-1)/2))
    
    
    R = np.zeros((coarse_n,fine_n))
    for i in range(coarse_n):
            col = 2 * i
            R[i,col] = 1 / 4
            R[i,col+1] = 1 / 2
            R[i,col+2] = 1 / 4
            
    coarse_A = np.dot(np.dot(R,A),2*np.transpose(R))
    
    P = np.transpose(R) * 2
            
    r = bmat - np.dot(Amat,xmat)
    coarse_b = np.dot(R,r)
    x0 = np.zeros((coarse_n))
    d_x = MultiGri(coarse_n, coarse_A, coarse_b,x0)
    
    xmat = xmat + np.dot(P,d_x)
    xmat = Gauss(Amat, bmat,xmat)
    
    return xmat
    
    
            
    

for i in range(n):
    b[i] = 1
    
exactAinv = np.linalg.inv(A)
exactx = np.dot(exactAinv,b)    

eps = 1e-6
norm = 1
x = np.zeros((n))
while(norm > eps):
    
    x = MultiGri(n,A,b,x)
    
    
    r = b - np.dot(A,x)
    norm = np.linalg.norm(r)
    print(norm)

