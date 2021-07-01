import numpy as np
import scipy.io as scio
import math
"""
参考论文：Matlab-Implementation of the Finite Element Method in Elasticity

本地参考代码：D:\FluidSim\FluidSim\FEMNEW\2002-AJ_CC_FS_KR-Matlab_Implementation_FEM_Elasticity\Software2\fem_lame2d

完成状态：主体部分完成，未后处理
"""
coordinates = scio.loadmat('coordinates.mat')['coordinates']
elements = scio.loadmat('elements.mat')['elements']
dirichlet = scio.loadmat('dirichlet.mat')['dirichlet']
dirichletNodes = scio.loadmat('dirichletNodes.mat')['dirichletnodes']

elements -= 1
dirichlet -= 1
dirichletNodes -= 1
Nx = coordinates.shape[0]

A = np.zeros((Nx * 3,Nx * 3))
b = np.zeros((Nx * 3))

E = 1e5
nu = 0.3
mu = E / (2*(1 + nu))
lam = E * nu / ((1 + nu)*(1 - 2*nu))

for k in range(Nx):
    R11 = np.array([[2,-2,-1,1],[-2,2,1,-1]
                    ,[-1,1,2,-2],[1,-1,-2,2]])/6
    R12 = np.array([[1,1,-1,-1],[-1,-1,1,1],
                    [-1,-1,1,1],[1,1,-1,-1]])/4
    R22 = np.array([[2,1,-1,-2],[1,2,-2,-1],
                    [-1,-2,2,1],[-2,-1,1,2]])/6
    
    x0 = coordinates[elements[k,0],0]
    y0 = coordinates[elements[k,0],1]
    z0 = coordinates[elements[k,0],2]
    x1 = coordinates[elements[k,1],0]
    y1 = coordinates[elements[k,1],1]
    z1 = coordinates[elements[k,1],2]
    x2 = coordinates[elements[k,2],0]
    y2 = coordinates[elements[k,2],1]
    z2 = coordinates[elements[k,2],2]
    x3 = coordinates[elements[k,3],0]
    y3 = coordinates[elements[k,3],1]
    z3 = coordinates[elements[k,3],2]
    
    G = np.array([[1,1,1,1],[x0,x1,x2,x3],[y0,y1,y2,y3],[z0,z1,z2,z3]])
    temp0 = np.linalg.inv(G)
    PhiGrad = np.dot(temp0,np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]]))
    R = np.zeros((6,12))
    
    idx0 = np.array([0,3,4],dtype = int)
    idx1 = np.array([0,3,6,9],dtype = int)
    for i in range(3):
        for j in range(4):
            R[idx0[i],idx1[j]] = PhiGrad[j,i]
            
    idx0 = np.array([3,1,5],dtype = int)
    idx1 = np.array([1,4,7,10],dtype = int)
    for i in range(3):
        for j in range(4):
            R[idx0[i],idx1[j]] = PhiGrad[j,i]
            
    idx0 = np.array([4,5,2],dtype = int)
    idx1 = np.array([2,5,8,11],dtype = int)
    for i in range(3):
        for j in range(4):
            R[idx0[i],idx1[j]] = PhiGrad[j,i]
            
            
    C = np.array([[lam+2*mu,lam,lam,0,0,0],
                  [lam,lam+2*mu,lam,0,0,0],
                  [lam,lam,lam+2*mu,0,0,0],
                  [0,0,0,mu,0,0],
                  [0,0,0,0,mu,0],
                  [0,0,0,0,0,mu]])
    
    temp = np.dot(np.transpose(R),C)
    stima = np.dot(temp,R)*np.linalg.det(G) / 6
    idx = np.zeros((12),dtype = int)
    idx[0] = elements[k,0] * 3
    idx[1] = elements[k,0] * 3 + 1
    idx[2] = elements[k,0] * 3 + 2
    idx[3] = elements[k,0] * 3
    idx[4] = elements[k,0] * 3 + 1
    idx[5] = elements[k,0] * 3 + 2
    idx[6] = elements[k,0] * 3
    idx[7] = elements[k,0] * 3 + 1
    idx[8] = elements[k,0] * 3 + 2
    idx[9] = elements[k,0] * 3
    idx[10] = elements[k,0] * 3 + 1
    idx[11] = elements[k,0] * 3 + 2
    for i in range(12):
        for j in range(12):
            idxi = idx[i]
            idxj = idx[j]
            A[idxi,idxj] += stima[i,j]
            
coor = coordinates[dirichletNodes[:,0],:]
num = coor.shape[0]
W = np.zeros((num * 3))
for i in range(num):
    if coor[i,1] < -50:
        W[i * 3] = 0.1
        
B = np.zeros((3 * num,3 * Nx))
for i in range(num):
    idxi = i * 3
    idxj = dirichletNodes[i] * 3
    B[idxi,idxj] = 1
    B[idxi + 1,idxj + 1] = 1
    B[idxi + 2,idxj + 2] = 1
    
nmax = 3 * num + 3 * Nx
A1 = np.zeros((nmax,nmax))
A1[0:3*Nx,0:3*Nx] = A[:,:]
A1[3*Nx:,0:3*Nx] = B
A1[0:3*Nx,3*Nx:] = np.transpose(B)
b0 = np.zeros((nmax))
b0[0:3*Nx] = b
b0[3*Nx:] = W