import numpy as np
import scipy.io as scio
import math
"""
参考论文：Matlab-Implementation of the Finite Element Method in Elasticity

本地参考代码：D:\FluidSim\FluidSim\FEMNEW\2002-AJ_CC_FS_KR-Matlab_Implementation_FEM_Elasticity\Software2\fem_lame2d

完成状态：主体部分完成，未后处理

"""
coordinates = scio.loadmat('coordinates.mat')['coordinates']
elements3 = scio.loadmat('elements3.mat')['elements3']
elements4 = scio.loadmat('elements4.mat')['elements4']
dirichlet = scio.loadmat('dirichlet.mat')['dirichlet']
elements3 -= 1
elements4 -= 1
dirichlet -= 1
Nx = coordinates.shape[0]
A = np.zeros((2*Nx,2*Nx))
b = np.zeros((2*Nx))

E = 1e5
nu = 0.3
mu = E / (2*(1 + nu))
lam = E * nu / ((1 + nu)*(1 - 2*nu))

for k in range(elements3.shape[0]):
    x0 = coordinates[elements3[k,0],0] # x0
    y0 = coordinates[elements3[k,0],1] # y0
    x1 = coordinates[elements3[k,1],0] # x1
    y1 = coordinates[elements3[k,1],1] # y1
    x2 = coordinates[elements3[k,2],0] # x2
    y2 = coordinates[elements3[k,2],1] # y2
    
    G = np.array([[1,1,1],[x0,x1,x2],[y0,y1,y2]])
    G0 = np.array([[0,0],[1,0],[0,1]])
    PhiGrad = np.dot(np.linalg.inv(G),G0)
    R = np.zeros((3,6))
    R[0,0] = R[2,1] = PhiGrad[0,0]
    R[2,0] = R[1,1] = PhiGrad[0,1]
    R[0,2] = R[2,3] = PhiGrad[1,0]
    R[2,2] = R[1,3] = PhiGrad[1,1]
    R[0,4] = R[2,5] = PhiGrad[2,0]
    R[2,4] = R[1,5] = PhiGrad[2,1]
            
    C = mu * np.array([[2,0,0],[0,2,0],[0,0,1]]) + lam * np.array([[1,1,0],[1,1,0],[0,0,0]])
    stima3 = np.zeros((6,6))
    temp = np.dot(np.transpose(R),C)
    stima3 = np.dot(temp,R)*np.det(G)/2
    
    idx = np.zeros((6),dtype = int)
    idx[0] = elements3[k,0] * 2
    idx[1] = elements3[k,0] * 2 + 1
    idx[2] = elements3[k,1] * 2
    idx[3] = elements3[k,1] * 2 + 1
    idx[4] = elements3[k,2] * 2
    idx[5] = elements3[k,2] * 2 + 1
    
    for i in range(6):
        for j in range(6):
            idxi = idx[i]
            idxj = idx[j]
            A[idxi,idxj] += stima3[i,j]
   
            
for k in range(elements4.shape[0]):
    R11 = np.array([[2,-2,-1,1],[-2,2,1,-1]
                    ,[-1,1,2,-2],[1,-1,-2,2]])/6
    R12 = np.array([[1,1,-1,-1],[-1,-1,1,1],
                    [-1,-1,1,1],[1,1,-1,-1]])/4
    R22 = np.array([[2,1,-1,-2],[1,2,-2,-1],
                    [-1,-2,2,1],[-2,-1,1,2]])/6
    
    x0 = coordinates[elements4[k,0],0] # x0
    y0 = coordinates[elements4[k,0],1] # y0
    x1 = coordinates[elements4[k,1],0] # x1
    y1 = coordinates[elements4[k,1],1] # y1
    x2 = coordinates[elements4[k,2],0] # x2
    y2 = coordinates[elements4[k,2],1] # y2
    x3 = coordinates[elements4[k,3],0] # x3
    y3 = coordinates[elements4[k,3],1] # y3
    
    F0 = np.array([[x1 - x0,y1 - y0],[x3 - x0,y3 - y0]])
    F = np.linalg.inv(F0)
    
    
    stima4 = np.zeros((8,8))
    idx0 = np.array([0,2,4,6])
    idx1 = np.array([1,3,5,7])
    
    E = np.dot(np.transpose(F),np.array([[lam+2*mu,0],[0,mu]]))
    E = np.dot(E,F)
    num = 4
    for i in range(num):
        for j in range(num):
            idxi = idx0[i]
            idxj = idx0[j]
            stima4[idxi,idxj] = E[0,0]*R11[i,j] + E[0,1]*R12[i,j] + E[1,0]*R12[j,i] + E[1,1]*R22[i,j]
    
    E = np.dot(np.transpose(F),np.array([[mu,0],[0,lam + 2*mu]]))
    E = np.dot(E,F)
    for i in range(num):
        for j in range(num):
            idxi = idx1[i]
            idxj = idx1[j]
            stima4[idxi,idxj] = E[0,0]*R11[i,j] + E[0,1]*R12[i,j] + E[1,0]*R12[j,i] + E[1,1]*R22[i,j]
            
    E = np.dot(np.transpose(F),np.array([[0,mu],[lam,0]]))
    E = np.dot(E,F)
    for i in range(num):
        for j in range(num):
            idxi = idx1[i]
            idxj = idx0[j]
            stima4[idxi,idxj] = E[0,0]*R11[i,j] + E[0,1]*R12[i,j] + E[1,0]*R12[j,i] + E[1,1]*R22[i,j]
    
    for i in range(num):
        for j in range(num):
            idxi = idx0[i]
            idxj = idx1[j]
            stima4[idxi,idxj] = stima4[idxj,idxi]
    
    idx = np.zeros((8),dtype = int)
    idx[0] = elements4[k,0] * 2
    idx[1] = elements4[k,0] * 2 + 1
    idx[2] = elements4[k,1] * 2
    idx[3] = elements4[k,1] * 2 + 1
    idx[4] = elements4[k,2] * 2
    idx[5] = elements4[k,2] * 2 + 1
    idx[6] = elements4[k,3] * 2
    idx[7] = elements4[k,3] * 2 + 1
    
    stima4 /= np.linalg.det(F)
    
    for i in range(8):
        for j in range(8):
            idxi = idx[i]
            idxj = idx[j]
            A[idxi,idxj] += stima4[i,j] 

dirichletNodes = scio.loadmat('dirichletNodes.mat')['DirichletNodes'] - 1
num = len(dirichletNodes)
coor = coordinates[dirichletNodes[:,0],:]
M = np.zeros((2 * num,2))
W = np.zeros((2 * num))
phi = np.zeros((num))
r = np.zeros((num))
for i in range(num):
    M[2 * i,0] = 1
    M[2*i+1,1] = 1
    
    phi[i] = math.atan2(coor[i,1],coor[i,0])
    r[i] = np.sqrt(coor[i,0]**2 + coor[i,1]**2)
    
alpha = 0.544483737
omega = np.pi * 3 / 4
C_1 = -np.cos((alpha+1)*omega)/np.cos((alpha-1)*omega)
C_2 = 2*(lam+2*mu)/(lam+mu)
ut = (1/(2*mu)) * r**alpha*((alpha+1)*np.sin((alpha+1)*phi)+(C_2+alpha-1)*C_1*np.sin((alpha-1)*phi))
ur = (1/(2*mu))*r**alpha* (-(alpha+1)*np.cos((alpha+1)*phi)+(C_2-(alpha+1))*C_1*np.cos((alpha-1)*phi));

value0 = ur * np.cos(phi) - ut * np.sin(phi)
value1 = ur * np.sin(phi) + ut * np.cos(phi)

for i in range(num):
    W[i * 2] = value0[i]
    W[i * 2 + 1] = value1[i]

B = np.zeros((2 * num,2 * Nx))
for i in range(num):
    idxi = i * 2
    idxj = dirichlet[i] * 2
    B[idxi,idxj] = 1
    idxi = i * 2 + 1
    idxj = dirichlet[i] * 2 + 1
    B[idxi,idxj] = 1
    
nmax = 2 * num + 2 * Nx
A1 = np.zeros((nmax,nmax))
A1[0:450,0:450] = A[:,:]
A1[450:,0:450] = B
A1[0:450,450:] = np.transpose(B)
b0 = np.zeros((nmax))
b0[0:450] = b
b0[450:] = W