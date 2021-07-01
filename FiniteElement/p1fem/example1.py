import numpy as np
'''
Efficient Implementation of Adaptive P1-FEM in Matlab


coor = np.array([[-1,-1],
                 [0,-1],
                 [-0.5,-0.5],
                 [-1,0],
                 [0,0],
                 [1,0],
                 [-0.5,0.5],
                 [0.5,0.5],
                 [-1,1],
                 [0,1],
                 [1,1]])

elements = np.array([[1,2,3],
                     [2,5,3],
                     [5,4,3],
                     [4,1,3],
                     [4,5,7],
                     [5,10,7],
                     [10,9,7],
                     [9,4,7],
                     [5,6,8],
                     [6,11,8],
                     [11,10,8],
                     [10,5,8]])

dirichlet = np.array([[1,2],[2,5],[5,6]])
neumann = np.array([[6,11],[11,10],[10,9],[9,4],[4,1]])

dirichlet -= 1
neumann -= 1
elements -= 1

# Adaptive Refine RGB
for k in range(4):
    dx = np.zeros((12,3))
    dy = np.zeros((12,3))
    reshape = np.zeros((12,3))
    
    I = np.zeros((44),dtype = int)
    J = np.zeros((44),dtype = int)
    for i in range(12):
        for j in range(3):
            j0 = int((j + 1)%3)
            dx[i,j] = coor[elements[i,j0],0] - coor[elements[i,j],0]
            dy[i,j] = coor[elements[i,j0],1] - coor[elements[i,j],1]
            reshape[i,j] = dx[i,j]**2 + dy[i,j]**2
            
            idx = int(j * 12 + i)
            I[idx] = elements[i,j]
            J[idx] = elements[i,j0]
    
    # Provide GeometricData
    I[36:39] = dirichlet[0:3,1]
    J[36:39] = dirichlet[0:3,0]
    I[39:44] = neumann[:,1]
    J[39:44] = neumann[:,0]
    
    idxIJ = np.zeros((44),dtype = int)
    idxJI = np.zeros((44))
    edgeNumber = np.zeros((44))
    idx0 = 0
    idx1 = 0
    Number2Edges0 = np.zeros((11,11))
    Number2Edges1 = np.zeros((11,11))
    for i in range(44):
        if I[i] < J[i]:
            idxIJ[idx0] = i
            edgeNumber[i] = idx0 + 1
            Number2Edges0[I[i],J[i]] = idx0 + 1 # 防止全是零
            idx0 = idx0 + 1
        else:
            idxJI[i] = i + 1
            Number2Edges1[J[i],I[i]] = i + 1 
            
    numberingIJ = np.zeros((22))
    idx = 0
    for j in range(11):
        for i in range(11):
            if Number2Edges0[i,j] != 0:
                numberingIJ[idx] = Number2Edges0[i,j]
                idx = idx + 1
    idx = 0
    for j in range(11):
        for i in range(11):
            if Number2Edges1[i,j] != 0:
                edgeNumber[int(Number2Edges1[i,j] - 1)] = numberingIJ[idx]
                idx = idx + 1
                
    element2edges = np.zeros((12,3))
    edge2nodes = np.zeros((22,2))
    for i in range(12):
        for j in range(3):
            idx = int(j * 12 + i)
            element2edges[i,j] = edgeNumber[idx]
    for i in range(22):
        edge2nodes[i,0] = I[idxIJ[i]]
        edge2nodes[i,1] = J[idxIJ[i]]
        
    tes = 1
'''

import scipy.io as scio
coordinate = scio.loadmat('coordinatesfinal.mat')['coordinates']
element = scio.loadmat('elementsfinal.mat')['elements']
neumann = scio.loadmat('neumannfinal.mat')['neumann']
dirichlet = scio.loadmat('dirichlet.mat')['dirichlet']
            
element -= 1
neumann -= 1
dirichlet -= 1
nE = element.shape[0]
nC = coordinate.shape[0]
x = np.zeros((nC))

c1 = coordinate[element[:,0],:]
d21 = coordinate[element[:,1],:] - c1
d31 = coordinate[element[:,2],:] - c1

area4 = 2 * (d21[:,0]*d31[:,1] - d21[:,1]*d31[:,0])
I = np.zeros((nE * 9),dtype = int)
J = np.zeros((nE * 9),dtype = int)

a = np.zeros((nE))
b = np.zeros((nE))
c = np.zeros((nE))
A = np.zeros((9,nE))
Aa = np.zeros((1601,1601))
for i in range(0,nE):
    idx = i * 9
    I[idx] = I[idx + 3] = I[idx + 6] = element[i,0]
    I[idx + 1] = I[idx + 4] = I[idx + 7] = element[i,1]
    I[idx + 2] = I[idx + 5] = I[idx + 8] = element[i,2]
    
    J[idx] = J[idx + 1] = J[idx + 2] = element[i,0]
    J[idx + 3] = J[idx + 4] = J[idx + 5] = element[i,1]
    J[idx + 6] = J[idx + 7] = J[idx + 8] = element[i,2]
    
    a[i] = (d21[i,0]*d31[i,0] + d21[i,1]*d31[i,1])/area4[i]
    b[i] = (d31[i,0]*d31[i,0] + d31[i,1]*d31[i,1])/area4[i]
    c[i] = (d21[i,0]*d21[i,0] + d21[i,1]*d21[i,1])/area4[i]
    A[0,i] = -2*a[i] + b[i] + c[i]
    A[1,i] = a[i] - b[i]
    A[2,i] = a[i] - c[i]
    A[3,i] = a[i] - b[i]
    A[4,i] = b[i]
    A[5,i] = -a[i]
    A[6,i] = a[i] - c[i]
    A[7,i] = -a[i]
    A[8,i] = c[i]
    
for i in range(9*nE):
    idx = int(i % 9)
    idy = int(i // 9)
    Aa[I[i],J[i]] = A[idx,idy]
    
    
