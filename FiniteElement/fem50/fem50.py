import numpy as np
import math 
"""
D:\FluidSim\FluidSim\FEMNEW\fem50-master\src

https://github.com/cpraveen/fem50

Remarks around 50 lines of Matlab: short finite element implementation

"""
coordinates = np.array([[0,0],[1,0],[1.59,0],
                        [2,1],[3,1.41],[3,2],
                        [3,3],[2,3],[1,3],
                        [0,3],[0,2],[0,1],
                        [1,1],[1,2],[2,2]])

dirichlet = np.array([[3,4],[4,5],[7,8],[8,9],
                      [9,10],[10,11],[11,12],[12,1]])
dirichlet -= 1
elements3 = np.array([[2,3,13],[3,4,13],
                      [4,5,15],[5,6,15]])
elements3 -= 1
elements4 = np.array([[1,2,13,12],
                      [12,13,14,11],
                      [13,4,15,14],
                      [11,14,9,10],
                      [14,15,8,9],
                      [15,6,7,8]])
elements4 -= 1
neumann = np.array([[5,6],[6,7],[1,2],[2,3]])
neumann -= 1

Nx = coordinates.shape[0]
A = np.zeros((Nx,Nx))
b = np.zeros((Nx))

for k in range(elements3.shape[0]):
    
    x0 = coordinates[elements3[k,0],0]
    y0 = coordinates[elements3[k,0],1]
    x1 = coordinates[elements3[k,1],0]
    y1 = coordinates[elements3[k,1],1]
    x2 = coordinates[elements3[k,2],0]
    y2 = coordinates[elements3[k,2],1]
    
    G0 = np.array([[1,1,1],[x0,x1,x2],[y0,y1,y2]])
    dG = np.linalg.det(np.array([[1,1,1],[0,2,0],[0,0,1]])) # 算面积
    G1 = np.array([[0,0],[1,0],[0,1]])
    G0inv = np.linalg.inv(G0)
    G = np.dot(G0inv,G1)
    M = np.dot(np.linalg.det(G0),G)
    M = np.dot(M,np.transpose(G)) / 2
    
    for i in range(3):
        for j in range(3):
            idx0 = elements3[k,i]
            idx1 = elements3[k,j]
            A[idx0,idx1] += M[i,j]
            
for k in range(elements4.shape[0]):
    
    x0 = coordinates[elements4[k,0],0]
    y0 = coordinates[elements4[k,0],1]
    x1 = coordinates[elements4[k,1],0]
    y1 = coordinates[elements4[k,1],1]
    x2 = coordinates[elements4[k,2],0]
    y2 = coordinates[elements4[k,2],1]
    x3 = coordinates[elements4[k,3],0]
    y3 = coordinates[elements4[k,3],1]
    
    Dphi = np.array([[x1 - x0,y1 - y0],
                     [x3 - x0,y3 - y0]])
    
    B = np.linalg.inv(np.dot(Dphi,np.transpose(Dphi)))

    C1 = (np.array([[2,-2],[-2,2]])*B[0,0]
        + np.array([[3,0],[0,-3]])*B[0,1]
        + np.array([[2,1],[1,2]])*B[1,1]) 
    
    C2 = (np.array([[-1,1],[1,-1]])*B[0,0]
        + np.array([[-3,0],[0,3]])*B[0,1]
        + np.array([[-1,-2],[-2,-1]])*B[1,1]) 
    
    M = np.zeros((4,4))
    M[0:2,0:2] = M[2:4,2:4] = C1
    M[0:2,2:4] = M[2:4,0:2] = C2
    M = M * np.linalg.det(Dphi) / 6
    
    for i in range(4):
        for j in range(4):
            idx0 = elements4[k,i]
            idx1 = elements4[k,j]
            A[idx0,idx1] += M[i,j]
            
# Volume Forces
for k in range(elements3.shape[0]):
    x0 = coordinates[elements3[k,0],0]
    y0 = coordinates[elements3[k,0],1]
    x1 = coordinates[elements3[k,1],0]
    y1 = coordinates[elements3[k,1],1]
    x2 = coordinates[elements3[k,2],0]
    y2 = coordinates[elements3[k,2],1]
    
    G0 = np.array([[1,1,1],[x0,x1,x2],[y0,y1,y2]])
    
    for i in range(3):
        b[elements3[k,i]] += np.linalg.det(G0)*0
        
for k in range(elements4.shape[0]):
    x0 = coordinates[elements4[k,0],0]
    y0 = coordinates[elements4[k,0],1]
    x1 = coordinates[elements4[k,1],0]
    y1 = coordinates[elements4[k,1],1]
    x2 = coordinates[elements4[k,2],0]
    y2 = coordinates[elements4[k,2],1]
    
    G0 = np.array([[1,1,1],[x0,x1,x2],[y0,y1,y2]])
    
    for i in range(4):
        b[elements4[k,i]] += np.linalg.det(G0)*0
    
for k in range(neumann.shape[0]):
    x0 = coordinates[neumann[k,0],0]
    y0 = coordinates[neumann[k,0],1]
    x1 = coordinates[neumann[k,1],0]
    y1 = coordinates[neumann[k,1],1]
    
    norm = np.sqrt((x0-x1)**2 + (y0 - y1)**2)
    for i in range(2):
        b[neumann[k,i]] += norm * 0

u = np.zeros((Nx))
for k in range(dirichlet.shape[0]):
    x0 = coordinates[dirichlet[k,0],0]
    y0 = coordinates[dirichlet[k,0],1]
    
    theta = math.atan2(y0,x0)
    if theta < 0:
        theta = theta + 2 * np.pi
    r = np.sqrt(x0**2+y0**2)
    value = r ** (2/3) * np.sin(2 * theta / 3)
    u[dirichlet[k,0]] = value
    
    x0 = coordinates[dirichlet[k,1],0]
    y0 = coordinates[dirichlet[k,1],1]
    
    theta = math.atan2(y0,x0)
    if theta < 0:
        theta = theta + 2 * np.pi
    r = np.sqrt(x0**2+y0**2)
    value = r ** (2/3) * np.sin(2 * theta / 3)
    u[dirichlet[k,1]] = value
    
b0 = np.zeros((Nx))
for i in range(Nx):
    for j in range(Nx):
        b0[i] += (b[i] - A[i,j]*u[j])
        
freenodes = np.array([1,5,12,13,14],dtype = int)
numfree = freenodes.shape[0]
A1 = np.zeros((numfree,numfree))
b1 = np.zeros((numfree))
for i in range(numfree):
    for j in range(numfree):
        A1[i,j] = A[freenodes[i],freenodes[j]]
    b1[i] = b0[freenodes[i]]
u1 = np.dot(np.linalg.inv(A1),b1)
idx = 0
for i in range(Nx):
    if freenodes[idx] == i:
        u[i] = u1[idx]
        idx += 1
        
        