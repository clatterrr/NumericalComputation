import numpy as np

numNodes = 891
coord = np.zeros((numNodes,2))
K = np.zeros((numNodes*2,numNodes*2))
for i in range(numNodes):
    coord[i,0] = int(i / 11) / 10
    coord[i,1] = i % 11 / 10
    
gaussxi = np.array([-1,1,1,-1]) / np.sqrt(3)
gausseta = np.array([-1,-1,1,1]) / np.sqrt(3)

numEle = 800
E = 1e5
nu = 0.25
D = E / (1 - nu**2)*np.array([[1,nu,0],[nu,1,0],[0,0,(1-nu)/2]])

for e in range(numEle):
    Ke = np.zeros((8,8))
    start = int((e // 10) * 11 + e % 10)
    xe = np.zeros((4,2))
    xe[0,:] = coord[start,:]
    xe[1,:] = coord[start+11,:]
    xe[2,:] = coord[start+12,:]
    xe[3,:] = coord[start+1,:]
    for igauss in range(4):
        xi = gaussxi[igauss]
        eta = gausseta[igauss]
        Jpar = np.zeros((2,4))
        Jpar[0,0] = -(1 - eta)/4
        Jpar[0,1] = (1 - eta)/4
        Jpar[0,2] = (1 + eta)/4
        Jpar[0,3] = -(1 + eta)/4
        
        Jpar[1,0] = -(1 - xi)/4
        Jpar[1,1] = -(1 + xi)/4
        Jpar[1,2] = (1 + xi)/4
        Jpar[1,3] = (1 - xi)/4
        
        J = np.dot(Jpar,xe)
        Jinv = np.linalg.inv(J)
        Npar = np.zeros((2,4))
        B = np.zeros((3,8))
        for i in range(4):
            
            Npar[0,i] = Jinv[0,0] * Jpar[0,i] + Jinv[0,1] * Jpar[1,i]  
            Npar[1,i] = Jinv[1,0] * Jpar[0,i] + Jinv[1,1] * Jpar[1,i]
            B[0,2*i] = Npar[0,i]
            B[1,2*i+1] = Npar[1,i]
            B[2,2*i] = Npar[1,i]
            B[2,2*i+1] = Npar[0,i]
            
        temp = np.dot(np.transpose(B),D)
        detJ = np.linalg.det(J)
        Ke = Ke + np.dot(temp,B)*detJ