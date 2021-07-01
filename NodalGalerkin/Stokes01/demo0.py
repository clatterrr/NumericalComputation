import numpy as np
import scipy.io as scio
'''
D:\FluidSim\FluidSim\FEMNEW\HDGlab-master\Stokes

https://git.lacan.upc.edu/hybridLab/HDGlab


'''


mesh = scio.loadmat('mesh.mat')['mx']
shape0 = scio.loadmat('shape0.mat')['shape0']
shape1 = scio.loadmat('shape1.mat')['shape1']
shape2 = scio.loadmat('shape2.mat')['shape2']
Xe = np.zeros((6,2))
nOfElements = 16
nu = 1
nOfFaces = 3
tau = np.array([3,3,3])
nsd = 2
nsd2 = nsd * nsd
nOfElementNodes = 6
ndofP = 6
ndofU = 12
ndofL = 24


nOfElements = 1
for i in range(0,nOfElements):
    
    Xe = mesh[i*6:i*6+6,:]
    
    J1 = np.zeros((7,2))
    J2 = np.zeros((7,2))
    Xg = np.zeros((7,2))
    for j in range(7):
        for k in range(6):
            J1[j,:] += shape1[k,j]*Xe[k,:]
            J2[j,:] += shape2[k,j]*Xe[k,:]
            Xg[j,:] += shape0[k,j]*Xe[k,:]
    detJ = J1[:,0]*J2[:,1] - J1[:,1]*J2[:,0]
    Jinv11 = J2[:,1] / detJ
    Jinv12 = -J1[:,1] / detJ
    Jinv21 = -J2[:,0] / detJ
    Jinv22 = J1[:,0] / detJ
    dNx1 = np.zeros((6,7))
    dNx2 = np.zeros((6,7))
    for j in range(6):
        for k in range(7):
            dNx1[j,k] = shape1[j,k]*Jinv11[k] + shape2[j,k]*Jinv12[k]
            dNx2[j,k] = shape1[j,k]*Jinv21[k] + shape2[j,k]*Jinv22[k]
    
    # Gauss Points
    wXY = np.zeros((7))
    gaussWeight = np.array([0.45,0.2519,0.2519,0.2519,0.2648,0.2648,0.2648])
    Nw = np.zeros((6,7))
    for j in range(7):
        wXY[j] = gaussWeight[j] * detJ[j]
    for j in range(6):
        Nw[j,:] = shape0[j,:] * wXY[:]
    NwN = np.dot(Nw,np.transpose(shape0))
    
    All = np.zeros((24,24))
    for j in range(4):
        for k0 in range(6):
            idx0 = j + k0 * 4
            for k1 in range(6):
                idx1 = j + k1 * 4
                All[idx0,idx1] = - NwN[k0,k1]
            