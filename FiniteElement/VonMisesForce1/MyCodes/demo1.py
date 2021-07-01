import numpy as np
import scipy.io as scio
'''
Efficient and flexible MATLAB implementation of 2D and 3D elastoplastic problems

https://github.com/matlabfem/matlab_fem_elastoplasticity

D:\FluidSim\FluidSim\FEMNEW\matlab_fem_elastoplasticity-master

'''

young = 206900
poisson = 0.29
shear = young / (2*(1 + poisson))
bulk = young / (3*(1 - 2*poisson))
# constant volume forces in each directions
volumeForce = np.array([0,-1])
# constant tranction on the back side of the body in each direction
tractionForce = np.array([0,450])
size_xy = 10
size_hole = 5

nx = 96
coord = np.zeros((nx,2))
for i in range(nx):
    if i < 30:
        coord[i,0] = i % 6 + 5
        coord[i,1] = i // 6
    else:
        coord[i,0] = (i - 30) % 11
        coord[i,1] = (i - 30) // 11 + 5
    

elem = scio.loadmat('elem.mat')['ELEM']
surf = scio.loadmat('surf.mat')['SURF']
neumann = scio.loadmat('neumann.mat')['NEUMANN']
dirichlet = scio.loadmat('dirichlet.mat')['DIRICHLET']


def quadrature_volume(elementType):
    if elementType == 1:
        # 1-point quadrature rule
        Xi = np.array([1/3,1/3])
        WF = 1/2
    elif elementType == 2:
        # 7-point quadrature rule
        Xi=np.array([[0.1012865073235, 0.7974269853531, 0.1012865073235,
                      0.4701420641051, 0.4701420641051, 0.0597158717898, 1/3]
                     ,[0.1012865073235, 0.1012865073235, 0.7974269853531, 
                       0.0597158717898, 0.4701420641051, 0.4701420641051, 1/3]])
        WF=np.array([0.1259391805448, 0.1259391805448, 0.1259391805448,
                     0.1323941527885, 0.1323941527885, 0.1323941527885, 0.225])/2
    elif elementType == 3:
        pt = 1 / np.sqrt(3)
        Xi = np.array([[-pt,-pt,pt,pt],[-pt,pt,-pt,pt]])
        WF = np.array([1,1,1,1])
    elif elementType == 4: 
        pt = np.sqrt(3/5)
        Xi = np.array([[-pt,pt,pt,-pt,0,pt,0,-pt,0],
                       [-pt,-pt,pt,pt,-pt,0,pt,0,0]])
        WF = np.array([25/81,25/81,25/81,25/81,40/81,40/81,40/81,40/81,64/81])
        
    return Xi,WF
        
def quadrature_surface(elementType):
    if elementType == 1:
        # 1-point quadrature rule
        Xi_s = 0
        WF_s = 2
    elif elementType == 2:
        pt = 1 / np.sqrt(3)
        Xi_s = np.array([-pt,pt])
        WF_s = np.array([1,1])
    return Xi_s,WF_s
        
def local_basis_volume(elementType,Xi):
    xi1 = Xi[0]
    xi2 = Xi[1]
    if elementType == 1:
        HatP = np.array([[1-xi1-xi2],[xi1],[xi2]])
        DHatP1 = np.array([[-1],[1],[0]])
        DHatP2 = np.array([[-1],[0],[1]])
    elif elementType == 2:
        xi0 = 1 - xi1 - xi2
        nq = len(xi1)
        HatP = np.array([[xi0*(2*xi0 - 1)],[xi1*(2*xi1 - 1)],[xi2*(2*xi2 - 1)],
                         [4*xi1*xi2],[4*xi0*xi2],[4*xi0*xi1]])
        DHatP1 = np.array([[-4*xi0+1],[-4*xi1-1],[0],[4*xi2],[-4*xi2],[4*(xi0-xi1)]])
        DHatP2 = np.array([[-4*xi0+1],[0],[4*xi2-1],[4*xi1],[4*(xi0-xi2)],[-4*xi1]])
    return HatP,DHatP1,DHatP2
        
def local_basis_surface(elementType,Xi_s):
    xi = Xi_s
    if (elementType == 1) | (elementType == 3):
        HatP_s = np.array([[1-xi],[1+xi]])/2
        DHatP1_s = np.array([[-1/2],[1/2]])
    elif (elementType == 2) | (elementType == 4):
        HatP_s = np.array([[xi*(xi-1)/2],
                           [xi*(xi+1)/2],
                           [(xi+1)*(1-xi)]])
        DHatP1_s = np.array([[xi-1/2],
                             [xi+1/2],
                             [-2*xi]])
    return HatP_s,DHatP1_s
    
Xi_s,WF_s = quadrature_surface(1)
Xi,WFf = quadrature_volume(1)
HatP,DHatP1,DHatP2 = local_basis_volume(1,Xi)
HatP_s,DHatP1_s = local_basis_surface(1, Xi_s)

n_e = elem.shape[1] # 元素的数量
n_q = 1 # number of quadratic points
n_int = n_e * n_q # total number of integration points

shearMat = shear * np.ones((n_e))
bulkMat = bulk * np.ones((n_int))

J11 = np.zeros((n_int))
J12 = np.zeros((n_int))
J21 = np.zeros((n_int))
J22 = np.zeros((n_int))

coordx = np.zeros((3,n_int))
coordy = np.zeros((3,n_int))
for i in range(n_int):
    coordx[:,i] = coord[elem[:,i]-1,0]
    coordy[:,i] = coord[elem[:,i]-1,1]
    J11[i] = np.dot(coordx[:,i],DHatP1[:])
    J12[i] = np.dot(coordy[:,i],DHatP1[:])
    J21[i] = np.dot(coordx[:,i],DHatP2[:])
    J22[i] = np.dot(coordy[:,i],DHatP2[:])

Det = J11 * J22 - J12 * J21
Jinv11 = J22 / Det
Jinv12 = -J12 / Det
Jinv21 = -J21 / Det
Jinv22 = J11 / Det

Dphi1 = np.zeros((3,n_int))
Dphi2 = np.zeros((3,n_int))
for i in range(n_int):
    Dphi1[:,i] = Jinv11[i] * DHatP1[:,0] + Jinv12[i] * DHatP2[:,0]
    Dphi2[:,i] = Jinv21[i] * DHatP1[:,0] + Jinv22[i] * DHatP2[:,0]

n_b = 6 * elem.shape[0]
# strain-displacement matrix B
vB = np.zeros((n_b,n_int))
for i in range(n_int):
    vB[0,i] = Dphi1[0,i]
    vB[6,i] = Dphi1[1,i]
    vB[12,i] = Dphi1[2,i]
    
    vB[5,i] = Dphi1[0,i]
    vB[11,i] = Dphi1[1,i]
    vB[17,i] = Dphi1[2,i]
    
    vB[4,i] = Dphi2[0,i]
    vB[10,i] = Dphi2[1,i]
    vB[16,i] = Dphi2[2,i]
    
    vB[2,i] = Dphi2[0,i]
    vB[8,i] = Dphi2[1,i]
    vB[14,i] = Dphi2[2,i]
    
# 96 个点
# 150 个元素
B = np.zeros((450,192))
for i in range(150):
    
    idyb = np.zeros((6))
    idyb[0] = elem[0,i]*2 - 2 # 顶点0的x轴位移编号
    idyb[1] = elem[0,i]*2 - 1 # 顶点0的y轴位移编号
    idyb[2] = elem[1,i]*2 - 2
    idyb[3] = elem[1,i]*2 - 1
    idyb[4] = elem[2,i]*2 - 2
    idyb[5] = elem[2,i]*2 - 1
    for j in range(18):
        
        idx = int(i * 3 + j % 3)
        idy = int(idyb[j // 3])
        B[idx,idy] = vB[j,i]
  
Vol = np.array([[1,1,0],[1,1,0],[0,0,0]]) # lambda
Dia = np.array([[1,0,0],[0,1,0],[0,0,1/2]]) # mu
Dev = Dia - Vol / 3
Elast = np.zeros((9,n_int))
D = np.zeros((450,450))
for i in range(150):
    for j in range(3):
        for k in range(3):
            idx = int(j * 3 + k)
            Elast[idx,i] = 2 * Dev[j,k] * shearMat[i] + Vol[j,k] * bulkMat[i]
for i in range(150):
    for j in range(9):
        idx = int(i * 3 + j % 3)
        idy = int(i * 3 + j // 3)
        D[idx,idy] = Elast[j,i]/2
        
k0 = np.dot(np.transpose(B),D)
K = np.dot(k0,B)