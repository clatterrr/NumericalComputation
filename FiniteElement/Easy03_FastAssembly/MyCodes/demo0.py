import numpy as np
import scipy.io as scio
import math
'''
Fast MATLAB assembly of FEM matrices in 2D and 3D: nodal elements

'''
coordinates = scio.loadmat('coordinates.mat')['coordinates']
elements = scio.loadmat('elements.mat')['elements'] - 1
dirichlet = scio.loadmat('dirichlet.mat')['dirichlet'] - 1


NE = elements.shape[0]
DIM = coordinates.shape[1]
NLB = 3 # number of local basis function
coord = np.zeros((DIM,NLB,NE))
for i in range(DIM):
    for j in range(NLB):
        coord[i,j,:] = coordinates[elements[:,j],i]