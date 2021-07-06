import numpy as np
"""
https://github.com/jnfran92/vibro-acus-fem
"""

numEle = 4
numDof = 2
numNodes = numEle + 1
numNodesEle = 2
rho = 2700
E = 7.1E10

L = 0.5
b = 0.02 # 宽度
h = 0.005 # 厚度
A = b * h # 截面积
Iz = (b * h**3)/12

# 质量和惯性矩阵
a = L / (numNodes - 1) / 2
M = np.zeros((numNodesEle * numDof))
K = np.zeros((numNodesEle * numDof))
mcons = (rho * A * a)/105
M = np.array([[78,22*a,27,-13*a],
              [22*a,8*a**2,13*a,-6*a**2],
              [27,13*a,78,-22*a],
              [-13*a,-6*a**2,-22*a,8*a**2]])*mcons
kcons = (E * Iz) / (2 * a**3)
K = np.array([[3,3*a,-3,3*a],
              [3*a,4*a**2,-3*a,2*a**2],
              [-3,-3*a,3,-3*a],
              [3*a,2*a**2,-3*a,4*a**2]])*kcons

Mg= np.zeros((numNodes * numDof))
Kg = np.zeros((numNodes * numDof))
p = numNodesEle * numDof
for e in range(numNodes-1):
    Ae = np.zeros((p,2*(e-1)+p+numNodes*numDof-p-2*(e-1)))
    for i in range(p):
        Ae[i,i+2*(e-1)] = 1
    Mg = Mg + np.dot(np.dot(np.transpose(Ae),M),Ae)
    Kg = Kg + np.dot(np.dot(np.transpose(Ae),K),Ae)
cons = np.array([0,1,8])
for i in range(len(cons)):
    Mg[:,i] = 0
    Mg[i,:] = 0
    Kg[:,i] = 0
    Kg[i,:] = 0
    
# A,B = np.linalg.eig(Kg,Mg)