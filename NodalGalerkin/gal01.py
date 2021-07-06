import numpy as np
"""
https://github.com/kaushikcfd/Discontinuous-Galerkin
"""

N = 5
def LegendrePolynomial():
    LegPoly = np.zeros((N+1))
    LegPoly_1 = np.zeros((N+1))
    LegPoly_2 = np.zeros((N+1))
    temp = np.zeros((N+1))
    
    LegPoly_1[0] = 0
    LegPoly_1[1] = 1
    LegPoly_2[0] = 1
    
    if N == 0:
        LegPoly[0] = 1
        return 
    elif N == 1:
        LegPoly[0] = 0
        LegPoly[1] = 1
        return 
    for i in range(2,N+1):
        for j in range(1,i+1):
            LegPoly[j] = ((2*i - 1)/i)*LegPoly_1[j-1]
        
        LegPoly[0] = -((i-1)/i)*LegPoly_2[0]
        temp[0:i+1] = LegPoly_1[0:i+1]
        LegPoly_1[0:i+2] = LegPoly[0:i+2]
        LegPoly_2[0:i+1] = temp[0:i+1]
    
    t = 1

def polyDeriv(GivenPoly):
    DerivedPoly = np.zeros((N + 1))    
    for i in range(1,N+1):
        DerivedPoly[i] = i * GivenPoly[i]

def lobattoNodes():
    LegendrePolynomial()
    polyDeriv()
    Nodes = np.zeros((N+1))
    for i in range(1,N-1):
        
    t = 1


def setDomain():
    lobattoNodes()
    t = 1