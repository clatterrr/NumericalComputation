import numpy as np
import random
"""
https://github.com/rlguy/FLIPViscosity3D
https://ending2015a.github.io/Ending2015a/52045/
"""
isize = 64
jsize = 64
ksize = 64
dx = 1 / max(max(isize,jsize),ksize)
phi = np.zeros((isize,jsize,ksize))


particleCount = 0
particlePos = np.zeros((0,3))
particleVelocity = np.zeros((0,3))

# 添加粒子,add liquid
for k in range(ksize):
    for j in range(jsize):
        for i in range(isize):
            x = random.random() * dx
            y = random.random() * dx
            z = random.random() * dx
            
            p0 = phi[i,j,k]
            p1 = phi[i+1,j,k]
            p2 = phi[i,j+1,k]
            p3 = phi[i,j,k+1]
            p4 = phi[i+1,j,k+1]
            p5 = phi[i,j+1,k+1]
            p6 = phi[i+1,j+1,k]
            p7 = phi[i+1,j+1,k+1]
            
            p8 = (p0*(1-x)*(1-y)*(1-z) + 
                  p1*x*(1-y)*(1-z) + 
                  p2*(1-x)*y*(1-z) + 
                  p3*(1-x)*(1-y)*z + 
                  p4*x*(1-y)*z + 
                  p5*(1-x)*y*z + 
                  p6*x*y*(1-z) + 
                  p7*x*y*z)
            
            if p8 < 0:
                solid_phi = p8
                if p8 >= 0:
                    particlePos.append(np.array([x,y,z]))
                    
def ComputeSignedDistanceFromParticles(radius):
    
    for pidx in range(particleCount):
        gridx = int(particlePos[pidx,0])
        gridy = int(particlePos[pidx,1])
        gridz = int(particlePos[pidx,2])
        gminx = max(gridx-1,0)
        gmaxx = min(gridx+2,isize)
        gminy = max(gridy-1,0)
        gmaxy = min(gridy+2,jsize)
        gminz = max(gridz-1,0)
        gmaxz = min(gridz+2,ksize)
        for k in range(gminz,gmaxz):
            for j in range(gminy,gmaxy):
                for i in range(gminx,gmaxx):
                    dist = np.sqrt(
                    ((i+0.5)*dx - particlePos[pidx,0])**2 + 
                    ((j+0.5)*dx - particlePos[pidx,1])**2 +
                    ((k+0.5)*dx - particlePos[pidx,2])**2)
                    if dist - radius < phi[i,j,k]:
                        phi[i,j,k] = dist
        
def getDistanceAtCellCenter():
    for k in range(ksize):
        for j in range(jsize):
            for i in range(isize):
                if phi[i,j,k] < 0.5 * dx:
                    
                    ph = (phi[i,j,k] + phi[i+1,j,k] + phi[i,j+1,k] + 
                          phi[i,j,k+1] + phi[i+1,j+1,k] + phi[i,j+1,k+1] +
                          phi[i+1,j,k+1] + phi[i+1,j+1,k+1])/8
                    if ph < 0:
                        phi[i,j,k] = -0.5 * dx
            
def ComputeVelocityScalarField(field,direction):
    weight = np.zeros((isize,jsize,ksize))
    r = dx
    rsq = r * r
    coef1 = 4 / 9 / (r**6)
    coef2 = 17 / 9 / (r**4)
    coef3 = 22 / 9 / r / r
    hdx = dx / 2
    if direction == 0:
        offset = np.array([0,hdx,hdx])
    elif direction == 1:
        offset = np.array([hdx,0,hdx])
    elif direction == 2:
        offset = np.array([hdx,hdx,0])
    
    for pidx in range(particleCount):
        p = particlePos[pidx,:] - offset
        velocityComponent = particleVelocity[pidx,direction]
        gridx = int(p[0])
        gridy = int(p[1])
        gridz = int(p[2])
        gminx = max(gridx-1,0)
        gmaxx = min(gridx+2,isize)
        gminy = max(gridy-1,0)
        gmaxy = min(gridy+2,jsize)
        gminz = max(gridz-1,0)
        gmaxz = min(gridz+2,ksize)
        for k in range(gminz,gmaxz):
            for j in range(gminy,gmaxy):
                for i in range(gminx,gmaxx):
                    v = np.array([i*dx,j*dx,k*dx]) - p
                    distsq = v * v
                    if distsq < rsq:
                        weight[i,j,k] += (1 - coef1 * distsq**2 + 
                            coef2 * distsq **2 - coef3 * distsq)
                        field[i,j,k] += weight[i,j,k] * velocityComponent
    for k in range(ksize):
        for j in range(jsize):
            for i in range(isize):
                value = field[i,j,k]
                
                if weight[i,j,k] < 1e-8:
                    continue
                
                field[i,j,k] = value / weight[i,j,k]
                
def applyBodyForce():
    _MACvelocityU = np.zeros((isize,jsize,ksize))
    _MACvelocityV = np.zeros((isize,jsize,ksize))
    _MACvelocityW = np.zeros((isize,jsize,ksize))
    _liquidSDF = np.zeros((isize,jsize,ksize))
    dt = 0
    for k in range(ksize):
        for j in range(jsize):
            for i in range(isize):
                flag = False
                if i == 0:
                    flag = _liquidSDF[i,j,k] == True
                elif i == isize - 1:
                    flag = _liquidSDF[i-1,j,k] == True
                else:
                    flag = (_liquidSDF[i,j,k] == True) | (_liquidSDF[i-1,j,k] == True)
                if flag == True:
                    _MACvelocityU[i,j,k] = 0 * dt
                
                flag = False
                if j == 0:
                    flag = _liquidSDF[i,j,k] == True
                elif j == jsize - 1:
                    flag = _liquidSDF[i,j-1,k] == True
                else:
                    flag = (_liquidSDF[i,j,k] == True) | (_liquidSDF[i,j-1,k] == True)
                if flag == True:
                    _MACvelocityV[i,j,k] = 0 * dt
                    
                flag = False
                if k == 0:
                    flag = _liquidSDF[i,j,k] == True
                elif k == ksize - 1:
                    flag = _liquidSDF[i,j,k-1] == True
                else:
                    flag = (_liquidSDF[i,j,k] == True) | (_liquidSDF[i,j,k-1] == True)
                if flag == True:
                    _MACvelocityW[i,j,k] = 0 * dt
                    
# Incomplete factorization
def IncompleteFactorization(A):
    matrixlen = A.shape[0]
    
    # L 是 A 的下三角矩阵
    L = np.zeros((matrixlen,matrixlen))
    for i in range(matrixlen):
        for j in range(i+1):
            L[i,j] = A[i,j]
    
    omega = 1
    for i in range(matrixlen):
        L[i,i] = np.sqrt(L[i,i])
        L[i+1:,i] = L[i+1,i] / L[i,i]
        for j in range(matrixlen):
            if (L[j,i] > 0) & (j > i):
                fullupdate = L[:,i] * L[j,i]
                incompleteupdate = np.zeros((matrixlen))
                for k in range(matrixlen):
                    if A[k,j] != 0:
                        incompleteupdate[k] = fullupdate[k]
                missing = sum(fullupdate - incompleteupdate)
                L[j:,j] = L[j:,j] - incompleteupdate[j:]
                L[j,j] = L[j,j] - omega * missing
                

def fractionInside2(phiLeft,phiRight):
    if (phiLeft < 0) & (phiRight < 0):
        return 1
    if (phiLeft < 0) & (phiRight >= 0):
        return phiLeft / (phiLeft - phiRight)
    if (phiLeft >= 0) & (phiRight < 0):
        return phiRight / (phiRight - phiLeft)
    return 0

def fractionInside4(phibl,phibr,phitl,phitr):
    insideCount = (phibl < 0) + (phibr < 0) + (phitl < 0) + (phitr < 0)
    if insideCount == 4:
        return 1
    elif insideCount == 3:
        # Rotate Cycle 没写
        side0 = 1 - fractionInside2(phibl, phitl)
        side1 = 1 - fractionInside2(phibl, phibr)
        return 1 - side0 * side1 / 2
    elif insideCount == 2:
        return 1
    elif insideCount == 1:
        side0 = 1 - fractionInside2(phibl, phitl)
        side1 = 1 - fractionInside2(phibl, phibr)
        return side0 * side1 / 2
    else:
        return 0
                
def ComputeWeight():
    _weightGridU = np.zeros((isize,jsize,ksize))
    for k in range(ksize):
        for j in range(jsize):
            for i in range(isize):
                weight = 1 - fractionInside4(phi[i,j,k],phi[i,j+1,k],phi[i,j,k+1],phi[i,j+1,k+1])
                if weight < 0:
                    weight = 0
                if weight > 1:
                    weight = 1
                _weightGridU[i,j,k] = weight

def CalculateMatrixCoefficients():
    n = 1
    Adiag = np.zeros((n,n))
    scale = 1
    Aplusi = np.zeros((n,n))
    Aplusj = np.zeros((n,n))
    phi = np.zeros((n,n))
    _weightGridU = np.zeros((isize,jsize,ksize))
    for k in range(ksize):
        for j in range(jsize):
            for i in range(isize):
                # if it is pressure grid
                term = _weightGridU[i,j,k]
                phiRight = phi[i+1,j,k]
                if phiRight < 0:
                    Adiag[i,j,k] += term
                    Aplusi[i,j,k] -= term
                else:
                    theta = max(fractionInside2(phi[i,j,k],phi[i+1,j,k]),0.01)
                    Adiag[i,j,k] += term / theta
                
                term = _weightGridU[i,j,k]
                phiLeft = phi[i-1,j,k]
                if phiRight < 0:
                    Adiag[i,j,k] += term
                else:
                    theta = max(fractionInside2(phi[i-1,j,k],phi[i,j,k]),0.01)
                    Adiag[i,j,k] += term / theta

def CalculatePreconditionerVector():
    tau = 0.97
    sigma = 0.25
    
time = 0
timeFinal = 1
while(time < timeFinal):
    
    # UpdateSignedSDF
    ComputeSignedDistanceFromParticles()
    getDistanceAtCellCenter()
    
    applyBodyForce()
    
    # applyVisocity
    
    # Project
    
    time = 10
    