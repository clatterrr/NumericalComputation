import numpy as np

sizeX = 128
sizeY = 128
densityAir = 0.1
densitySoot = 0.25
diffusion = 0.01
timestep = 0.0025

phi = np.zeros((sizeX+1,sizeY+1))
volume = np.ones((sizeX,sizeY))
normalX = np.zeros((sizeX,sizeY))
normalY = np.zeros((sizeX,sizeY))
cell = np.zeros((sizeX,sizeY))
body = np.zeros((sizeX,sizeY))
mask = np.zeors((sizeX,sizeY))

particleCount = 64
# 0 for d,1 for t,2 for u,3 for v
quantities = np.zeros((4,particleCount)) 

def gridToParticles(alpha):
    for i in range(quantities.shape[0]):
        for j in range(particleCount):
            quantities[i,j] *= 1 - alpha
            # quantities[i,j] = posX[j] + posY[j]
            
def buildHeatDiffusionMatrix(timestep):
    _aDiag = np.zeros((sizeX,sizeY))
    for i in range(sizeX):
        for j in range(sizeY):
            _aDiag[i,j] = 1
    _aPlusX = np.zeros((sizeX,sizeY))
    _aPlusY = np.zeros((sizeX,sizeY))
    cell = np.zeros((sizeX,sizeY))
    _diffusion = 1
    _hx = 1
    scale = _diffusion * timestep / (_hx * _hx) 
    for j in range(sizeY):
        for i in range(sizeX):
            if cell[i,j] != 1: # e.g. not fluid cell
                continue
            if (i < sizeX - 1) & (cell[i+1,j] == 1):
                _aDiag[i,j] += scale
                _aDiag[i+1,j] += scale
                _aPlusX[i,j] = -scale
            if (j < sizeY - 1) & (cell[i,j+1] == 1):
                _aDiag[i,j] += scale
                _aDiag[i,j+1] += scale
                _aPlusY[i,j] = -scale
                
def buildPreconditioner():
    tau = 0.97
    sigma = 0.25
    cell = np.zeros((sizeX,sizeY))
    _precon = np.zeros((sizeX,sizeY))
    _aDiag = np.zeros((sizeX,sizeY))
    _aPlusX = np.zeros((sizeX,sizeY))
    _aPlusY = np.zeros((sizeX,sizeY))
    for j in range(sizeY):
        for i in range(sizeX):
            if cell[i,j] != 1: # e.g. not fluid cell
                continue
        e = _aDiag[i,j]
        if (i > 0) & (cell[i-1,j] == 1):
            px = _aPlusX[i-1,j] * _precon[i-1,j]
            py = _aPlusY[i-1,j] * _precon[i-1,j]
            e = e - (px*px + tau*px*py)
        if (j > 0) & (cell[i,j-1] == 1):
            px = _aPlusX[i,j-1] * _precon[i,j-1]
            py = _aPlusY[i,j-1] * _precon[i,j-1]
            e = e - (py*py + tau*px*py)
        if e < sigma * _aDiag[i,j]:
            e = _aDiag[i,j]
        _precon[i,j] = 1 / np.sqrt(e)
        
def fromParticleToGrid():
    posX = np.zeros((particleCount))
    posY = np.zeros((particleCount))
    weight = np.zeros((sizeX,sizeY))
    _src = np.zeros((sizeX,sizeY))
    propertity = np.zeros((particleCount)) # 粒子的某种属性
    _ox = 0
    _oy = 0
    for i in range(particleCount):
        x = posX[i] - _ox
        y = posY[i] - _ox
        x = max(0.5,min(sizeX - 1.5,x))
        y = max(0.5,min(sizeY - 1.5,x))
        
        x0 = int(x)
        y0 = int(y)
        k = (1 - abs(x0 - x))*(1 - abs(y0 - y))
        weight[x0,y0] += k
        _src[x0,y0] += k * propertity[i]
        
def applyPreconditioner():
    _precon = np.zeros((sizeX,sizeY))
    dst = np.zeros((sizeX,sizeY))
    a = np.zeros((sizeX,sizeY))
    _aPlusX = np.zeros((sizeX,sizeY))
    _aPlusY = np.zeros((sizeX,sizeY))
    cell = np.zeros((sizeX,sizeY))
    
    for j in range(sizeY):
        for i in range(sizeX):
            if cell[i,j] != 1: # e.g. not fluid cell
                continue
            t = a[i,j]
            if (i > 0) & (cell[i-1,j] == 1):
                t -= _aPlusX[i-1,j]*_precon[i-1,j]*dst[i-1,j]
            if (j > 0) & (cell[i,j-1] == 1):
                t -= _aPlusY[i,j-1]*_precon[i,j-1]*dst[i,j-1]
            dst[i,j] = t * _precon[i,j]
            
    for j in range(sizeY):
        for i in range(sizeX):
            if cell[i,j] != 1: # e.g. not fluid cell
                continue
            t = a[i,j]
            if (i < sizeX - 1) & (cell[i+1,j] == 1):
                t -= _aPlusX[i+1,j]*_precon[i+1,j]*dst[i+1,j]
            if (j < sizeY - 1) & (cell[i,j+1] == 1):
                t -= _aPlusY[i,j+1]*_precon[i,j+1]*dst[i,j+1]
            dst[i,j] = t * _precon[i,j]
            
def project():
    a = np.zeros((sizeX,sizeY))
    b = np.zeros((sizeX,sizeY))
    cell = np.zeros((sizeX,sizeY))
    res = 0
    for j in range(sizeY):
        for i in range(sizeX):
            if cell[i,j] == 1: # fluid cell
               res += a[i,j] * b[i,j]            
    return res

_z = np.zeros((sizeX,sizeY)) # Auxiliary vector
_s = np.zeros((sizeX,sizeY)) # Search vector
_precon = np.zeros((sizeX,sizeY)) # Preconditioner
_aDiag = np.zeros((sizeX,sizeY)) # Matrix diagonal
_aPlusX = np.zeros((sizeX,sizeY)) # Matrix off_diagonals
_aPlusY = np.zeros((sizeX,sizeY)) 

def matrixVectorProduct():
    dst = np.zeros((sizeX,sizeY)) # _z
    b = np.zeros((sizeX,sizeY)) # _s
    _aDiag = np.zeros((sizeX,sizeY))
    _aPlusX = np.zeros((sizeX,sizeY))
    _aPlusY = np.zeros((sizeX,sizeY))
    for j in range(sizeY):
        for i in range(sizeX):
            t = _aDiag[i,j] * b[i,j]
            if i > 0:
                t += _aPlusX[i-1,j] * b[i-1,j]
            if j > 0:
                t += _aPlusY[i,j-1] * b[i,j-1]
            if i < sizeX - 1:
                t -= _aPlusX[i,j] * b[i+1,j]
            if i < sizeY - 1:
                t -= _aPlusY[i,j] * b[i,j+1]
            dst[i,j] = t
            
    