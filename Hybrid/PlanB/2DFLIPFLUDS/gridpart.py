import numpy as np
"""
 代码库计划——只抄一些实用代码，不复现
 2DFLIPFLUID
"""
Nx = 32
Ny = 32
dx = 1 / Nx

def getGridVelocity(i,j):
    return 0

def getGridPressure(i,j):
    return 0

def getGridDensity(i,j):
    return 0

divergence = np.zeros((Nx,Ny))
def computeGridDivergence():
    for j in range(Ny):
        for i in range(Nx):
            divergence[i,j] = ((getGridVelocity(i+1,j)[0] - getGridVelocity(i-1,j)[0])/(2*dx) 
                              + (getGridVelocity(i,j+1)[1] - getGridVelocity(i,j-1)[1])/(2*dx))
            

forcex = np.zeros((Nx,Ny))
forcey = np.zeros((Nx,Ny))
velocityx = np.zeros((Nx,Ny))
velocityy = np.zeros((Nx,Ny))
def computeGridPressureForces():
    for j in range(Ny):
        for i in range(Nx):
            forcex[i,j] = (getGridPressure(i+1, j) - getGridPressure(i-1, j)) / (2 * dx)
            forcey[i,j] = (getGridPressure(i, j+1) - getGridPressure(i, j-1)) / (2 * dx)
def computeVelocityBasedPressureForces():
    for j in range(Ny):
        for i in range(Nx):   
            computeGridPressureForces()
            velocityx[i,j] -= forcex[i,j]
            velocityy[i,j] -= forcey[i,j]
    
def bilnearlyInterpolate(x,y):
    posx = x / dx
    posy = y / dx
    i = int(posx)
    j = int(posy)
    fx = posx - i
    fy = posy - j
    density =  (getGridDensity(i, j) * (1-fx)*(1-fy) +  
            getGridDensity(i+1, j) * fx*(1-fy) +  
            getGridDensity(i, j+1) * (1-fx)*fy +  
            getGridDensity(i+1, j+1) * fx*fy)
    return density

density_temp = np.zeros((Nx,Ny))
density = np.zeros((Nx,Ny))
dt = 1
def advectDensitySemiLagrange():
    
    
    for j in range(Ny):
        for i in range(Nx):
            x = i * dx - getGridVelocity(i, j) * dt
            y = j * dx - getGridVelocity(i, j) * dt
            density_temp[i,j] = bilnearlyInterpolate(x, y)
            
def computeVelocity():
    for j in range(Ny):
        for i in range(Nx):
            velocityx[i,j] = forcex[i,j] * density[i,j] * dt
            velocityy[i,j] = forcey[i,j] * density[i,j] * dt
            
def UpdateGrid():
    
    computeVelocity()
    advectDensitySemiLagrange()
    for loop in range(10):
        computeGridDivergence()
        # computePressure()
        computeVelocityBasedPressureForces
        
def GridToParticle(particlePos):
    posx = particlePos[0] / dx
    posy = particlePos[1] / dx
    i = int(posx)
    j = int(posy)
    fx = posx - i
    fy = posy - j
    velo = (getGridVelocity(i, j) * (1-fx)*(1-fy) + 
            getGridVelocity(i+1, j) * fx*(1-fy) + 
            getGridVelocity(i, j+1) * (1-fx)*fy + 
            getGridVelocity(i+1,j+1) * fx*fy)
    return velo
    
    