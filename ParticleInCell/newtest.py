import numpy as np
import matplotlib.pyplot as plt
nx = 32


gridu = np.zeros((nx+1,nx))
gridv = np.zeros((nx,nx+1))
griduTemp = np.zeros((nx+1,nx))
gridvTemp = np.zeros((nx,nx+1))
gridUAcc = np.zeros((nx+1,nx))
gridVAcc = np.zeros((nx,nx+1))
gridULast = np.zeros((nx+1,nx))
gridVLast = np.zeros((nx,nx+1))

divergence = np.zeros((nx,nx))
pressure = np.zeros((nx,nx))
pressureTemp = np.zeros((nx+2,nx+2))
n_sqrt = 32
n_particle = n_sqrt * n_sqrt
particlePos = np.zeros((n_particle,2))
particleVel = np.zeros((n_particle,2))
particleFraction = np.zeros((n_particle))

marker = np.zeros((nx,nx)) # 0 for solid ,1 for fluid,2 for air

valid = np.zeros((nx+1,nx+1))
validTemp = np.zeros((nx+1,nx+1))

for i in range(nx):
    for j in range(nx):
        
        if i < nx//2 and j < nx//2:
            marker[i,j] = 1
        else:
            marker[i,j] = 2
            
        if i == 0 or j == 0 or i == nx-1 or j == nx-1:
            marker[i,j] = 0
            
for k in range(n_particle):
    x = int(k % n_sqrt) / n_sqrt * (nx//2 - 1)
    y = int(k / n_sqrt) / n_sqrt * (nx//2 - 1)
    particlePos[k,0] = x + 5 / 4
    particlePos[k,1] = y + 5 / 4
    
dt = 1 / 40.0
gravity = - 1

def applyGravity():
    for i in range(nx):
        for j in range(nx):
            if marker[i,j] == 1:
                gridv[i,j] += dt * gravity * (nx - j)

def boundary():
    for i in range(nx):
        for j in range(nx):
            if marker[i,j] == 0:
                gridu[i,j] = 0
                gridu[i+1,j] = 0
                gridv[i,j] = 0
                gridv[i,j+1] = 0
                
    # for i in range(1,nx):
    #     for j in range(nx):
    #         if marker[i-1,j] == 0 and marker[i,j] == 1:
    #             gridu[i,j] = -gridu[i+1,j]
    #         if marker[i,j] == 1 and marker[i+1,j] == 1:
    #             gridu[i+1,j] = -gridu[i,j]
    # for i in range(nx):
    #     for j in range(1,nx):
    #         if marker[i,j-1] == 0 and marker[i,j] == 1:
    #             gridv[i,j] = -gridv[i,j+1]
    #         if marker[i,j] == 1 and marker[i,j+1] == 1:
    #             gridv[i,j+1] = -gridv[i,j]
    
def inDomain(i,j):
    if i >= 0 and i < nx and j >= 0 and j < nx:
        return 1
    else:
        return 0
    
# 这玩意从来没想过
def diffusion():
    
    diffTime = 10
    
    # 这个时候valid数组最上面一行是没有意义的，因为要配合gridu
    # 该配合你的演出我却视而不见？
    for i in range(nx):
        for j in range(nx):
            if (inDomain(i, j) and marker[i,j] == 1) or (inDomain(i-1, j) and marker[i-1,j] == 1):
                valid[i,j] = 1
            else:
                valid[i,j] = 0
                
    for k in range(diffTime):
        griduTemp = gridu.copy()
        validTemp = valid.copy()
        
        for i in range(nx+1):
            for j in range(nx):
                if validTemp[i,j] == 0:
                    sum = 0
                    count = 0
                    for i0 in range(-1,2):
                        for j0 in range(-1,2):
                            n = i + i0
                            m = j + j0
                            if inDomain(n, m) and validTemp[n,m] == 1:
                                sum += griduTemp[i,j]
                                count += 1
                    if count > 0:
                        gridu[i,j] = sum / count
                        valid[i,j] = 1
                                
                
    
    # 这个时候valid数组最右面一列是没有意义的，因为要配合gridv
    for i in range(nx):
        for j in range(nx):
            if (inDomain(i, j) and marker[i,j] == 1) or ( inDomain(i, j-1) and marker[i,j-1] == 1):
                valid[i,j] = 1
            else:
                valid[i,j] = 0
                
    for k in range(diffTime):
        gridvTemp = gridv.copy()
        validTemp = valid.copy()    
    
        for i in range(nx):
            for j in range(nx+1):
                if validTemp[i,j] == 0:
                    sum = 0
                    count = 0
                    for i0 in range(-1,2):
                        for j0 in range(-1,2):
                            n = i + i0
                            m = j + j0
                            if inDomain(n, m) and validTemp[n,m] == 1:
                                sum += gridvTemp[i,j]
                                count += 1
                    if count > 0:
                        gridv[i,j] = sum / count
                        valid[i,j] = 1
        
def calcDivergence():
    for i in range(nx):
        for j in range(nx):
            if marker[i,j] == 1:
                divergence[i,j] = (gridu[i+1,j] - gridu[i,j]) + (gridv[i,j+1] - gridv[i,j])
            else:
                divergence[i,j] = 0
    
# 耗时间的垃圾解法，下次换preconditioned conjugate gradient，也就是PCG
def solvePressure():
    divflat = np.zeros((nx * nx))
    pflat = np.zeros((nx * nx))
    A = np.zeros((nx*nx,nx*nx))
    for k in range(nx*nx):
        A[k,k] = 4
        i = int(k % nx)
        j = int(k / nx)
        if inDomain(i-1, j):
            A[k,k-1] = -1
        if inDomain(i+1, j):
            A[k,k+1] = -1
        if inDomain(i, j-1):
            A[k,k-nx] = -1
        if inDomain(i, j+1):
            A[k,k+nx] = -1
    for i in range(nx):
        for j in range(nx):
            idx = j * nx + i
            divflat[idx] = divergence[i,j]
    pflat = np.dot(np.linalg.inv(A),-divflat)
    for i in range(nx):
        for j in range(nx):
            idx = j * nx + i
            pressure[i,j] = pflat[idx]
            
    pressureTemp[:,:] = 0
    for k in range(1000):
        for i in range(1,nx+1):
            for j in range(1,nx+1):
                if marker[i-1,j-1] == 2:
                    continue
                pressureTemp[i,j] = (pressureTemp[i+1,j] + pressureTemp[i-1,j]
                                  + pressureTemp[i,j+1] + pressureTemp[i,j-1] - divergence[i-1,j-1])/4
        pressureTemp[0,:] = pressureTemp[1,:]
        pressureTemp[:,0] = pressureTemp[:,1]
        pressureTemp[nx+1,:] = pressureTemp[nx,:]
        pressureTemp[:,nx+1] = pressureTemp[:,nx]
            
    
def applyPressure():
    for i in range(nx):
        for j in range(nx):
            if (inDomain(i, j) and marker[i,j] == 1) or (inDomain(i-1, j) and marker[i-1,j] == 1):
                # gridu[i,j] += (pressure[i-1,j] - pressure[i,j])
                gridu[i,j] += (pressureTemp[i,j+1] - pressureTemp[i+1,j+1])
            if (inDomain(i, j) and marker[i,j] == 1) or (inDomain(i, j-1) and marker[i,j-1] == 1):
                # gridv[i,j] += (pressure[i,j-1] - pressure[i,j]) 
                gridv[i,j] += (pressureTemp[i+1,j] - pressureTemp[i+1,j+1])
    


def interpolate(pos,u,v):
    x = pos[0]
    y = pos[1]
    gx = int(x)
    gy = int(y)
    fx = x - gx
    fy = y - gy
    
    term1 = u[gx,gy] * (1 - fx) * (1 - fy) + u[gx+1,gy] * fx * (1 - fy)
    term2 = u[gx,gy+1] * (1 - fx) * fy + u[gx+1,gy+1] * fx * fy
    
    term3 = v[gx,gy] * (1 - fx) * (1 - fy) + v[gx+1,gy] * fx * (1 - fy)
    term4 = v[gx,gy+1] * (1 - fx) * fy + v[gx+1,gy+1] * fx * fy
    
    return np.array([term1+term2,term3+term4])
    

def gridToParticle():
    flip = 0.95 # Amazaing
    for i in range(nx+1):
        for j in range(nx):
            gridUAcc[i,j] = gridu[i,j] - gridULast[i,j]
    for i in range(nx):
        for j in range(nx+1):
            gridVAcc[i,j] = gridv[i,j] - gridVLast[i,j]
    for k in range(n_particle):
        vel_now = interpolate(particlePos[k,:],gridu,gridv)
        acc_now = interpolate(particlePos[k,:],gridUAcc,gridVAcc)
        particleVel[k,:] = flip * vel_now[:] + (1 - flip) * (particleVel[k,:] + acc_now[:])
    
# 很神奇，居然不用particleVelocity?
def particleAdvect():
    for k in range(n_particle):
        mid = particlePos[k,:] + dt * interpolate(particlePos[k,:],gridu,gridv) * 0.5
        particlePos[k,:] = mid + dt * interpolate(mid,gridu,gridv)
        # particlePos[k,:] += dt * particleVel[k,:]
    
def updateMarker():
    
    # 先将所有水和所有空气统一设置成空气
    for i in range(nx):
        for j in range(nx):
            if marker[i,j] == 1:
                marker[i,j] = 2
    
    # 如果网格里至少有一个粒子，则设置成水
    for k in range(n_particle):
        i = int(particlePos[k,0])
        j = int(particlePos[k,1])
        if marker[i,j] == 2:
            marker[i,j] = 1

def particleToGrid():
    gridu[:,:] = 0
    gridv[:,:] = 0
    griduw = np.zeros((nx+1,nx))
    gridvw = np.zeros((nx,nx+1))
    for k in range(n_particle):
        x = particlePos[k,0]
        y = particlePos[k,1]
        gx = int(x)
        gy = int(y)
        fx = x - gx
        fy = y - gy
        
        gridu[gx,gy] += (1 - fx) * (1 - fy) * particleVel[k,0]
        gridu[gx+1,gy] += fx * (1 - fy) * particleVel[k,0]
        gridu[gx,gy+1] += (1 - fx) * fy * particleVel[k,0]
        gridu[gx+1,gy+1] += fx * fy * particleVel[k,0]
        
        griduw[gx,gy] += (1 - fx) * (1 - fy)
        griduw[gx+1,gy] += fx * (1 - fy)
        griduw[gx,gy+1] += (1 - fx) * fy
        griduw[gx+1,gy+1] += fx * fy
        
        gridv[gx,gy] += (1 - fx) * (1 - fy) * particleVel[k,1]
        gridv[gx,gy+1] += (1 - fx) * fy * particleVel[k,1]
        gridv[gx+1,gy] += fx * (1 - fy) * particleVel[k,1]
        gridv[gx+1,gy+1] += fx * fy * particleVel[k,1]
        
        gridvw[gx,gy] += (1 - fx) * (1 - fy)
        gridvw[gx,gy+1] += (1 - fx) * fy
        gridvw[gx+1,gy] += fx * (1 - fy)
        gridvw[gx+1,gy+1] += fx * fy
        
    for i in range(nx+1):
        for j in range(nx):
            if griduw[i,j] == 0:
                gridu[i,j] = 0
            else:
                gridu[i,j] /= griduw[i,j]
                
    for i in range(nx):
        for j in range(nx+1):
            if gridvw[i,j] == 0:
                gridv[i,j] = 0
            else:
                gridv[i,j] /= gridvw[i,j]

def saveVelocity():
    gridULast = gridu.copy()
    gridVLast = gridv.copy()

time = 0
timeFinal = 1000
while(time < timeFinal):
    time += 1
    
    applyGravity()
    boundary()
    
    diffusion()
    boundary()
    
    calcDivergence()
    solvePressure()
    applyPressure()
    
    diffusion()
    boundary()
    
    
    gridToParticle()
    particleAdvect()
    
    
    
    updateMarker()
    particleToGrid()
    
    saveVelocity()
    
    if (time % 5) != 0:
        continue
    
    fig, ax = plt.subplots()
    plt.title('f model: T=%f' %time)
    ax.set_xlim(0.8,nx-0.2)
    ax.set_ylim(0.8,nx-0.2)
    for i in range(n_particle):
        circle1 = plt.Circle((particlePos[i,0], particlePos[i,1]), 0.1)
        ax.add_artist(circle1)
    plt.show()
    
