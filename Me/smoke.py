import numpy as np
Nx = 32
Ny = 32
Nz = 32

density = np.zeros((Nx,Ny,Nz))
densityTemp = np.zeros((Nx,Ny,Nz))
densityTemp2 = np.zeros((Nx,Ny,Nz))
Velocityx = np.zeros((Nx,Ny,Nz))
Velocityy = np.zeros((Nx,Ny,Nz))
Velocityz = np.zeros((Nx,Ny,Nz))
Velocityxt = np.zeros((Nx,Ny,Nz))
Velocityyt = np.zeros((Nx,Ny,Nz))
Velocityzt = np.zeros((Nx,Ny,Nz))
Velocityxt2 = np.zeros((Nx,Ny,Nz))
Velocityyt2 = np.zeros((Nx,Ny,Nz))
Velocityzt2 = np.zeros((Nx,Ny,Nz))
Pressure = np.zeros((Nx,Ny,Nz))
Divergence = np.zeros((Nx,Ny,Nz))



emiiterSize = 4
density[Nx//2-emiiterSize//2:Nx//2+emiiterSize//2,Ny//2-emiiterSize//2:Ny//2+emiiterSize//2,emiiterSize//2:emiiterSize] = 1
Velocityy[:,:,:] = 3
dt = 0.1

def clampSize(v):
    if v < 0:
        return int(0)
    elif v > Nx - 2:
        return int(Nx - 2)
    else:
        return int(v)

def advectionDensity():
    densityTemp[:,:,:] = 0
    for z in range(1,Nz-1):
        for y in range(1,Ny-1):
            for x in range(1,Nx-1):
                oldx = x - dt * Velocityx[x,y,z]
                oldy = y - dt * Velocityy[x,y,z]
                oldz = z - dt * Velocityz[x,y,z]
                gx = clampSize(oldx)
                gy = clampSize(oldy)
                gz = clampSize(oldz)
                fx = oldx - gx
                fy = oldy - gy
                fz = oldz - gz
                
                if (x == 15) & (y == 15) & (z == 3):
                    test = 1
                
                v000 = density[gx,gy,gz]
                v100 = density[gx+1,gy,gz]
                v010 = density[gx,gy+1,gz]
                v001 = density[gx,gy,gz+1]
                v110 = density[gx+1,gy+1,gz]
                v011 = density[gx,gy+1,gz+1]
                v101 = density[gx+1,gy,gz+1]
                v111 = density[gx+1,gy+1,gz+1]
                
                densityTemp[x,y,z] = (v000 * (1 - fx) * (1 - fy) * (1 - fz) + 
                                      v100 * fx * (1 - fy) * (1 - fz) +
                                      v010 * (1 - fx) * fy * (1 - fz) +
                                      v001 * (1 - fx) * (1 - fy) * fz +
                                      v110 * fx * fy * (1 - fz) +
                                      v011 * (1 - fx) * fy * fz +
                                      v101 * fx * (1 - fy) * fz +
                                      v111 * fx * fy * fz)
                
def advectionVelocity():
    for z in range(1,Nz-1):
        for y in range(1,Ny-1):
            for x in range(1,Nx-1):
                oldx = x - dt * Velocityx[x,y,z]
                oldy = y - dt * Velocityy[x,y,z]
                oldz = z - dt * Velocityz[x,y,z]
                gx = clampSize(oldx)
                gy = clampSize(oldy)
                gz = clampSize(oldz)
                fx = oldx - gx
                fy = oldy - gy
                fz = oldz - gz
                
                v000 = Velocityx[gx,gy,gz]
                v100 = Velocityx[gx+1,gy,gz]
                v010 = Velocityx[gx,gy+1,gz]
                v001 = Velocityx[gx,gy,gz+1]
                v110 = Velocityx[gx+1,gy+1,gz]
                v011 = Velocityx[gx,gy+1,gz+1]
                v101 = Velocityx[gx+1,gy,gz+1]
                v111 = Velocityx[gx+1,gy+1,gz+1]
                
                Velocityxt[x,y,z] = (v000 * (1 - fx) * (1 - fy) * (1 - fz) + 
                                      v100 * fx * (1 - fy) * (1 - fz) +
                                      v010 * (1 - fx) * fy * (1 - fz) +
                                      v001 * (1 - fx) * (1 - fy) * fz +
                                      v110 * fx * fy * (1 - fz) +
                                      v011 * (1 - fx) * fy * fz +
                                      v101 * fx * (1 - fy) * fz +
                                      v111 * fx * fy * fz)   
                
                v000 = Velocityy[gx,gy,gz]
                v100 = Velocityy[gx+1,gy,gz]
                v010 = Velocityy[gx,gy+1,gz]
                v001 = Velocityy[gx,gy,gz+1]
                v110 = Velocityy[gx+1,gy+1,gz]
                v011 = Velocityy[gx,gy+1,gz+1]
                v101 = Velocityy[gx+1,gy,gz+1]
                v111 = Velocityy[gx+1,gy+1,gz+1]
                
                Velocityyt[x,y,z] = (v000 * (1 - fx) * (1 - fy) * (1 - fz) + 
                                      v100 * fx * (1 - fy) * (1 - fz) +
                                      v010 * (1 - fx) * fy * (1 - fz) +
                                      v001 * (1 - fx) * (1 - fy) * fz +
                                      v110 * fx * fy * (1 - fz) +
                                      v011 * (1 - fx) * fy * fz +
                                      v101 * fx * (1 - fy) * fz +
                                      v111 * fx * fy * fz)  
                
                v000 = Velocityz[gx,gy,gz]
                v100 = Velocityz[gx+1,gy,gz]
                v010 = Velocityz[gx,gy+1,gz]
                v001 = Velocityz[gx,gy,gz+1]
                v110 = Velocityz[gx+1,gy+1,gz]
                v011 = Velocityz[gx,gy+1,gz+1]
                v101 = Velocityz[gx+1,gy,gz+1]
                v111 = Velocityz[gx+1,gy+1,gz+1]
                
                Velocityzt[x,y,z] = (v000 * (1 - fx) * (1 - fy) * (1 - fz) + 
                                      v100 * fx * (1 - fy) * (1 - fz) +
                                      v010 * (1 - fx) * fy * (1 - fz) +
                                      v001 * (1 - fx) * (1 - fy) * fz +
                                      v110 * fx * fy * (1 - fz) +
                                      v011 * (1 - fx) * fy * fz +
                                      v101 * fx * (1 - fy) * fz +
                                      v111 * fx * fy * fz)  
                
                
                
                
time = 0
timeFinal = 1
while(time < timeFinal):
    densityTemp2 = density.copy()
    Velocityxt2 = Velocityx.copy()
    Velocityyt2 = Velocityy.copy()
    Velocityzt2 = Velocityz.copy()
    advectionDensity()
    advectionVelocity()
    density = densityTemp.copy()
    Velocityx = Velocityxt.copy()
    Velocityy = Velocityyt.copy()
    Velocityz = Velocityzt.copy()
    # advectionDensity()
    # advectionVelocity()
    # density = (densityTemp + densityTemp2) / 2
    # Velocityx = (Velocityxt + Velocityxt2) / 2
    # Velocityy = (Velocityyt + Velocityyt2) / 2
    # Velocityz = (Velocityzt + Velocityzt2) / 2
    time = 5
    
    

