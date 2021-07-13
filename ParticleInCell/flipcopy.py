import numpy as np
import matplotlib.pyplot as plt

particlePos = np.load("pos32.npy")
n_particles = particlePos.shape[0]
particleVel = np.zeros((n_particles,2))
C = np.zeros((n_particles,2,2))
J = np.zeros((n_particles))
n_grid = 16
gridVel = np.zeros((2,n_grid,n_grid))
gridMass = np.zeros((n_grid,n_grid))
dx = 1 / n_grid
dt = 2e-4

p_rho = 1
p_vol = (dx * 0.5)**2
p_mass = p_vol * p_rho
gravity = 9.8
bound = 3
E = 400

for i in range(n_particles):
    J[i] = 1
    particleVel[i,0] = 0
    particleVel[i,1] = -1
    

time = 0
timeFinal = 10
while(time < timeFinal):
    
    for i in range(n_grid):
        for j in range(n_grid):
            gridVel[0,i,j] = gridVel[1,i,j] = 0
            gridMass[i,j] = 0
            
    for p in range(n_particles):
        px = particlePos[p,0] / dx
        py = particlePos[p,1] / dx
        gx = int(px - 0.5)
        gy = int(py - 0.5)
        fx = px - gx
        fy = py - gy
        wx = [0.5 * (1.5 - fx)**2, 0.75 - (fx - 1)**2, 0.5 * (fx - 0.5)**2]
        wy = [0.5 * (1.5 - fy)**2, 0.75 - (fy - 1)**2, 0.5 * (fy - 0.5)**2]
        stress = -dt * 4 * E * p_vol * (J[p] - 1) / dx**2
        affine = np.array([[stress,0],[0,stress]]) + p_mass * C[p,:,:]
        for i in range(3):
            for j in range(3):
                weight = wx[i] * wy[j]
                dpos = [(i - fx)*dx,(j - fy)*dx]
                aff = affine[0,0] * dpos[0] + affine[0,1] * dpos[1]
                gridVel[0,gx+i,gy+j] += weight * (p_mass * particleVel[p,0] + aff)
                aff = affine[1,0] * dpos[0] + affine[1,1] * dpos[1]
                gridVel[1,gx+i,gy+j] += weight * (p_mass * particleVel[p,1] + aff)
                gridMass[gx+i,gy+j] += weight * p_mass
                
    rest = 1
                
    for i in range(n_grid):
        for j in range(n_grid):
            if gridMass[i,j] > 0:
                gridVel[:,i,j] /= gridMass[i,j]
            gridVel[1,i,j] -= dt * gravity
            if i < bound and gridVel[0,i,j] < 0:
                gridVel[0,i,j] = 0
            if i > n_grid - bound and gridVel[1,i,j] > 0:
                gridVel[0,i,j] = 0
            if j < bound and gridVel[1,i,j] < 0:
                gridVel[1,i,j] = 0
            if j > bound and gridVel[1,i,j] > 0:
                gridVel[1,i,j] = 0
                
    rest = 1
    
    for p in range(n_particles):
        px = particlePos[p,0] / dx
        py = particlePos[p,1] / dx
        gx = int(px - 0.5)
        gy = int(py - 0.5)
        fx = px - gx
        fy = py - gy
        wx = [0.5 * (1.5 - fx)**2, 0.75 - (fx - 1)**2, 0.5 * (fx - 0.5)**2]
        wy = [0.5 * (1.5 - fy)**2, 0.75 - (fy - 1)**2, 0.5 * (fy - 0.5)**2]
        new_vel = np.zeros((2))
        new_C = np.zeros((2,2))
        for i in range(3):
            for j in range(3):
                weight = wx[i] * wy[j]
                dpos = np.array([(i - fx)*dx,(j - fy)*dx])
                new_vel += weight * gridVel[:,gx+i,gy+j]
                vx = gridVel[0,gx+i,gy+j]
                vy = gridVel[1,gx+i,gy+j]
                new_C += 4 * weight  /dx / dx * np.array([[vx*dpos[0],vx*dpos[0]],[vy*dpos[0],vy*dpos[1]]])
        particleVel[p] = new_vel
        particlePos += dt * particleVel
        J[p] *= 1 + dt * np.trace(new_C)
        C[p] = new_C
    
    plt.figure(1)
    fig, ax = plt.subplots()          
    for i in range(n_particles):
        circle1 = plt.Circle((particlePos[i,0], particlePos[i,1]), 0.01)
        ax.add_artist(circle1)
    plt.show()
    time += 1



