import numpy as np
import matplotlib.pyplot as plt
N = 12
dt = 5e-5
dx = 1 / N
rho = 4e1
num_face = 2 * N ** 2   # number of faces
num_vert = (N + 1) ** 2 # number of vertices
E, nu = 4e4, 0.2  # Young's modulus and Poisson's ratio
mu, lam = E / 2 / (1 + nu), E * nu / (1 + nu) / (1 - 2 * nu) # Lame parameters
ball_pos, ball_radius = np.array([0.5, 0.0]), 0.31
damping = 14.5

pos = np.zeros((num_vert,2))
posold = np.zeros((num_vert,2))
vel = np.zeros((num_vert,2))
element = np.zeros((num_face,3),dtype = int)
Bmat = np.zeros((num_face,2,2))
Binv = np.zeros((num_face,2,2))
Fmat = np.zeros((num_face,2,2))
V = np.zeros((num_face))
phi = np.zeros((num_face)) # 每个面的势能

gravity = np.array([0,-1])
attractor_pos = np.array([0,0])
attractor_strength = 0

def init_mesh():
    for j in range(N):
        for i in range(N):
            idx = (j * N + i) * 2
            a = j * (N + 1) + j
            b = a + 1
            c = a + N + 2
            d = a + N + 1
            '''
            d ---- c
            
            
            a ---- b
            '''
            element[idx,:] = [a,b,c]
            element[idx+1,:] = [c,d,a]
            
def init_pos():
    for j in range(N+1):
        for i in range(N+1):
            idx = j * (N + 1) + i
            pos[idx,0] = i / N / 4 +  0.45
            pos[idx,1] = j / N / 4 + 0.45
            vel[idx,0] = 0
            vel[idx,1] = 0
    for i in range(num_face):
        ia = element[i,0]
        ib = element[i,1]
        ic = element[i,2]
        ax = pos[ia,0]
        ay = pos[ia,1]
        bx = pos[ib,0]
        by = pos[ib,1]
        cx = pos[ic,0]
        cy = pos[ic,1]
        B_i_inv = np.array([[ay-cy,by-cy],[ax-cx,bx-cx]])
        B_i_inv = np.array([[ax-cx,bx-cx],[ay-cy,by-cy]])
        Binv[i,:,:] = B_i_inv 
        Bmat[i] = np.linalg.inv(B_i_inv)
        
def update():
    U = 0 # 总势能
    for i in range(num_face):
        ia = element[i,0]
        ib = element[i,1]
        ic = element[i,2]
        ax = pos[ia,0]
        ay = pos[ia,1]
        bx = pos[ib,0]
        by = pos[ib,1]
        cx = pos[ic,0]
        cy = pos[ic,1]
        V[i] = (ay-cy)*(ax-cx) + (by-cy)*(bx-cx)
        D_i = np.array([[ay-cy,by-cy],[ax-cx,bx-cx]])
        D_i = np.array([[ax-cx,bx-cx],[ay-cy,by-cy]])
        Fmat[i,:,:] = np.dot(D_i,Bmat[i,:,:])
    for i in range(num_face):
        F_i = Fmat[i,:,:]
        log_J_i = np.log(np.linalg.det(F_i))
        phi_i = mu / 2 * (np.trace(np.dot(np.transpose(F_i),F_i))-2)
        phi_i -= mu * log_J_i
        phi_i += lam / 2 * log_J_i ** 2
        phi[i] = phi_i
        U += V[i] * phi_i

def advance():
    for i in range(num_vert):
        acc = 0
        i0 = int(i % (N + 1))
        j0 = int(i // (N + 1))
        if i0 > 0:
            acc += (vel[i-1,:] - vel[i,:]) / (rho * dx ** 2)
        if i0 < N:
            acc += (vel[i+1,:] - vel[i,:]) / (rho * dx ** 2)
        if j0 > 0:
            acc += (vel[i-N-1,:] - vel[i,:]) / (rho * dx ** 2)
        if j0 < N:
            acc += (vel[i+N+1,:] - vel[i,:]) / (rho * dx ** 2)
        acc += (vel[i,:] - posold[i,:]) / (rho * dx ** 2)
        g = gravity * 0.8 
        vel[i,:] += dt * (acc + g * 40)
        vel[i,:] *= np.exp(-dt * damping)
    for i in range(num_vert):
        disp = pos[i,:] - ball_pos
        disp2 = np.sqrt(disp[0]**2 + disp[1]**2)
        if disp2 <= ball_radius:
            Nov = vel[i,0] * disp[0] + vel[i,1] *disp[1]
            if Nov < 0:
                vel[i,:] -= Nov * disp * 10
        if pos[i,0] < 0 or pos[i,1] < 0 or pos[i,0] > 1 or pos[i,1] > 1:
            vel[i,:] = 0
        pos[i,:] += dt * vel[i,:]
    
        
init_mesh()
init_pos()
time = 0
timeFinal = 10000
while(time < timeFinal):
    posold = vel.copy()
    update()
    advance()
    time += 1
        
    if time % 500 != 0:
        continue
    fig, ax = plt.subplots()
    plt.title('f model: T=%f' %time)
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    for i in range(num_vert):
        circle1 = plt.Circle((pos[i,0], pos[i,1]), 0.01)
        ax.add_artist(circle1)
    plt.show()
    