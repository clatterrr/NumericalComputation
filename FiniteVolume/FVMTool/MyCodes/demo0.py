import numpy as np
'''
FVD Tool ConvectionTVDexample

'''
nmax = 500
phi = np.zeros((nmax))
for i in range(nmax):
    if (i >= 19) & (i <= 120):
        phi[i] = 1
    if (i >= 179) & (i <= 400):
        phi[i] = np.sin((i * 0.02 + 0.01) * np.pi)
dt = 0.0005
dx = 0.002
tmax = 1000
dphi = np.zeros((nmax))
rp = np.zeros((nmax))
psiplus = np.zeros((nmax))
psiminus = np.zeros((nmax))
flux = np.zeros((nmax))
eps = 1e-10
u = np.zeros((nmax + 1))
u[:] = 0.3
M = np.zeros((nmax+2,nmax+2))
for i in range(2,nmax+1):
        M[i,i] = 2000 + 150
        M[i,i-1] = - 150
        
M[0,0] = 1
M[0,1] = 1
M[1,0] = -75
M[1,1] = 2075

M[501,0] = 1
M[0,1] = 1
M[501,1] = -1
M[0,500] = -1
M[501,500] = -1
M[0,501] = -1
M[501,501] = 1
rhst = np.zeros((nmax + 2))
rhs = np.zeros((nmax + 2))
phit = np.zeros((nmax,tmax+1))
phit[:,0] = phi[:]
for t in range(tmax):
    
    rhst[1:nmax+1] = phi[:] / dt
    
    for i in range(nmax-1):
        dphi[i] = (phi[i+1] - phi[i])/dx
        
    for i in range(nmax-2):
        phiout = (abs(dphi[i+1]) >= eps)*dphi[i+1] + eps*(dphi[i+1] == 0) + eps*(abs(dphi[i+1]) < eps)*np.sign(dphi[i+1])
        rp[i] = dphi[i] / phiout
        
        flux[i] = max(0,max(min(2*rp[i],1),min(rp[i],2)))
        
        psiplus[i+1] = 0.5 * flux[i] * (phi[i+2] - phi[i+1])
        
        phiout = (abs(dphi[i]) >= eps)*dphi[i] + eps*(dphi[i] == 0) + eps*(abs(dphi[i]) < eps)*np.sign(dphi[i])
        rp[i] = dphi[i+1] / phiout
        
        flux[i] = max(0,max(min(2*rp[i],1),min(rp[i],2)))
        
        psiminus[i] = 0.5 * flux[i] * (phi[i] - phi[i+1])
        

    for i in range(1,nmax):
        rhs[i+1] = - 1 / dx * u[i]*(psiplus[i] - psiplus[i-1])
        
    rhs = rhs + rhst
    Minv = np.linalg.inv(M)
    for i in range(1,nmax+1):
        summ = 0
        for j in range(nmax+2):
            summ += Minv[i,j] * rhs[j]
        phi[i-1] = summ
    phi[0] = 0
    phit[:,t+1] = phi[:]
    if t % 100 == 0:
        print(t)


from matplotlib import pyplot as plt
from matplotlib import animation
fig = plt.figure()
ax = plt.axes(xlim=(0,500), ylim=(-2, 2))
line, = ax.plot([], [], lw=2)
def init():
    line.set_data([], [])
    return line,
def animate(i):
    x = np.linspace(0, 500, 500)
    y = phit[:,i]
    line.set_data(x, y)
    return line,
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=tmax, interval=20, blit=True)
anim.save('basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

plt.show()