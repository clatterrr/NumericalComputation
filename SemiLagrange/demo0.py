import numpy as np
'''
https://github.com/iCFD/SemiLagrangian

'''
CFL = 8.5
a = -1
b = 1
Lx = b - a
Nx = 200
dx = Lx / Nx
x = np.zeros((Nx))
h = np.zeros((Nx))

def Init0(x0,x1,x2):
    return np.exp(-x1 * (x0 - x2)**2)

beta = np.log(2)/(36*0.005**2)
z = -0.7
delta = 0.005
for i in range(Nx):
    x[i] = a + i*dx + dx/2
    if (x[i] >= -0.8) & (x[i] <= -0.6):
        h[i] = (Init0(x[i],beta,z-delta) + Init0(x[i],beta,z+delta) + 4*Init0(x[i],beta,z))/6


tmax = 100
ht = np.zeros((Nx,tmax+1))
xo =  np.zeros((Nx))
ho = np.zeros((Nx))
dt = 0.005
u = -2
ht[:,0] = h[:]
for t in range(tmax):
    for i in range(tmax):
        xo[i] = x[i] + u * dt
        if xo[i] < a:
            xo[i] += Lx
        if xo[i] > b:
            xo[i] -= Lx
        
        idx0 = int((xo[i] - a - dx/2)/dx)
        idx1 = int(idx0 + 1)
        if idx1 >= Nx:
            idx1 = 0
        d = (xo[i] - a - dx)/dx - idx0
        
        ho[i] = (1 - d)*h[idx0] + d*h[idx1]
    
    h = ho.copy()
    x = xo.copy()
    
    ht[:,t+1] = h[:]
    








from matplotlib import pyplot as plt
from matplotlib import animation

# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = plt.axes(xlim=(a,b), ylim=(-2, 2))
line, = ax.plot([], [], lw=2)

# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    return line,

# animation function.  This is called sequentially
def animate(i):
    x = np.linspace(-1, 1, 200)
    y = ht[:,i]
    line.set_data(x, y)
    return line,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=tmax, interval=20, blit=True)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
anim.save('basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

plt.show()