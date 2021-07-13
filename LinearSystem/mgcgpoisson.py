import numpy as np
# taichi mgpcg
nx = 32
n_levels = 3
r = np.zeros((n_levels,nx,nx))
z = np.zeros((n_levels,nx,nx)) # M ^-1 r
x = np.zeros((nx,nx)) # 未知量
p = np.zeros((nx,nx)) # 共轭梯度，方向
Ap = np.zeros((nx,nx))
alpha = np.zeros((nx,nx))
beta = np.zeros((nx,nx))
sum = np.zeros((nx,nx))

def inDomain(i,j):
    if i < 0 or j < 0 or i >= nx  or j >= nx:
        return 0
    return 1

def ComputeAp(Ap,p):
    size = Ap.shape[0]
    for i in range(size):
        for j in range(size):
            Ap[i,j] = 4 * p[i,j] - p[i-1,j] - p[i+1,j] - p[i,j-1] - p[i,j+1]

def restrict(l):
    size = 0
    for i in range(size):
        for j in range(size):
            res = r[l,i,j] - (4 * z[l,i,j] - z[l,i-1,j] - z[l,i+1,j] - z[l,i,j-1] - z[l,i,j+1])
         
def smooth(l,phase):
    # phase = red black Gauss-Seidel phase
    for i in range(nx):
        for j in range(nx):
            if (i+j+k) & 1 == phase:
                z[l,i,j] = (r[l,i,j] + z[l,i+1,j] + z[l,i,j+1] + z[l,i-1,j] + z[l,i,j-1])/4
                
            
         
def applyPreconditioner():
    for l in range(n_levels):
        

    