import numpy as np
import scipy.io as scio
u = scio.loadmat('linwave_u.mat')['u']
N = 128
timeFinal = 4
CFL = 0.2
xmin = -1
xmax = 1
m = 1
h = (xmax - xmin) / N
x = np.zeros((2,N))
for i in range(0,N):
    x[0,i] = xmin + i * h
    x[1,i] = x[0,i] + h
    
def LegendreGL(m):
    localx = np.zeros((m+1))
    localw = np.zeros((m+1))
    if m == 1:
        localx = np.array([-1,1])
        localw = np.array([1,1])
    elif m == 2:
        localx = np.array([-1,0,1])
        localw = np.array([1/3,4/3,1/3])
    else:
        J = np.zeros((m-1))
        h1 = np.linspace(0,m-2,m-1)*2 + 2
    return localx,localw


def LegendreGQ(m):
    localx = np.zeros((m+1))
    localw = np.zeros((m+1))
    if m == 0:
        localx[0] = 0
        localw[0] = 2
    else:
        J = np.zeros((m+1,m+1))
        h1 = np.linspace(0,m,m+1)*2
        for i0 in range(m):
            J[i0,i0+1] = J[i0+1,i0] = 2/(h1[i0] + 2)*np.sqrt(
                (i0+1)**4/(h1[i0] + 1)/(h1[i0] + 3))
        localx,localV = np.linalg.eig(J)
        localw = 2 * (localV[0,:])**2
    return localx,localw
        
    return localx,localw

def LegendreP(localx,localm):
    PL = np.zeros((localm+1,len(localx)))
    PL[0,:] = np.sqrt(1 / 2)
    if localm == 0:
        return PL[0,:]
    PL[1,:] = np.sqrt(3 / 2)*localx
    if localm == 1:
        return PL[1,:]
    aold = np.sqrt(1 / 3)
    for i in range(localm - 1):
        anew = 2/(2*i+2)*np.sqrt((i+1)**4/(2*i+1)/(2*i+3));
        PL[i+2,:] = 1/anew*(-aold*PL[i,:] + localx*PL[i+1,:]);
        aold = anew
    return PL[localm,:]
    
def GradLegendreP(r0,m0):
    dP = np.zeros((len(r0)))
    if m0 > 0:
        Ph = -m0*r0*LegendreP(r0,m0) + m0*np.sqrt((2*m0+1)/(2*m0-1))*LegendreP(r0,m0-1)
        dPe = r0**(m0+1)*m0*(m0+1)/2*np.sqrt((2*m0+1)/2)
        endP = (abs(abs(r)-1)>10*1e-10)
        rh = r0 * endP
        dP = ~endP * dPe + endP * Ph / (1 - rh**2)
    return dP

def VandermondeDG(m,r):
    V = np.zeros((len(r),m+1))
    for j in range(m+1):
        V[:,j] = LegendreP(r,j)
    return V

def GradVandermondeDG(m,r):
    Vr = np.zeros((len(r),m+1))
    for j in range(m+1):
        Vr[:,j] = GradLegendreP(r,j)
    return Vr

r,waste = LegendreGL(m)
V = VandermondeDG(m,r)
Vr = GradVandermondeDG(m,r)
D = np.dot(Vr,np.linalg.inv(V))
Ma = np.linalg.inv(np.dot(V,np.transpose(V)))
S = np.dot(Ma,D)
Vinv = np.linalg.inv(V)
Q = np.zeros((m+1,m+1))
Pmat = np.zeros((m+1,m+1))
Xm = np.zeros((m+1,m+1))
Xp = np.zeros((m+1,m+1))
x,w = LegendreGQ(m)
x = -x
for i in range(m+1):
    Pmat[i,:] = LegendreP(x,i)
    Xm[i,:] = LegendreP(x-2,i)
    Xp[i,:] = LegendreP(x+2,i)
for l in range(m):
    lamb = np.zeros((m+1,m+1))
    for i in range(m+1):
        lamb[i,i] = w[i]
    # Set up operator to recover derivaties
    A = np.zeros((m+1-l,m+1-l))
    A[0,0] = 1 / np.sqrt((2 * l + 3) * (2 * l + 1))
    A[m-l,m-l] = 1 / np.sqrt(2 * (m + 2) + 1) / np.sqrt(2 * (m+2) - 1)
    for i in range(1,m-l):
        Ah = 1 / np.sqrt(2*(l + i) + 1) / np.sqrt(2*(l + i) - 1)
        A[i,i] = Ah
        A[i+1,i-1] = - Ah
    # Recover derivatives at quadrature points
    Ph1 = np.dot(np.linalg.inv(A),Pmat[l:m+1,:])    
    Pmat[0:l+1,:] = 0
    Pmat[l+1:m+1,:] = Ph1[0:m-l,:]
    
    # Compute smoothness operator for order l and update
    Qtemp = np.dot(Pmat,lamb)
    Qh = np.dot(Qtemp,np.transpose(Pmat))
    Q = Q + 2**(2*l+1)*Qh
    
# Initialize operator for smoothness indicator in nodal space
Q = np.dot(np.dot(np.transpose(Vinv),Q),Vinv)
Xp = np.dot(np.transpose(Vinv),Xp)
Xm = np.dot(np.transpose(Vinv),Xm)

# Initialize extraction vector
VtoE = np.zeros((2,N))
for j in range(N):
    VtoE[0,j] = j * (m + 1)
    VtoE[1,j] = (j + 1) * (m + 1) - 1
    
# Compute smallest spatial scale timestep
rLGLmin = abs(r[0] - r[1]) 
time = 0
tstep = 0

# Set timestep
k = CFL*rLGLmin*h;
    
time = 0
timeFinal = 4
while(time < timeFinal):
    time = 5
    
    