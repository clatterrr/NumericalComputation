import numpy as np

nmax = 20

frac = np.zeros((nmax,nmax))
momx = np.zeros((nmax,nmax))
momy = np.zeros((nmax,nmax))
u = np.zeros((nmax,nmax))
v = np.zeros((nmax,nmax))
p = np.zeros((nmax,nmax))
volume = np.zeros((nmax,nmax))
normalx = np.zeros((nmax,nmax))
normaldir = np.zeros((nmax,nmax))
frac[0:10,0:10] = 1
volume[0:10,0:10] = 1
normalx[0:10,0:10] = 1


for j in range(nmax//2):
    p[0:10,j] = (nmax//2 - j) / 10
for i in range(1,nmax-1):
    for j in range(0,nmax):
        u[i,j] = (p[i-1,j] - p[i+1,j])/2
for i in range(0,nmax):
    for j in range(1,nmax-1):
        v[i,j] = (p[i,j-1] - p[i,j+1])/2 - 0.1


u[0,:] = u[1,:]
u[nmax-1,:] = u[nmax-2,:]
v[:,0] = v[:,1]
v[:,nmax-1] = v[:,nmax-2]
uface = np.zeros((nmax-1,nmax))
vface = np.zeros((nmax,nmax-1))
v[:,:] = 0
for i in range(nmax-1):
    for j in range(nmax):
        uface[i,j] = (p[i,j] - p[i+1,j])
for i in range(nmax):
    for j in range(nmax-1):
        vface[i,j] = (p[i,j] - p[i,j+1])
        vface[i,j] = 0

momx = frac * u
momy = frac * v        



def RungeKutta(q,n):
    
    qdx = np.zeros((nmax,nmax))
    qdy = np.zeros((nmax,nmax))
    
    fluxfaceleft = np.zeros((nmax-1,nmax))
    fluxfaceright = np.zeros((nmax-1,nmax))
    fluxfacetop = np.zeros((nmax,nmax-1))
    fluxfacebottom = np.zeros((nmax,nmax-1))
    
    fluxfacex = np.zeros((nmax-1,nmax))
    fluxfacey = np.zeros((nmax,nmax-1))
    
    res = np.zeros((nmax,nmax))
    
    for i in range(1,nmax-1):
        for j in range(0,nmax):
            qc = (q[i-1,j] - q[i+1,j])/2
            qr = q[i,j] - q[i+1,j]
            ql = q[i-1,j] - q[i,j]
            qtemp = 2 * qr
            if abs(ql) < abs(qr):
                qtemp = 2 * ql
            if abs(qc) < abs(qtemp):
                qtemp = qc
            if ql * qr <= 0:
                qtemp = 0
            qdx[i,j] = qtemp
    qdx[0,:] = qdx[1,:]
    qdx[nmax-1,:] = qdx[nmax-2,:]
    for i in range(0,nmax):
        for j in range(1,nmax-1):
            qc = (q[i,j-1] - q[i,j+1])/2
            qr = q[i,j] - q[i,j+1]
            ql = q[i,j-1] - q[i,j]
            qtemp = 2 * qr
            if abs(ql) < abs(qr):
                qtemp = 2 * ql
            if abs(qc) < abs(qtemp):
                qtemp = qc
            if ql * qr <= 0:
                qtemp = 0
            qdy[i,j] = qtemp
    qdy[:,0] = qdy[:,1]
    qdy[:,nmax-1] = qdy[:,nmax-2]
    
 
    
    
    for i in range(1,nmax-1):
        for j in range(0,nmax):
            fluxfaceleft[i,j] = q[i,j] + qdx[i,j]/2
            if (normalx[i,j] > 0) & (normalx[i,j] < 1) :
                if normaldir[i,j] == 1:
                    fluxfaceleft[i,j] = 0
            fluxfaceright[i,j] = q[i+1,j] - qdx[i+1,j]/2
            if (normalx[i+1,j] > 0) & (normalx[i+1,j] < 1) :
                if normaldir[i+1,j] == 1:
                    fluxfaceright[i,j] = q[i-1,j] + qdx[i-1,j]*3/2
            fluxfacex[i,j] = uface[i,j]*(fluxfaceleft[i,j] + fluxfaceright[i,j])/2
            if (normalx[i,j] < 0.99) & (normalx[i+1,j] == 0):
                fluxfacex[i,j] = 0
            if(normalx[i,j] == 1) & (normalx[i+1,j] == 0):
                fluxfacex[i,j] = q[i-1,j] + qdx[i-1,j]*3/2
            
    for j in range(nmax):
        fluxfacex[0,j] = 0
        if (q[0,j] <= 1) & (q[1,j] > 0) & (q[1,j] < 1) & (q[0,j] > 0):
            fluxfacex[0,j] = 1
        fluxfacex[nmax-2,j] = 0
        if q[nmax-2,j] == 1:
            fluxfacex[i,j] = 1
    for i in range(0,nmax):
        for j in range(1,nmax-1):
            fluxfacetop[i,j] = q[i,j] + qdy[i,j]/2
            fluxfacebottom[i,j] = q[i,j+1] - qdy[i,j+1]/2
            fluxfacey[i,j] = vface[i,j]*(fluxfacetop[i,j] + fluxfacebottom[i,j])/2
    for i in range(nmax):
        fluxfacey[i,0] = 0
        if q[i,0] == 1:
            fluxfacey[i,0] = vface[i,0] * 1
    
    if n == 31:
        n = 31
    
    for i in range(nmax-1):
        for j in range(nmax):
            res[i,j] -= fluxfacex[i,j]
            res[i+1,j] += fluxfacex[i,j]
    for i in range(nmax):
        for j in range(nmax-1):
            res[i,j] -= fluxfacey[i,j]
            res[i,j+1] += fluxfacey[i,j]
            
    return res
    


dt = 0.1
time = 0
timeFinal = 200

ft = np.zeros((timeFinal+1,nmax,nmax))
ft[0,:,:] = frac[:,:]
nt = np.zeros((timeFinal+1,nmax,nmax))
nt[0,:,:] = normalx[:,:]
ut = np.zeros((timeFinal+1,nmax-1,nmax))
ut[0,:,:] = uface[:,:]
while(time < timeFinal):
    
    
    dfrac = RungeKutta(frac,time)
    frac += dt * dfrac
    test = 1
    for i in range(nmax-1):
        for j in range(nmax-1):
            if frac[i,j] < 1:
                lack = min(1 - frac[i,j],frac[i,j+1])
                lack = min(lack,0.3)
                frac[i,j] += lack
                frac[i,j+1] -= lack
            normalx[i,j] = frac[i,j]
            normaldir[i,j] = 1
            
    for i in range(nmax):
        for j in range(nmax-1):
            if frac[i,j] < 1e-8:
                frac[i,j] = 0
            if frac[i,j] > 1 - 1e-8:
                frac[i,j] = 1
            
    for i in range(nmax-1):
        for j in range(nmax):
            flag0 = (frac[i,j] - frac[i+1,j]) > 1e-8
            if flag0:
                uface[i,j] = 1
            else:
                uface[i,j] = 0
            
            
    vface[nmax-1,:] = 1
            
    
            
    
    ft[time+1,:,:] = frac[:,:]
    nt[time+1,:,:] = normalx[:,:]
    ut[time+1,:,:] = uface[:,:]
    time = time + 1
        

        
    