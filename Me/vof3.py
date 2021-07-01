import numpy as np

nmax = 19

frac = np.zeros((nmax,nmax))
momx = np.zeros((nmax,nmax))
momy = np.zeros((nmax,nmax))
p = np.zeros((nmax+2,nmax+2))
div = np.zeros((nmax,nmax))
volume = np.zeros((nmax,nmax))
normalx = np.zeros((nmax,nmax))
normaldir = np.zeros((nmax,nmax))
for i in range(10):
    for j in range(10):
        frac[i,j] = 1
        p[i,j] = (10 - j)/10 + (10 - i) / 100
volume[0:10,0:10] = 1
normalx[0:10,0:10] = 1

uface = np.zeros((nmax+2,nmax+1))
vface = np.zeros((nmax+1,nmax+2))
for i in range(nmax-1):
    for j in range(nmax):
        uface[i+1,j] = (p[i,j] - p[i+1,j])
for i in range(nmax):
    for j in range(nmax-1):
        vface[i,j+1] = (p[i,j] - p[i,j+1])  - 0.1
        
        
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
            # fluxfaceleft[i,j] = q[i,j] + qdx[i,j]/2
            # fluxfaceright[i,j] = q[i+1,j] - qdx[i+1,j]/2
            # fluxfacex[i,j] = uface[i+1,j]*(fluxfaceleft[i,j] + fluxfaceright[i,j])/2
            
            
            
            if (normalx[i,j] == 1) & (normalx[i+1,j] < 1):
                if normaldir[i,j] == 0:
                    fluxfacex[i,j] = uface[i+1,j] * q[i,j]
            elif(normalx[i,j] < 1) & (normalx[i,j] > 0) & (normalx[i+1,j] < 1) & (normalx[i+1,j] > 0):
                if normaldir[i,j] == 0:
                    fluxfacex[i,j] = uface[i+1,j] * (q[i,j] + q[i+1,j])/2
                    
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
            fluxfacey[i,j] = vface[i,j+1]*(fluxfacetop[i,j] + fluxfacebottom[i,j])/2
    for i in range(nmax):
        fluxfacey[i,0] = 0
        if q[i,0] == 1:
            fluxfacey[i,0] = vface[i,0] * 1
    
    if n == 1:
        test = 30
    
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
timeFinal = 500

ft = np.zeros((timeFinal+1,nmax,nmax))
ut = np.zeros((timeFinal+1,nmax+2,nmax+1))
ft[0,:,:] = frac[:,:]
ut[0,:,:] = uface[:,:]
while(time < timeFinal):
    
    
    dfrac = RungeKutta(frac,time)
    frac += dt * dfrac
    test = 1
    
    
    for i in range(nmax-1):
        for j in range(nmax-1):
            if frac[i,j] < 1:
                lack = min(1 - frac[i,j],frac[i,j+1])
                lack = min(lack,0.8)
                frac[i,j] += lack
                frac[i,j+1] -= lack
            
            
    eps = 1e-8
    for i in range(nmax):
        for j in range(nmax-1):
            if frac[i,j] < eps:
                frac[i,j] = 0
            if frac[i,j] > 1 - eps:
                frac[i,j] = 1
            normalx[i,j] = frac[i,j]
            
    for i in range(nmax-1):
        for j in range(nmax):
            if (normalx[i,j] == 1) & (normalx[i+1,j] < 1):
                uface[i+1,j] = 1
            elif(normalx[i,j] < 1) & (normalx[i,j] > 0) & (normalx[i+1,j] < 1) & (normalx[i+1,j] > 0):
                uface[i+1,j] = (frac[i,j] + frac[i+1,j])/2
    
    # for i in range(nmax):
    #     for j in range(nmax):
    #         div[i,j] = uface[i,j] - uface[i+1,j] + vface[i,j] - vface[i,j+1]
    # vface[nmax,:] = 1
    # for k in range(0,100):
    #     for i in range(1,nmax+1):
    #         for j in range(1,nmax+1):
    #             term = (p[i+1,j] + p[i-1,j]) + (p[i,j+1] + p[i,j-1]) - div[i-1,j-1]
    #             p[i,j] = p[i,j] + (term/4 - p[i,j])
    #     p[0,:] = p[1,:]
    #     p[nmax+1,:] = p[nmax,:]
    #     p[:,0] = p[:,1]
    #     p[:,nmax+1] = p[:,nmax]
    # for i in range(nmax+1):
    #     for j in range(nmax):
    #         uface[i,j] = uface[i,j] + (p[i,j] - p[i+1,j]) / 10
    uface[:,9] = 0
    
    ft[time+1,:,:] = frac[:,:]
    ut[time+1,:,:] = uface[:,:]
    time = time + 1
        

        
    