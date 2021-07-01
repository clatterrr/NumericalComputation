import numpy as np
# pyro compressible rt instability
# D:\FluidSim\MiscOpenSource\pyro2-main\pyro2-main
dens = np.load('dens.npy')
ener = np.load('ener.npy')
ymom = np.load('ymom.npy')
xmom = np.zeros((72,200))
grav = -1
dt = 3.331667915625911e-05
e = np.zeros((72,200))
U = np.zeros((72,200,4))
U[:,:,0] = dens[:,:]
U[:,:,1] = ener[:,:]
U[:,:,2] = xmom[:,:]
U[:,:,3] = ymom[:,:]
Q = np.zeros((4,72,200))
ldx = np.zeros((4,72,200))
ldy = np.zeros((4,72,200))
ldxtemp = np.zeros((4,72,200))
ldytemp = np.zeros((4,72,200))
nmax = 72
mmax = 200
gamma = 1.4
for t in range(0,1):
    # cfl = 0.8
    # xtmp = griddx/(abs(xmom) + cs)
    # ytmp = griddx/(abs(xmom) + cs)
    # dt = cfl*float(min(xtmp.min(), ytmp.min()))
    Q[0,:,:] = U[:,:,0]
    Q[1,:,:] = xmom[:,:] / dens[:,:]
    Q[2,:,:] = ymom[:,:] / dens[:,:]
    e[:,:] = (ener[:,:] - 0.5*Q[0,:,:]*(Q[1,:,:]**2 + Q[2,:,:]**2))/dens[:,:]
    Q[3,:,:] = dens[:,:]*e[:,:]*(gamma - 1)
    
    for k in range(4):
        
        for i in range(1,nmax-1):
            for j in range(0,mmax):
                dc = (Q[k,i+1,j] - Q[k,i-1,j])/2
                dl = Q[k,i,j] - Q[k,i-1,j]
                dr = Q[k,i+1,j] - Q[k,i,j]
                temp = 2*dl
                if abs(dr) < abs(dl):
                    temp  = 2*dr
                if abs(dc) < abs(temp):
                    temp = dc
                if dl * dr <= 0:
                    temp = 0
                ldxtemp[k,i,j] = temp
        ldxtemp[k,0,:] = ldxtemp[k,1,:]
        ldxtemp[k,nmax-1,:] = ldxtemp[k,nmax-2,:]
        for i in range(1,nmax-1):
            for j in range(0,mmax):
                dc = (Q[k,i+1,j] - Q[k,i-1,j] - (ldxtemp[k,i+1,j] + ldxtemp[k,i-1,j])/4)*2/3
                dl = Q[k,i,j] - Q[k,i-1,j]
                dr = Q[k,i+1,j] - Q[k,i,j]
                temp = 2*dl
                if abs(dr) < abs(dl):
                    temp  = 2*dr
                if abs(dc) < abs(temp):
                    temp = dc
                if dl * dr <= 0:
                    temp = 0
                ldx[k,i,j] = temp
                
                
        for i in range(0,nmax):
            for j in range(1,mmax-1):
                dc = (Q[k,i,j+1] - Q[k,i,j-1])/2
                dl = Q[k,i,j] - Q[k,i,j-1]
                dr = Q[k,i,j+1] - Q[k,i,j]
                temp = 2*dl
                if abs(dr) < abs(dl):
                    temp  = 2*dr
                if abs(dc) < abs(temp):
                    temp = dc
                if dl * dr <= 0:
                    temp = 0
                ldytemp[k,i,j] = temp
        ldytemp[k,:,0] = ldytemp[k,:,1]
        ldytemp[k,:,nmax-1] = ldytemp[k,:,nmax-2]
        for i in range(0,nmax):
            for j in range(1,mmax-1):
                dc = (Q[k,i,j+1] - Q[k,i,j-1] - (ldytemp[k,i,j+1] + ldytemp[k,i,j-1])/4)*2/3
                dl = Q[k,i,j] - Q[k,i,j-1]
                dr = Q[k,i,j+1] - Q[k,i,j]
                temp = 2*dl
                if abs(dr) < abs(dl):
                    temp  = 2*dr
                if abs(dc) < abs(temp):
                    temp = dc
                if dl * dr <= 0:
                    temp = 0
                ldy[k,i,j] = temp

    nvar, qx, qy = Q.shape

    q_l = np.zeros_like(Q)
    q_r = np.zeros_like(Q)
    
    dt = 3.331667915625911e-05
    dx = 0.015625
    
    ng = 4
    nspec = 0
    nx = qx - 2 * ng
    ny = qy - 2 * ng
    ilo = ng
    ihi = ng + nx
    jlo = ng
    jhi = ng + ny

    ns = nvar - nspec

    dtdx = dt / dx
    dtdx4 = 0.25 * dtdx

    lvec = np.zeros((nvar, nvar))
    rvec = np.zeros((nvar, nvar))
    e_val = np.zeros(nvar)
    betal = np.zeros(nvar)
    betar = np.zeros(nvar)
    
    for i in range(ilo - 2,ihi + 2):
        for j in range(jlo-2,jhi+2):
            
            dq = ldx[:,i,j]
            q = Q[:,i,j]
            
            cs = np.sqrt(gamma * q[3] / q[0])
            lvec[:,:] = 0
            rvec[:,:] = 0
            e_val[:] = 0
            
            e_val[:] = np.array([q[1] - cs, q[1], q[1], q[1] + cs])

            lvec[0, :ns] = [0.0, -0.5 *
                                 q[0] / cs, 0.0, 0.5 / (cs * cs)]
            lvec[1, :ns] = [1.0, 0.0,
                                 0.0, -1.0 / (cs * cs)]
            lvec[2, :ns] = [0.0, 0.0,             1.0, 0.0]
            lvec[3, :ns] = [0.0, 0.5 *
                                 q[0] / cs,  0.0, 0.5 / (cs * cs)]

            rvec[0, :ns] = [1.0, -cs / q[0], 0.0, cs * cs]
            rvec[1, :ns] = [1.0, 0.0,       0.0, 0.0]
            rvec[2, :ns] = [0.0, 0.0,       1.0, 0.0]
            rvec[3, :ns] = [1.0, cs / q[0],  0.0, cs * cs]

            # now the species -- they only have a 1 in their corresponding slot
            e_val[ns:] = q[1]
            
            ix = -1
            for n in range(ix, ix + nspec):
                lvec[n, n] = 1.0
                rvec[n, n] = 1.0
                
            # this is one the right face of the current zone,
            # so the fastest moving eigenvalue is e_val[3] = u + c
            factor = 0.5 * (1.0 - dtdx * max(e_val[3], 0.0))
            q_l[:,i + 1, j] = q + factor * dq
            
            # compute the Vhat functions
            for m in range(nvar):
                sum = np.dot(lvec[m, :], dq)

                betal[m] = dtdx4 * (e_val[3] - e_val[m]) * \
                    (np.copysign(1.0, e_val[m]) + 1.0) * sum
                betar[m] = dtdx4 * (e_val[0] - e_val[m]) * \
                    (1.0 - np.copysign(1.0, e_val[m])) * sum

            # construct the states
            for m in range(nvar):
                sum_l = np.dot(betal, rvec[:, m])
                sum_r = np.dot(betar, rvec[:, m])

                q_l[m,i + 1, j] = q_l[m,i + 1,j] + sum_l
                q_r[m,i,j] = q_r[m,i,j] + sum_r
            
       