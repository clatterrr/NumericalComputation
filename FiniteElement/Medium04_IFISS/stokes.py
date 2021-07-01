import numpy as np
import scipy.io as scio
# https://personalpages.manchester.ac.uk/staff/david.silvester/ifiss/
nngpt = 9 # 每个网格几个点
if nngpt == 4:
    gpt = np.sqrt(1 / 3)
    gaussxi = np.array([-gpt,gpt,gpt,-gpt])
    gausseta = np.array([-gpt,-gpt,gpt,gpt])
    weight = np.array([1,1,1,1])
elif nngpt == 9:
    gpt = np.sqrt(3 / 5)
    gaussxi = np.array([-gpt,gpt,gpt,-gpt,0,gpt,0,-gpt,0])
    gausseta = np.array([-gpt,-gpt,gpt,gpt,-gpt,0,gpt,0,0])
    weight = np.array([25,25,25,25,40,40,40,40,64]) / 81
    
numgrid = 32 # 每行多少个格子
numgrid2 = numgrid*numgrid
numNodes = 65
gridindex = np.zeros((numgrid2,9))
gridmp = np.zeros((numgrid2,4))
temp = scio.loadmat('x.mat')['x']
x = temp[:,0]
y = temp[:,0]
gridcoordx = np.zeros((numgrid2,9))
gridcoordy = np.zeros((numgrid2,9))
for j in range(numgrid):
    for i in range(numgrid):
        st = j * numNodes * 2 + i * 2
        idx = i * numgrid + j
        gridindex[idx,:] = np.array([st,st+2,
                    st+numNodes*2+2,
                    st+numNodes*2,st+1,st+numNodes+2,
                    st+numNodes*2+1,st+numNodes,
                    st+numNodes+1])
        st = j * (numgrid + 1) + i
        gridmp[idx,:] = np.array([st,st+1,st+numgrid+2,st+numgrid+1])
        stx = j * 2
        sty = i * 2
        gridcoordx[idx,:] = np.array([x[stx],x[stx+2],x[stx+2],x[stx],
                               x[stx+1],x[stx+2],x[stx+1],x[stx],x[stx+1]])
        gridcoordy[idx,:] = np.array([y[sty],y[sty],y[sty+2],y[sty+2],
                               y[sty],y[sty+1],y[sty+2],y[sty+1],y[sty+1]])

ae = np.zeros((numgrid2,9,9))
re = np.zeros((numgrid2,9,9))
bbxe = np.zeros((numgrid2,9,9))
bbye = np.zeros((numgrid2,9,9))
bxe = np.zeros((numgrid2,4,9))
bye = np.zeros((numgrid2,4,9))
mpe = np.zeros((numgrid2,4,4))
ge = np.zeros((numgrid2,4))

for igpt in range(nngpt):
    xi = gaussxi[igpt]
    eta = gausseta[igpt]
    w = weight[igpt]
    
    # Deriv
    
    # evaluate bilinear shape functions
    phie = np.array([(xi - 1)*(eta - 1)/4,-(xi + 1)*(eta - 1)/4,
                    (xi + 1)*(eta + 1)/4,-(xi - 1)*(eta + 1)/4])
    phidxi = np.array([(eta-1)/4,-(eta-1)/4,(eta+1)/4,-(eta+1)/4])
    phideta = np.array([(xi-1)/4,-(xi+1)/4,(xi+1)/4,-(xi-1)/4])

    ivtx = 4
    dxdxi = np.zeros((numgrid2))
    dxdeta = np.zeros((numgrid2))
    dydxi = np.zeros((numgrid2))
    dydeta = np.zeros((numgrid2))
    for i in range(ivtx):
        dxdxi[:] = dxdxi[:] + gridcoordx[:,i] * phidxi[i]
        dxdeta[:] = dxdeta[:] + gridcoordx[:,i] * phideta[i]
        dydxi[:] = dydxi[:] + gridcoordy[:,i] * phidxi[i]
        dydeta[:] = dydeta[:] + gridcoordy[:,i] * phideta[i]
    jac = dxdxi * dydeta - dxdeta * dydxi
    
    jacinv = 1 / jac
    phi = np.zeros((numgrid2,ivtx))
    dphidx = np.zeros((numgrid2,ivtx))
    dphidy = np.zeros((numgrid2,ivtx))
    for i in range(ivtx):
        phi[:,i] = phie[i]
        dphidx[:,i] = phidxi[i]*dydeta[:] - phideta[i]*dydxi[:]
        dphidy[:,i] = (-phidxi[i]*dxdeta[:] + phideta[i]*dxdxi[:])
        
    # evaluate biquadratic shape functions
    ellx = np.array([xi*(xi-1)/2,1-xi*xi,xi*(xi+1)/2])
    elly = np.array([eta*(eta-1)/2,1-eta*eta,eta*(eta+1)/2])
    dellx = np.array([xi-1/2,-2*xi,xi+1/2])
    delly = np.array([eta-1/2,-2*eta,eta+1/2])
    psie = np.array([ellx[0]*elly[0],ellx[2]*elly[0],ellx[2]*elly[2],
                    ellx[0]*elly[2],ellx[1]*elly[0],ellx[2]*elly[1],
                    ellx[1]*elly[2],ellx[0]*elly[1],ellx[1]*elly[1]])
    dpsidxi = np.array([dellx[0]*elly[0],
                        dellx[2]*elly[0],
                        dellx[2]*elly[2],
                        dellx[0]*elly[2],
                        dellx[1]*elly[0],
                        dellx[2]*elly[1],
                        dellx[1]*elly[2],
                        dellx[0]*elly[1],
                        dellx[1]*elly[1]])
    dpsideta = np.array([ellx[0]*delly[0],
                        ellx[2]*delly[0],
                        ellx[2]*delly[2],
                        ellx[0]*delly[2],
                        ellx[1]*delly[0],
                        ellx[2]*delly[1],
                        ellx[1]*delly[2],
                        dellx[0]*delly[1],
                        ellx[1]*delly[1]])
    
    psi = np.zeros((numgrid2,9))
    dpsidx = np.zeros((numgrid2,9))
    dpsidy = np.zeros((numgrid2,9))
    for i in range(9):
        psi[:,i] = psie[i]
        dpsidx[:,i] = dpsidxi[i]*dydeta[:] - dpsideta[i]*dydxi[:]
        dpsidy[:,i] = (-dpsidxi[i]*dxdeta[:] + dpsideta[i]*dxdxi[:])
    

    for j in range(9):
        for i in range(9):
            ae[:,i,j] += w * dpsidx[:,i] * dpsidx[:,j] * jacinv
            ae[:,i,j] += w * dpsidy[:,i] * dpsidy[:,j] * jacinv
            re[:,i,j] += w * psi[:,i] * psi[:,j] * jac
            bbxe[:,i,j] += - w * psi[:,i] * dpsidx[:,j]
            bbye[:,i,j] += - w * psi[:,i] * dpsidy[:,j]
        for i in range(4):
            bxe[:,i,j] += -w * phi[:,i] * dpsidx[:,j]
            bye[:,i,j] += -w * phi[:,i] * dpsidy[:,j]
    for i in range(4):
        for j in range(4):
            mpe[:,i,j] += w * phi[:,i] * phi[:,j] * jac
            
nu = 8450
nvtx = nu / 2
a = np.zeros((nu,nu))
r = np.zeros((nu,nu))
bbx = np.zeros((nvtx,nvtx))
bby = np.zeros((nvtx,nvtx))
bx = np.zeros((nvtx,nvtx))
by = np.zeros((nvtx,nvtx))
for krow in range(9):
    nrow = gridindex[:,krow]
    for kcol in range(9):
        ncol = gridindex[:,kcol]
        for i0 in range(numgrid2):
            a[nrow[i0],ncol[i0]] += ae[i0]
            a[nrow[i0] + nvtx,ncol[i0] + nvtx] += ae[i0]
            r[nrow[i0],ncol[i0]] += re[i0]
            r[nrow[i0] + nvtx,ncol[i0] + nvtx] += re[i0]
            bbx[nrow[i0],ncol[i0]] += bbxe[i0]
            bby[nrow[i0],ncol[i0]] += bbye[i0]
    for kcol in range(4):
        ncol = gridmp[:,kcol]
        for i0 in range(numgrid2):
            bx[nrow[i0],ncol[i0]] += bxe[i0]
            by[nrow[i0],ncol[i0]] += bye[i0]
m = np.zeros((1089,1089))
for krow in range(4):
    nrow = gridmp[:,krow]
    for kcol in range(4):
        ncol = gridmp[:,kcol]
        for i0 in range(numgrid2):
            m[nrow[i0],ncol[i0]] += mpe[i0]
        
time = 0
timeFinal = 1
while(time < timeFinal):
    
    flowsol = np.zeros((nu))
    usol = flowsol[0:nu//2]
    vsol = flowsol[0:nu//2]
    gpt = np.sqrt(3 / 5)
    gaussxi = np.array([-gpt,gpt,gpt,-gpt,0,gpt,0,-gpt,0])
    gausseta = np.array([-gpt,-gpt,gpt,gpt,-gpt,0,gpt,0,0])
    weight = np.array([25,25,25,25,40,40,40,40,64]) / 81
    nnpgt = 9
    
    
    ne = np.zeros((numgrid2,9,9))
    for igpt in range(nnpgt):
        xi = gaussxi[igpt]
        eta = gausseta[igpt]
        w = weight[igpt]
        
        # Deriv
        
        # evaluate bilinear shape functions
        phie = np.array([(xi - 1)*(eta - 1)/4,-(xi + 1)*(eta - 1)/4,
                        (xi + 1)*(eta + 1)/4,-(xi - 1)*(eta + 1)/4])
        phidxi = np.array([(eta-1)/4,-(eta-1)/4,(eta+1)/4,-(eta+1)/4])
        phideta = np.array([(xi-1)/4,-(xi+1)/4,(xi+1)/4,-(xi-1)/4])
    
        ivtx = 4
        dxdxi = np.zeros((numgrid2))
        dxdeta = np.zeros((numgrid2))
        dydxi = np.zeros((numgrid2))
        dydeta = np.zeros((numgrid2))
        for i in range(ivtx):
            dxdxi[:] = dxdxi[:] + gridcoordx[:,i] * phidxi[i]
            dxdeta[:] = dxdeta[:] + gridcoordx[:,i] * phideta[i]
            dydxi[:] = dydxi[:] + gridcoordy[:,i] * phidxi[i]
            dydeta[:] = dydeta[:] + gridcoordy[:,i] * phideta[i]
        jac = dxdxi * dydeta - dxdeta * dydxi
        
        jacinv = 1 / jac
        phi = np.zeros((numgrid2,ivtx))
        dphidx = np.zeros((numgrid2,ivtx))
        dphidy = np.zeros((numgrid2,ivtx))
        for i in range(ivtx):
            phi[:,i] = phie[i]
            dphidx[:,i] = phidxi[i]*dydeta[:] - phideta[i]*dydxi[:]
            dphidy[:,i] = (-phidxi[i]*dxdeta[:] + phideta[i]*dxdxi[:])
            
        # evaluate biquadratic shape functions
        ellx = np.array([xi*(xi-1)/2,1-xi*xi,xi*(xi+1)/2])
        elly = np.array([eta*(eta-1)/2,1-eta*eta,eta*(eta+1)/2])
        dellx = np.array([xi-1/2,-2*xi,xi+1/2])
        delly = np.array([eta-1/2,-2*eta,eta+1/2])
        psie = np.array([ellx[0]*elly[0],ellx[2]*elly[0],ellx[2]*elly[2],
                        ellx[0]*elly[2],ellx[1]*elly[0],ellx[2]*elly[1],
                        ellx[1]*elly[2],ellx[0]*elly[1],ellx[1]*elly[1]])
        dpsidxi = np.array([dellx[0]*elly[0],
                            dellx[2]*elly[0],
                            dellx[2]*elly[2],
                            dellx[0]*elly[2],
                            dellx[1]*elly[0],
                            dellx[2]*elly[1],
                            dellx[1]*elly[2],
                            dellx[0]*elly[1],
                            dellx[1]*elly[1]])
        dpsideta = np.array([ellx[0]*delly[0],
                            ellx[2]*delly[0],
                            ellx[2]*delly[2],
                            ellx[0]*delly[2],
                            ellx[1]*delly[0],
                            ellx[2]*delly[1],
                            ellx[1]*delly[2],
                            dellx[0]*delly[1],
                            ellx[1]*delly[1]])
        
        psi = np.zeros((numgrid2,9))
        dpsidx = np.zeros((numgrid2,9))
        dpsidy = np.zeros((numgrid2,9))
        for i in range(9):
            psi[:,i] = psie[i]
            dpsidx[:,i] = dpsidxi[i]*dydeta[:] - dpsideta[i]*dydxi[:]
            dpsidy[:,i] = (-dpsidxi[i]*dxdeta[:] + dpsideta[i]*dxdxi[:])
        
        u_x = np.zeros((numgrid2))
        u_y = np.zeros((numgrid2))
        for k in range(9):
            u_x[:] += psi[:] * gridcoordx[:,0:4]
            u_y[:] += psi[:] * gridcoordy[:,0:4]
        for j in range(9):
            for i in range(9):
                ne[:,i,j] += w * u_x[:] * psi[:,i] * dpsidx[:,j]
                ne[:,i,j] += w * u_y[:] * psi[:,i] * dpsidy[:,j]
        nn = np.zeros((nvtx,nvtx))
        for krow in range(9):
            nrow = gridindex[:,krow]
            for kcol in range(9):
                ncol = gridindex[:,kcol]
                for i0 in range(numgrid2):
                    nn[nrow[i0],ncol[i0]] += ne[i0]
            
        G = np.zeros((nu,nu))
        A = np.zeros((nu,nu)) # vector diffusion matrix
        B = np.zeros((nu,nu)) # divergence matrix
        u = np.zeros((nu,nu))
        udot = np.zeros((nu,nu))
        N = np.zeros((nu//2,nu//2))
        sparseN = np.zeros((nu,nu))
        sparseN[0:nu//2,0:nu//2] = sparseN[nu//2:nu,nu//2:nu] = N[:,:]
        dt = 1
        viscosity = 1
        Anst = 2*G + dt*viscosity*A + dt*sparseN
        fnst = G*udot -(viscosity*A + sparseN)*u;
        gnst = -(B*u);
        
        # w = udot + (.5*dt)*udd;
        # udiff = v - w;
        # uAB2  = u + dt*w;
        # upTR  = u + dt*v;
        # % local truncation error estimate
        # d = (dt^2/(3*(dt+dt0)))*sqrt((udiff'*(G*udiff)));
        # % acceleration
        # acc = sqrt((udot'*(G*udot)));
        # time step rejection