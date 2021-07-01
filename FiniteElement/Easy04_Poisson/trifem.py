import numpy as np
# https://github.com/EnigmaHuang/Poisson2D_FEM
nmax = 11
numgrid = nmax - 1
coord = np.zeros((nmax*nmax,2))
element = np.zeros((numgrid * numgrid * 2,3),dtype = int)
for i in range(nmax):
    for j in range(nmax):
        idx = int(i * nmax + j)
        coord[idx,0] = j / 10
        coord[idx,1] = i / 10
for i in range(numgrid):
    for j in range(numgrid):
        idx = int(i * numgrid + j)*2
        st = int(i * nmax + j)
        element[idx,0] = st
        element[idx,1] = st+1
        element[idx,2] = st+1+nmax
        element[idx+1,0] = st
        element[idx+1,1] = st+1+nmax
        element[idx+1,2] = st+nmax
Kmat = np.zeros((nmax*nmax,nmax*nmax))

def femalpha(x,i0,j0):
    qx = np.array([-0.3399810435848563, 0.3399810435848563, -0.8611363115940526, 0.8611363115940526])
    qw = np.array([ 0.6521451548625461, 0.6521451548625461,  0.3478548451374538, 0.3478548451374538])
    n_quadrature = 4
    dx = x[j0,0] - x[i0,0]
    dy = x[j0,1] - x[i0,1]
    le = np.sqrt(dx*dx + dy*dy)
    semibma = le / 2
    semibpa = le / 2
    k2 = np.zeros((3,3))
    fii = 0
    fij = 0
    fjj = 0
    for iq in range(n_quadrature):
        t = semibma * qx[iq] + semibpa
        toverl = t / le
        alphax = pos[i0,0] + toverl * dx
        alphay = pos[i0,1] + toverl * dy;
        alpha   = 1
        
        fii = fii + ((1 - toverl) * (1 - toverl) * alpha) * qw[iq];
        fij = fij + ((1 - toverl) * toverl       * alpha) * qw[iq];
        fjj = fjj + (toverl * toverl * alpha) * qw[iq];
    
    k2[i0,i0] = semibma * fii
    k2[i0,j0] = semibma * fij
    k2[j0,i0] = semibma * fij
    k2[j0,j0] = semibma * fjj
    return k2

def fembeta(x,i0,j0):
    qx = np.array([-0.3399810435848563, 0.3399810435848563, -0.8611363115940526, 0.8611363115940526])
    qw = np.array([ 0.6521451548625461, 0.6521451548625461,  0.3478548451374538, 0.3478548451374538])
    n_quadrature = 4
    dx = x[j0,0] - x[i0,0]
    dy = x[j0,1] - x[i0,1]
    le = np.sqrt(dx*dx + dy*dy)
    semibma = le / 2
    semibpa = le / 2
    b2 = np.zeros((3))
    bi = 0
    bj = 0
    for iq in range(n_quadrature):
        t = semibma * qx[iq] + semibpa
        toverl = t / le
        gx = x[i0,0] + toverl * dx
        gy = x[i0,1] + toverl * dy
        gxy = gx + gy
        
        bi = bi + ((1-toverl)*gxy)*qw[iq]
        bj = bj + (toverl*gxy)*qw[iq]
    
    b2[i0] = semibma * bi
    b2[j0] = semibpa * bj
    return b2


for k in range(element.shape[0]):
    
    ele = element[k,:]
    pos = coord[element[k,:],:]
    
    ktemp = np.zeros((3,3))
    
    gaussxi = np.array([1/3,0.5,0.2,0.2])
    gausseta = np.array([1/3,0.2,0.2,0.6])
    w = np.array([-27/96,25/96,25/96,25/96])
    n_quadrature = 4
    
    for i in range(n_quadrature):
            xi = gaussxi[i]
            eta = gausseta[i]
            d_N1_d_xi = -1
            d_N1_d_eta = -1
            d_N2_d_xi = 1
            d_N2_d_eta = 0
            d_N3_d_xi = 0
            d_N3_d_eta = 1
            
            
            d_x_d_xi  = d_N1_d_xi  * pos[0,0] + d_N2_d_xi  * pos[1,0] + d_N3_d_xi  * pos[2,0]
            d_x_d_eta = d_N1_d_eta * pos[0,0] + d_N2_d_eta * pos[1,0] + d_N3_d_eta * pos[2,0]
            d_y_d_xi  = d_N1_d_xi  * pos[0,1] + d_N2_d_xi  * pos[1,1] + d_N3_d_xi  * pos[2,1]
            d_y_d_eta = d_N1_d_eta * pos[0,1] + d_N2_d_eta * pos[1,1] + d_N3_d_eta * pos[2,1]
            
            jac = d_x_d_xi * d_y_d_eta - d_y_d_xi * d_x_d_eta;
            
            dN = np.array([[(d_N1_d_xi *  d_y_d_eta + d_N1_d_eta * -d_y_d_xi) / jac,
                            (d_N2_d_xi *  d_y_d_eta + d_N2_d_eta * -d_y_d_xi) / jac,
                            (d_N3_d_xi *  d_y_d_eta + d_N3_d_eta * -d_y_d_xi) / jac],
                           [(d_N1_d_xi * -d_x_d_eta + d_N1_d_eta *  d_x_d_xi) / jac,
                            (d_N2_d_xi * -d_x_d_eta + d_N2_d_eta *  d_x_d_xi) / jac,
                            (d_N3_d_xi * -d_x_d_eta + d_N3_d_eta *  d_x_d_xi) / jac]])
            
            ktemp = ktemp + np.dot(np.transpose(dN),dN) * jac * w[i]
    ktemp2 = np.zeros((3,3))
    if(((pos[0,0] == 0) & (pos[1,0] == 0)) |
       ((pos[0,0] == 1) & (pos[1,0] == 1)) | 
       ((pos[0,1] == 0) & (pos[1,1] == 0)) | 
       ((pos[0,1] == 1) & (pos[1,1] == 1))):
        ktemp2 = ktemp2 + femalpha(pos,0,1)
                
    if(((pos[1,0] == 0) & (pos[2,0] == 0)) |
       ((pos[1,0] == 1) & (pos[2,0] == 1)) | 
       ((pos[1,1] == 0) & (pos[2,1] == 0)) | 
       ((pos[1,1] == 1) & (pos[2,1] == 1))):
        ktemp2 = ktemp2 + femalpha(pos,1,2)
                
    if(((pos[2,0] == 0) & (pos[0,0] == 0)) |
       ((pos[2,0] == 1) & (pos[0,0] == 1)) | 
       ((pos[2,1] == 0) & (pos[0,1] == 0)) | 
       ((pos[2,1] == 1) & (pos[0,1] == 1))):
        ktemp2 = ktemp2 + femalpha(pos,2,0)
            
    ktemp = ktemp + ktemp2
    test = 1
    
    for i in range(3):
        for j in range(3):
            Kmat[ele[i],ele[j]] += ktemp[i,j] 
            
BMat = np.zeros((nmax * nmax))
for k in range(element.shape[0]):
    ele = element[k,:]
    pos = coord[element[k,:],:]
    
    gaussxi = np.array([1/3,0.5,0.2,0.2])
    gausseta = np.array([1/3,0.2,0.2,0.6])
    w = np.array([-27/96,25/96,25/96,25/96])
    n_quadrature = 4
    btemp = np.zeros((3))
    for i in range(n_quadrature):
            xi = gaussxi[i]
            eta = gausseta[i]
            
            fx = pos[0,0] * (1 - xi - eta) + pos[1,0] * xi + pos[2,0] * eta
            fy = pos[0,1] * (1 - xi - eta) + pos[1,1] * xi + pos[2,1] * eta 
            
            fxy = 5 * fx * fy # 泊松方程右边项
            
            d_x_d_xi = pos[1,0] - pos[0,0]
            d_x_d_eta = pos[2,0] - pos[0,0]
            d_y_d_xi = pos[1,1] - pos[0,1]
            d_y_d_eta = pos[2,1] - pos[0,1]
            
            jac = d_x_d_xi * d_y_d_eta - d_x_d_eta * d_y_d_xi
            
            btemp = btemp + np.array([
                jac * fxy * (1 - xi - eta) * w[i],
                jac * fxy * xi * w[i],
                jac * fxy * eta * w[i]])
    
    btemp2 = np.zeros((3))
    if(((pos[0,0] == 0) & (pos[1,0] == 0)) |
       ((pos[0,0] == 1) & (pos[1,0] == 1)) | 
       ((pos[0,1] == 0) & (pos[1,1] == 0)) | 
       ((pos[0,1] == 1) & (pos[1,1] == 1))):
        btemp2 = btemp2 + fembeta(pos,0,1)
                
    if(((pos[1,0] == 0) & (pos[2,0] == 0)) |
       ((pos[1,0] == 1) & (pos[2,0] == 1)) | 
       ((pos[1,1] == 0) & (pos[2,1] == 0)) | 
       ((pos[1,1] == 1) & (pos[2,1] == 1))):
        btemp2 = btemp2 + fembeta(pos,1,2)
                
    if(((pos[2,0] == 0) & (pos[0,0] == 0)) |
       ((pos[2,0] == 1) & (pos[0,0] == 1)) | 
       ((pos[2,1] == 0) & (pos[0,1] == 0)) | 
       ((pos[2,1] == 1) & (pos[0,1] == 1))):
        btemp2 = btemp2 + fembeta(pos,2,0)
            
    btemp = btemp + btemp2
    test = 1
    
    for i in range(3):
        BMat[ele[i]] += btemp[i]