import numpy as np
# https://github.com/EnigmaHuang/Poisson2D_FEM
nmax = 11
numgrid = nmax - 1
coord = np.zeros((nmax*nmax,2))
element = np.zeros((numgrid * numgrid,4),dtype = int)
for i in range(nmax):
    for j in range(nmax):
        idx = int(i * nmax + j)
        coord[idx,0] = j / 10
        coord[idx,1] = i / 10
for i in range(numgrid):
    for j in range(numgrid):
        idx = int(i * numgrid + j)
        st = int(i * nmax + j)
        element[idx,0] = st
        element[idx,1] = st+1
        element[idx,2] = st+1+nmax
        element[idx,3] = st+nmax
Kmat = np.zeros((nmax*nmax,nmax*nmax))
Bmat = np.zeros((nmax*nmax))
def femalpha(x,i0,j0):
    qx = np.array([-0.3399810435848563, 0.3399810435848563, -0.8611363115940526, 0.8611363115940526])
    qw = np.array([ 0.6521451548625461, 0.6521451548625461,  0.3478548451374538, 0.3478548451374538])
    n_quadrature = 4
    dx = x[j0,0] - x[i0,0]
    dy = x[j0,1] - x[i0,1]
    le = np.sqrt(dx*dx + dy*dy)
    semibma = le / 2
    semibpa = le / 2
    k2 = np.zeros((4,4))
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
    b2 = np.zeros((4))
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
    
    ktemp = np.zeros((4,4))
    
    q = np.array([-1/np.sqrt(3),1/np.sqrt(3)])
    w = np.array([1,1])
    n_quadrature = 2
    
    for i in range(n_quadrature):
        for j in range(n_quadrature):
            qx = q[i]
            qy = q[j]
            
            d_N1_d_xi  = -0.25 * (1.0 - qy)
            d_N1_d_eta = -0.25 * (1.0 - qx)
            d_N2_d_xi  =  0.25 * (1.0 - qy)
            d_N2_d_eta = -0.25 * (1.0 + qx)
            d_N3_d_xi  =  0.25 * (1.0 + qy)
            d_N3_d_eta =  0.25 * (1.0 + qx)
            d_N4_d_xi  = -0.25 * (1.0 + qy)
            d_N4_d_eta =  0.25 * (1.0 - qx)
            
            d_x_d_xi  = d_N1_d_xi  * pos[0,0] + d_N2_d_xi  * pos[1,0] + d_N3_d_xi  * pos[2,0] + d_N4_d_xi  * pos[3,0]
            d_x_d_eta = d_N1_d_eta * pos[0,0] + d_N2_d_eta * pos[1,0] + d_N3_d_eta * pos[2,0] + d_N4_d_eta * pos[3,0]
            d_y_d_xi  = d_N1_d_xi  * pos[0,1] + d_N2_d_xi  * pos[1,1] + d_N3_d_xi  * pos[2,1] + d_N4_d_xi  * pos[3,1]
            d_y_d_eta = d_N1_d_eta * pos[0,1] + d_N2_d_eta * pos[1,1] + d_N3_d_eta * pos[2,1] + d_N4_d_eta * pos[3,1]
            
            jac = d_x_d_xi * d_y_d_eta - d_y_d_xi * d_x_d_eta
            
            dN = np.array([[(d_N1_d_xi *  d_y_d_eta + d_N1_d_eta * -d_y_d_xi) / jac,
                            (d_N2_d_xi *  d_y_d_eta + d_N2_d_eta * -d_y_d_xi) / jac,
                            (d_N3_d_xi *  d_y_d_eta + d_N3_d_eta * -d_y_d_xi) / jac,
                            (d_N4_d_xi *  d_y_d_eta + d_N4_d_eta * -d_y_d_xi) / jac],
                           [(d_N1_d_xi * -d_x_d_eta + d_N1_d_eta *  d_x_d_xi) / jac,
                            (d_N2_d_xi * -d_x_d_eta + d_N2_d_eta *  d_x_d_xi) / jac,
                            (d_N3_d_xi * -d_x_d_eta + d_N3_d_eta *  d_x_d_xi) / jac,
                            (d_N4_d_xi * -d_x_d_eta + d_N4_d_eta *  d_x_d_xi) / jac]])
            
            ktemp = ktemp + np.dot(np.transpose(dN),dN) * jac * w[i] * w[j]
            
    ktemp2 = np.zeros((4,4))
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
                
    if(((pos[2,0] == 0) & (pos[3,0] == 0)) |
       ((pos[2,0] == 1) & (pos[3,0] == 1)) | 
       ((pos[2,1] == 0) & (pos[3,1] == 0)) | 
       ((pos[2,1] == 1) & (pos[3,1] == 1))):
        ktemp2 = ktemp2 + femalpha(pos,2,3)
        
    if(((pos[3,0] == 0) & (pos[0,0] == 0)) |
       ((pos[3,0] == 1) & (pos[0,0] == 1)) | 
       ((pos[3,1] == 0) & (pos[0,1] == 0)) | 
       ((pos[3,1] == 1) & (pos[0,1] == 1))):
        ktemp2 = ktemp2 + femalpha(pos,3,0)
            
    ktemp = ktemp + ktemp2
    test = 1
    
    for i in range(4):
        for j in range(4):
            Kmat[ele[i],ele[j]] += ktemp[i,j] 
            
for k in range(element.shape[0]):
    
    ele = element[k,:]
    pos = coord[element[k,:],:]
    
    btemp = np.zeros((4))
    
    qx = np.array([-0.3399810435848563, 0.3399810435848563, -0.8611363115940526, 0.8611363115940526])
    qw = np.array([ 0.6521451548625461, 0.6521451548625461,  0.3478548451374538, 0.3478548451374538])
    n_quadrature = 2
    
    for i in range(n_quadrature):
        for j in range(n_quadrature):
            xi = qx[i]
            eta = qx[j]
            
            fx = pos[0,0] * 0.25 * (1.0 - xi) * (1.0 - eta)
            fx = fx + pos[1,0] * 0.25 * (1.0 + xi) * (1.0 - eta)
            fx = fx + pos[2,0] * 0.25 * (1.0 + xi) * (1.0 + eta)
            fx = fx + pos[3,0] * 0.25 * (1.0 - xi) * (1.0 + eta)
            fy = pos[0,1] * 0.25 * (1.0 - xi) * (1.0 - eta)
            fy = fy + pos[1,1] * 0.25 * (1.0 + xi) * (1.0 - eta)
            fy = fy + pos[2,1] * 0.25 * (1.0 + xi) * (1.0 + eta)
            fy = fy + pos[3,1] * 0.25 * (1.0 - xi) * (1.0 + eta);
            
            fxy = 5 * fx * fy # Poisson rhs
            
            d_N1_d_xi  = -0.25 * (1.0 - eta)
            d_N1_d_eta = -0.25 * (1.0 - xi)
            d_N2_d_xi  =  0.25 * (1.0 - eta)
            d_N2_d_eta = -0.25 * (1.0 + xi)
            d_N3_d_xi  =  0.25 * (1.0 + eta)
            d_N3_d_eta =  0.25 * (1.0 + xi)
            d_N4_d_xi  = -0.25 * (1.0 + eta)
            d_N4_d_eta =  0.25 * (1.0 - xi)
            
            d_x_d_xi  = d_N1_d_xi  * pos[0,0] + d_N2_d_xi  * pos[1,0] + d_N3_d_xi  * pos[2,0] + d_N4_d_xi  * pos[3,0]
            d_x_d_eta = d_N1_d_eta * pos[0,0] + d_N2_d_eta * pos[1,0] + d_N3_d_eta * pos[2,0] + d_N4_d_eta * pos[3,0]
            d_y_d_xi  = d_N1_d_xi  * pos[0,1] + d_N2_d_xi  * pos[1,1] + d_N3_d_xi  * pos[2,1] + d_N4_d_xi  * pos[3,1]
            d_y_d_eta = d_N1_d_eta * pos[0,1] + d_N2_d_eta * pos[1,1] + d_N3_d_eta * pos[2,1] + d_N4_d_eta * pos[3,1]
            
            
            jac = d_x_d_xi * d_y_d_eta - d_y_d_xi * d_x_d_eta
            
            N = np.array([(1 - xi)*(1 - eta)/4,
                          (1 + xi)*(1 - eta)/4,
                          (1 + xi)*(1 + eta)/4,
                          (1 - xi)*(1 + eta)/4])
            
            btemp[:] = btemp[:] + jac * fxy * N[:] * w[i] * w[j]
    btemp2 = np.zeros((4))
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
                
    if(((pos[2,0] == 0) & (pos[3,0] == 0)) |
       ((pos[2,0] == 1) & (pos[3,0] == 1)) | 
       ((pos[2,1] == 0) & (pos[3,1] == 0)) | 
       ((pos[2,1] == 1) & (pos[3,1] == 1))):
        btemp2 = btemp2 + fembeta(pos,2,3)
        
    if(((pos[3,0] == 0) & (pos[0,0] == 0)) |
       ((pos[3,0] == 1) & (pos[0,0] == 1)) | 
       ((pos[3,1] == 0) & (pos[0,1] == 0)) | 
       ((pos[3,1] == 1) & (pos[0,1] == 1))):
        btemp2 = btemp2 + fembeta(pos,3,0)
            
    btemp = btemp + btemp2
    test = 1
    
    for i in range(4):
        Bmat[ele[i]] += btemp[i] 