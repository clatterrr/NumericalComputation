import numpy as np
# CubbyFlow ENO
def upwind1(d,dx,direction):
    if direction == -1:
        return (d[0] - d[1])/dx
    elif direction == 1:
        return (d[2] - d[0])/dx
    
def centreDifference(d,dx):
    return (d[2] - d[0])/(2*dx)

def ENO3(d,dx):
    d1 = [d[1] - d[0],d[2] - d[1],d[3] - d[2],
          d[4] - d[3],d[5] - d[4],d[6] - d[5]]/dx
    d2 = [d1[1] - d1[0],d1[2] - d1[1],d1[3] - d1[2],
          d1[4] - d1[3],d1[5] - d1[4]] / dx / 2
    dfx = np.zeros((2))
    for k in range(2):
        d3 = np.zeros((2))
        if abs(d2[k+1]) < abs(d2[k+2]):
            c = d2[k+1]
            kstar = k - 1
            d3[0] = (d2[k+1] - d2[k]) / 3 / dx
            d3[1] = (d2[k+2] - d2[k+1]) / 3 / dx
        else:
            c = d2[k+2]
            kstar = k
            d3[0] = (d2[k+2] - d2[k+1]) / 3 / dx
            d3[1] = (d2[k+3] - d2[k+2]) / 3 / dx
            
        if abs(d3[0]) < abs(d3[1]):
            cstar = d3[0]
        else:
            cstar = d3[1]
        
        dq1 = d1[k+2]
        dq2 = c * (2 * (1 - k) - 1) * dx
        dq3 = cstar*(3*np.sqrt(1-k)-6*(1-kstar)+2)*dx*dx
        dfx[k] = dq1 + dq2 + dq3
    
    return dfx
        
def WENO5(v):
    vdev = np.zeros((5))
    dfx = np.zeros((2))
    dx = 1
    for k in range(2):
        if k == 0:
            for m in range(5):
                vdev[m] = (v[m+1] - v[m]) / dx
        else:
            for m in range(5):
                vdev[m] = (v[6-m] - v[5-m]) / dx
        phix1 = vdev[0]/3 - vdev[1]*7/6 + vdev[2]*11/6
        phix2 = -vdev[1]/6 + vdev[2]*5/6 + vdev[3]/3
        phix3 = vdev[2]/3 + vdev[3]*5/6 - vdev[4]/6
        
        s1 = (vdev[0]-2*vdev[1]+vdev[2])**2*13/12 + (vdev[0]-4*vdev[1]+3*vdev[2])**2/4
        s2 = (vdev[1]-2*vdev[2]+vdev[3])**2*13/12 + (vdev[1]-vdev[3])**2/4
        s3 = (vdev[2]-2*vdev[3]+vdev[4])**2*13/12 + (3*vdev[2]-4*vdev[3]+vdev[4])**2/4
        
        eps = 1e-10
        alpha1 = (s1 + eps)/10
        alpha2 = (s2 + eps)*3/5
        alpha3 = (s3 + eps)*3/10
        
        dfx[k] = (alpha1 * phix1 + alpha2 * phix2 + alpha3 * phix3) / (alpha1 + alpha2 + alpha3)
    
    return dfx
    
                
