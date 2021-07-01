import numpy as np
import scipy.io as scio

# 本地地址：D:\FluidSim\FluidSim\FEMNEW\Navier-stokes-Satsit14master



tmax = 100
celem = scio.loadmat('celem.mat')['celem'] - 1
node = scio.loadmat('node.mat')['node']
nmax = 341
u0 = np.ones((nmax))
v0 = np.zeros((nmax))
    
numOfElements = 100
for t in range(0,tmax):
    
    # Assembly matrix
    
    for e in range(numOfElements):
        xp = np.array([-np.sqrt(3/5),0,np.sqrt(3/5)])
        wp = np.array([5/9,8/9,5/9])
        w = np.zeros((9))
        gp = np.zeros((9,2))
        ngp = 9
        for i in range(3):
            for j in range(3):
                idx = int(3 * i + j)
                w[idx] = wp[i] * wp[j]
                gp[idx,0] = xp[j]
                gp[idx,1] = xp[i]
        
        N = np.zeros((8,9))
        dN = np.zeros((2,8,9))
        M = np.zeros((4,9))
        dM = np.zeros((2,4,9))
        A1 = np.zeros((8,8))
        A2 = np.zeros((8,4))
        A8 = np.zeros((8,4))
        A4 = np.zeros((4,8))
        A6 = np.zeros((4,8))
        
        
        for k in range(ngp):
            xi = gp[k,0]
            eta = gp[k,1]
            N[0,k] = -(1 - xi) * (1 - eta) * (1 + xi + eta) / 4
            N[1,k] = (1 - 2*xi)*(1 - eta) / 2
            N[2,k] = -(1 + xi)*(1 - eta)*(1 - xi + eta) / 4
            N[3,k] = (1 + xi)*(1 - eta**2) / 2
            N[4,k] = -(1 + xi)*(1 + eta)*(1 - xi - eta) / 4
            N[5,k] = (1 - xi**2)*(1 + eta)/2
            N[6,k] = -(1 - xi)*(1 + eta)*(1 + xi - eta)/4                
            N[7,k] = (1 - xi)*(1 - eta**2)/2
            
            dN[0,0,k] = (1 - eta) * (2 * xi + eta) / 4
            dN[0,1,k] = - xi * (1 - eta)
            dN[0,2,k] = (eta - 1)*(eta - 2*xi)/4
            dN[0,3,k] = (1 - eta**2)/2
            dN[0,4,k] = (eta + 1)*(eta + 2*xi)/4
            dN[0,5,k] = - xi * (eta + 1)
            dN[0,6,k] = -(1 + eta)*(eta - 2*xi)/4
            dN[0,7,k] = (-1 + eta**2) / 2
            
            dN[1,0,k] = (1 - xi)*(2 * eta + xi)/4
            dN[1,1,k] = (-1 + xi**2)/2
            dN[1,2,k] = (1 + xi)*(2*eta - xi)/4
            dN[1,3,k] = -eta*(1 + xi)
            dN[1,4,k] = (1 + xi)*(xi + 2*eta)/4
            dN[1,5,k] = (1 - xi**2)/2
            dN[1,6,k] = (1 - xi)*(2*eta - xi)/4
            dN[1,7,k] = eta * (xi - 1)
            
            M[0,k] = (1 - xi)*(1 - eta)/4
            M[1,k] = (1 + xi)*(1 - eta)/4
            M[2,k] = (1 + xi)*(1 + eta)/4
            M[3,k] = (1 - xi)*(1 + eta)/4
            
            dM[0,0,k] = - (1 - eta)/4
            dM[0,1,k] = (1 - eta)/4
            dM[0,2,k] = (1 + eta)/4
            dM[0,3,k] = - (1 + eta)/4
            
            dM[1,0,k] = - (1 + eta)/4
            dM[1,1,k] = -(1 + xi)/4
            dM[1,2,k] = (1 + xi)/4
            dM[1,3,k] = (1 - xi)/4
            
        
                    
        
        for k in range(ngp):
            ubar = 0
            vbar = 0
            for i in range(8):
                ubar = ubar + N[i,k] * u0[i]
                vbar = vbar + N[i,k] * v0[i]
            coord8 = np.zeros((8,2))
            for i in range(8):
                coord8[i,:] = node[celem[e,i],:]
            coord4 = np.zeros((4,2))
            for i in range(4):
                coord4[i,:] = node[celem[e,i+8],:]
            Jacob8 = np.dot(dN[:,:,k],coord8[:,:])
            detJac = abs(np.linalg.det(Jacob8))
            Jinv8 = np.linalg.inv(Jacob8)
            gdN = np.zeros((2,8))
            for i in range(2):
                for j in range(2):
                    gdN[i,:] += Jinv8[i,j]*dN[j,:,k]
            Re = 1
            
            for i in range(8):
                for j in range(8):
                    A1[i,j] = A1[i,j] + w[k] * (N[i,k] * ubar * gdN[0,j] + N[i,k] * vbar * gdN[1,j] + 
                                 (gdN[0,i] * gdN[0,j] + gdN[1,i] * gdN[1,j])/Re) * detJac
                    
            A9 = A1.copy()
            Jacob4 = np.dot(dM[:,:,k],coord4[:,:])
            Jinv4 = np.linalg.inv(Jacob8)
            gdM = np.zeros((2,4))
            for i in range(2):
                for j in range(2):
                    gdM[i,:] += Jinv4[i,j]*dM[j,:,k]
            detJac4 = abs(np.linalg.det(Jacob4))
            
            for i in range(8):
                for j in range(4):
                    A2[i,j] = A2[i,j] + w[k]*N[i,k]*gdM[0,j]*detJac4
                    A8[i,j] = A8[i,j] + w[k]*N[i,k]*gdM[1,j]*detJac4
            for i in range(4):
                for j in range(8):
                    A4[i,j] = A4[i,j] + w[k]*M[i,k]*gdN[0,j]*detJac4
                    A6[i,j] = A6[i,j] + w[k]*M[i,k]*gdN[1,j]*detJac4
                    
        tes = 1
                
            
                
        