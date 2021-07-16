import numpy as np
import math
import random
volsize = [4,2,2]
volres = [64,32,32]
Nx = 64
Ny = 32
Nz = 32
hbar = 0.1 # 普朗克常量
dt = 1 / 48
tmax = 50

jet_velocity = [1,0,0]
dx = volsize[0] / volres[0]
dy = volsize[1] / volres[1]
dz = volsize[2] / volres[2]


norm = np.sqrt(1*1 + 0.01*0.01)
psi1 = np.ones((volres[0],volres[1],volres[2]),dtype = complex)/norm
psi2 = np.ones((volres[0],volres[1],volres[2]),dtype = complex)/100/norm

mask = np.zeros((volres[0],volres[1],volres[2]),dtype = complex)

def Prepare(psi1,psi2,phase):
    for i in range(64):
        for j in range(32):
            for k in range(32):
                if ((((j*dy - nozzle_cen[1])**2 + (k*dz - nozzle_cen[2])**2) <= nozzle_rad**2) and
                         (abs(i*dx - nozzle_cen[0] ) <= nozzle_len / 2 )):
                    amp1 = abs(psi1[i,j,k])
                    amp2 = abs(psi2[i,j,k])
                    psi1[i,j,k] = amp1 * np.exp(1j * phase[i,j,k])
                    psi2[i,j,k] = amp2 * np.exp(1j * phase[i,j,k])
    return psi1,psi2

def VelocityOneForm(psi1,psi2):
    vx = np.zeros((64,32,32))
    vy = np.zeros((64,32,32))
    vz = np.zeros((64,32,32))
    for i in range(64):
        for j in range(32):
            for k in range(32):
                    ixp = int((i+1)%64)
                    jyp = int((j+1)%32)
                    kzp = int((k + 1) % 32)
                    term = np.conj(psi1[i,j,k])*psi1[ixp,j,k] + np.conj(psi2[i,j,k])*psi2[ixp,j,k]
                    vx[i,j,k] = math.atan2(term.imag,term.real) 
                    term = np.conj(psi1[i,j,k])*psi1[i,jyp,k] + np.conj(psi2[i,j,k])*psi2[i,jyp,k]
                    vy[i,j,k] = math.atan2(term.imag,term.real) 
                    term = np.conj(psi1[i,j,k])*psi1[i,j,kzp] + np.conj(psi2[i,j,k])*psi2[i,j,kzp]
                    vz[i,j,k] = math.atan2(term.imag,term.real)
    return vx,vy,vz

def PressureProject(psi1,psi2):
    vx,vy,vz = VelocityOneForm(psi1,psi2)
    div = np.zeros((64,32,32))
    for i in range(64):
        for j in range(32):
            for k in range(32):
                ixm = int((i-1 + 64) % 64)
                jym = int((j-1 + 32) % 32)
                kzm = int((k-1 + 32) % 32)
                div[i,j,k] = (vx[i,j,k] - vx[ixm,j,k])/dx/dx + (vy[i,j,k]
                                - vy[i,jym,k])/dy/dy + (vz[i,j,k] - vz[i,j,kzm])/dz/dz
    
    # 解算压力泊松方程
    # fftn 算得有问题，numpy算的和matlab的精度不同
    f = np.fft.fftn(div)
    fac = np.zeros((Nx,Ny,Nz))
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                sx = np.sin(np.pi * i / Nx) / dx
                sy = np.sin(np.pi * j / Ny) / dy
                sz = np.sin(np.pi * k / Nz) / dz
                denom = sx**2 + sy**2 + sz**2
                if denom == 0:
                    fac[i,j,k] = 0
                else:
                    fac[i,j,k] = - 1 / denom / 4
    f = f * fac
    p = np.fft.ifftn(f)
    
    # GaugeTransform
    eip = np.exp(-1j * p)
    psi1 = psi1 * eip
    psi2 = psi2 * eip
    return psi1,psi2

def VelocityRK(vx,vy,vz,newx,newy,newz):
    ux = np.zeros((n_particles))
    uy = np.zeros((n_particles))
    uz = np.zeros((n_particles))
    for i in range(n_particles):
        px = newx[i] % 4
        py = newy[i] % 2
        pz = newz[i] % 2

        ix = int((px / dx + 1) % Nx)
        iy = int((py / dy + 1) % Ny)
        iz = int((pz / dz + 1) % Nz)
        
        ixp = int((ix + 1) % Nx)
        iyp = int((iy + 1) % Ny)
        izp = int((iz + 1) % Nz)
        
        wx = (px / dx + 1) % Nx - ix
        wy = (py / dy + 1) % Ny - iy
        wz = (pz / dz + 1) % Nz - iz
        
        ux[i] = (1 - wz)*((1 - wy)*vx[ix,iy,iz] + wy * vx[ix,iyp,iz])
        ux[i] += wz*((1 - wy)*vx[ix,iy,izp] + wy * vx[ix,iyp,izp])
        
        uy[i] = (1 - wz)*((1 - wx)*vy[ix,iy,iz] + wx * vy[ixp,iy,iz])
        uy[i] += wz*((1 - wx)*vy[ix,iy,izp] + wx * vy[ixp,iy,izp])
        
        uz[i] = (1 - wy)*((1 - wx)*vz[ix,iy,iz] + wx * vx[ixp,iy,iz])
        uz[i] += wy*((1 - wx)*vz[ix,iyp,iz] + wx * vx[ixp,iyp,iz])
    return ux,uy,uz
    

vec = 10
for i in range(64):    
    for j in range(32):
        for k in range(volres[2]):
            kx = (i - volres[0] / 2) / volsize[0]
            ky = (j - volres[1] / 2) / volsize[1]
            kz = (k - volres[2] / 2) / volsize[2]
            fac = -4 * np.pi**2 * hbar
            lam = fac *(kx**2 + ky**2 + kz**2)
            mask[i,j,k] = np.exp(1j * lam * dt / 2)
        
    
nozzle_cen = [2-1.7,1-0.034,1+0.066]
nozzle_len = 0.5
nozzle_rad = 0.5

n_particles = 50
particle_x = np.zeros((n_particles))
particle_y = np.zeros((n_particles))
particle_z = np.zeros((n_particles))

for ite in range(10):
    phase = np.zeros((Nx,Ny,Nz))
    for i in range(Nx):
        phase[i,:,:] = dx * i * vec
    psi1,psi2 = Prepare(psi1, psi2,phase)
    psi1,psi2 = PressureProject(psi1, psi2)
               
iter_max = 1
for ite in range(iter_max):
    
    # SchroedingerFlow
    p1 = np.fft.fftshift(np.fft.fftn(psi1))
    p2 = np.fft.fftshift(np.fft.fftn(psi2))
    p1 = p1 * mask
    p2 = p2 * mask
    psi1 = np.fft.ifftn(np.fft.fftshift(p1))
    psi2 = np.fft.ifftn(np.fft.fftshift(p2))
    
    norm = np.sqrt(abs(psi1)**2 + abs(psi2)**2)
    psi1 /= norm
    psi2 /= norm
    
    psi1,psi2 = PressureProject(psi1, psi2)
    
    omega = 5
    phase = np.zeros((Nx,Ny,Nz))
    for i in range(Nx):
        phase[i,:,:] = i * vec * dx - omega * ite * dt
    psi1,psi2 = Prepare(psi1, psi2,phase)
    psi1,psi2 = PressureProject(psi1,psi2)
    
    newx = np.zeros((n_particles))
    newy = np.zeros((n_particles))
    newz = np.zeros((n_particles))
    for i in range(n_particles):
        rt = random.random()
        newx[i] = nozzle_cen[0]
        newy[i] = nozzle_cen[1] + 0.9 * nozzle_rad * np.cos(rt)
        newz[i] = nozzle_cen[2] + 0.9 * nozzle_rad * np.sin(rt)
        particle_x[i] = newx[i]
        particle_y[i] = newy[i]
        particle_z[i] = newz[i]
        
        
        
    vx,vy,vz = VelocityOneForm(psi1, psi2)
    vx = vx / dx / 10
    vy = vy / dy / 10
    vz = vz / dz / 10
    
    ux1,uy1,uz1 = VelocityRK(vx,vy,vz,particle_x,particle_y,particle_z)
    ux2,uy2,uz2 = VelocityRK(vx,vy,vz,particle_x+ux1*dt/2, particle_y + uy1*dt/2, particle_z+uz1*dt/2)
    ux3,uy3,uz3 = VelocityRK(vx,vy,vz,particle_x+ux2*dt/2, particle_y + uy2*dt/2, particle_z+uz2*dt/2)
    ux4,uy4,uz4 = VelocityRK(vx,vy,vz,particle_x+ux3*dt, particle_y + uy3*dt, particle_z+uz3*dt)
        
    particle_x_new = particle_x + dt/6*(ux1 + 2*ux2 + 2*ux3 + ux4)
    particle_y_new = particle_y + dt/6*(uy1 + 2*uy2 + 2*uy3 + uy4)
    particle_z_new = particle_z + dt/6*(uz1 + 2*uz2 + 2*uz3 + uz4)
    
    # particle_x = particle_x_new.copy()
    # particle_y = particle_y_new.copy()
    # particle_z = particle_z_new.copy()
    