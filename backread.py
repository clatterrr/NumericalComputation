import csv
import numpy as np

# 读取csv至字典
csvFile = open("back.csv", "r")
reader = csv.reader(csvFile)

# 建立空字典
result = {}
idx = 0
nx = 125
ny = 50
size = nx * ny
u = np.zeros((nx,ny))
v = np.zeros((nx,ny))
rho = np.zeros((nx,ny))
ru = np.zeros((nx,ny))
rv = np.zeros((nx,ny))
for item in reader:
    # 忽略第一行
    if reader.line_num == 1:
        continue
    if idx == 5250:
        break
    if idx < 1250:
        i = int(idx % 25)
        j = int(idx // 25)
    else:
        i = int((i - idx) % 100)
        j = int((idx - 1250) / 100)
    u[i,j] = item[6]
    v[i,j] = item[7]
    rho[i,j] = item[10]
    ru[i,j] = rho[i,j] * u[i,j]
    rv[i,j] = rho[i,j] * v[i,j]
    idx += 1
csvFile.close()

time = 0
timeFinal = 1
dx = 4.1 / 60
dy = 1 / 30
drho_dx = np.zeros((nx,ny))
drho_dy = np.zeros((nx,ny))
vdx = np.zeros((nx,ny))
vdy = np.zeros((nx,ny))
def p(a0,a1,a2,d):
    al = (a1 - a0)/d
    ar = (a2 - a1)/d
    small=1e-12
    aa=np.sqrt(small+ar*ar)
    bb=np.sqrt(small+al*al)
    return ar*al*(np.sign(ar)+np.sign(al))/(aa+bb)
while(time < timeFinal):
    
    for i in range(nx-1):
        for j in range(ny):
            drho_dx[i,j] = (ru[i+1,j] - ru[i,j])/dx
    
    for i in range(nx):
        for j in range(ny-1):
            drho_dy[i,j] = (rv[i,j+1] - rv[i,j])/dy
            
    for i in range(1,nx-1):
        for j in range(1,ny-1):
            vdx[i,j] = p(ru[i-1,j],ru[i,j],ru[i+1,j],dx)
            vdy[i,j] = p(rv[i,j-1],rv[i,j],rv[i,j+1],dy)
            
            
    time = 10
            
    