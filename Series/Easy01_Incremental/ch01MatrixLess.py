import numpy as np
sizeX = 128
sizeY = 128
hx = 1 / min(sizeX,sizeY)
density = 0.1
timestep = 0.005

_d = np.zeros((sizeX,sizeY)) # 密度
_u = np.zeros((sizeX,sizeY)) # x轴速度
_v = np.zeros((sizeX,sizeY)) # y轴速度
_p = np.zeros((sizeX,sizeY)) # 压力
_r = np.zeros((sizeX,sizeY)) # right hand side

time = 0
timeFinal = 1

def lerp(x,y,value):
    x = min(max(x - 0.5,0),sizeX - 1.001)
    y = min(max(y - 0.5,0),sizeY - 1.001)
    ix = int(x)
    iy = int(y)
    x = x - ix
    y = y - iy
    
    v00 = value[ix,iy]
    v10 = value[ix+1,iy]
    v01 = value[ix,iy+1]
    v11 = value[ix+1,iy+1]
    
    v0 = v00 * (1 - x) + v10 * x
    v1 = v01 * (1 - x) + v11 * x
    v2 = v0 * (1 - y) + v1 * y
    return v2
    
def advect(value):
    for j in range(sizeY):
        for i in range(sizeX):
            x = i + 0.5 # 网格中心点
            y = j + 0.5 # 网格中心点
            
            uvel = lerp(x, y, _u) / hx
            vvel = lerp(x, y, _v) / hx
            
            x -= uvel * timestep
            y -= vvel * timestep
            
            value[i,j] = lerp(x, y, value)
    

while(time < timeFinal):
    time = 10
    
    # 初始化固定值
    for i in range(sizeX):
        for j in range(sizeY):
            if ((i * hx >= 0.45) & (i * hx <= 0.45 + 0.1)
                & (j * hx >= 0.2) & (j * hx <= 0.2 + 0.01)):
                _d[i,j] = 1
                _u[i,j] = 0
                _v[i,j] = 3
                
    # 计算右手项
    scale = 1 / hx
    for j in range(sizeY):
        for i in range(sizeX):
            _r[i,j] =scale * (_u[i,j] - _u[i+1,j]  + _v[i,j] - _v[i,j+1])
            
    # 解算压力泊松方程
    scale = timestep / (density * hx * hx)
    for it in range(100):
        maxDelta = 0
        for j in range(sizeY):
            for i in range(sizeX):
                diag = 0
                offDiag = 0
                if i > 0:
                    diag += scale
                    offDiag -= scale * _p[i-1,j]
                if j > 0:
                    diag += scale
                    offDiag -= scale * _p[i,j-1]
                if i < sizeX - 1:
                    diag += scale
                    offDiag -= scale * _p[i+1,j]
                if j < sizeY - 1:
                    diag += scale
                    offDiag -= scale * _p[i,j+1]
                newP = (_r[i,j] - offDiag) / diag
                maxDelta = max(maxDelta,abs(_p[i,j] - newP))
                _p[i,j] = newP
        if maxDelta < 1e-5:
            break
        
    # 更新速度
    scale = timestep / (density * hx)
    for j in range(sizeY):
        for i in range(sizeX):
            _u[i,j] -= scale * _p[i,j]
            _u[i+1,j] += scale * _p[i,j]
            _v[i,j] -= scale * _p[i,j]
            _v[i,j+1] += scale * _p[i,j]
            
    for i in range(sizeX):
        _v[i,0] = _v[i,sizeY - 1] = 0
    for j in range(sizeY):
        _u[0,i] = _u[sizeX-1,i] = 0
    
    # 对流,应该用备份的，下面的代码是错的
    advect(_d)
    advect(_u)
    advect(_v)
    
    

            
            