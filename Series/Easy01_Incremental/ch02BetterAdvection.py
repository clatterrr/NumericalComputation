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

def cubicInterpolate(x,y,value):
    x = min(max(x - 0.5,0),sizeX - 1.001)
    y = min(max(y - 0.5,0),sizeY - 1.001)
    ix = int(x)
    iy = int(y)
    x = x - ix
    y = y - iy
    
    x0 = max(ix - 1,0)
    x1 = ix
    x2 = ix + 1
    x3 = min(ix + 2,sizeX - 1)
    
    y0 = max(iy - 1,0)
    y1 = iy
    y2 = iy + 1
    y3 = min(iy + 2,sizeY - 1)
    
    q0 = cerp(value[x0,y0],value[x1,y0],value[x2,y0],value[x3,y0],x)
    q1 = cerp(value[x0,y1],value[x1,y1],value[x2,y1],value[x3,y1],x)
    q2 = cerp(value[x0,y2],value[x1,y2],value[x2,y2],value[x3,y2],x)
    q3 = cerp(value[x0,y3],value[x1,y3],value[x2,y3],value[x3,y3],x)
    
    return cerp(q0,q1,q2,q3,y)
    
def cerp(a,b,c,d,x):
    xsq = x * x
    xcu = xsq * x
    minV = min(a,min(b,min(c,d)))
    maxV = max(a,max(b,max(c,d)))
    t0 = a * (0 - 0.5*x + 1 * xsq - 0.5 * xcu)
    t1 = b * (1 + 0*x - 2.5 * xsq + 1.5 * xcu)
    t2 = c * (0 + 0.5*x + 2 * xsq - 1.5 * xcu)
    t3 = d * (0 + 0*x - 0.5 * xsq + 0.5 * xcu)
    return min(max(t0+t1+t2+t3,minV),maxV)
    
def advect(value):
    for j in range(sizeY):
        for i in range(sizeX):
            x = i + 0.5 # 网格中心点
            y = j + 0.5 # 网格中心点
            
            # RungeKutta 3 order
            firstU = lerp(x, y, _u) / hx
            firstV = lerp(x, y, _v) / hx
            
            midX = x - 0.5 * timestep * firstU
            midY = y - 0.5 * timestep * firstV
            
            midU = lerp(x, y, _u) / hx
            midV = lerp(x, y, _v) / hx
            
            lastX = x - 0.75 * timestep * firstU
            lastY = x - 0.75 * timestep * firstV
            
            lastU = lerp(x, y, _u) / hx
            lastV = lerp(x, y, _v) / hx
            
            x -= timestep * (firstU*2/9 + midU*3/9 + lastU*4/9)
            y -= timestep * (firstU*2/9 + midU*3/9 + lastU*4/9)
            
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
    
    

            
            