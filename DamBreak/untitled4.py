import numpy as np

n = 8
b = np.zeros((n*n)) # Ax = b
x = np.zeros((n*n)) # 未知量
d = np.zeros((n*n)) # 方向direction
r = np.zeros((n*n)) # 残差residual

def domain(var,i,j):
    idx = j * n + i
    if i < 0 or j < 0 or i >= n  or j >= n:
        return 0
    return var[idx]

def A(var,i,j):
    idx = j * n + i
    return (var[idx]*4 - domain(i-1, j) 
            - domain(i+1, j)- domain(i, j-1)- domain(i, j+1))

b[8//2*8+8//2] = 1

for j in range(n):
    for i in range(n):
        idx = j * n + i
        d[idx] = b[idx] - A(x,i,j)
        r[idx] = d[idx]
    
    


