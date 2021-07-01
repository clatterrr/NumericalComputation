本篇对应的代码为

本篇笔记有对应的论文，为"M. Cermak, S. Sysala, J. Valdman: Efficient and flexible MATLAB implementation of 2D and 3D elastoplastic problems. Applied Mathematics and Computation 355, (2019) pp. 595-614, DOI: 10.1016/j.amc.2019.02.054" available also at https://arxiv.org/pdf/1805.04155.pdf

当应力超过yield stress，材料就开始yields，或者说是从弹性区进入塑性区。应力和分为volumetric (or hydrostatic)和deviatoric 部分，也就是正应力和斜应力部分。

当材料做功的时候，能量被作为应变能量(Strain Energy)存储在体内。应变能量密度为
$$
W = \frac{1}{2}\sigma \varepsilon
$$
这部分能量，一部分是正应力能量，用于改变物体的体积，另一部分用于改变的物体的形状。而Von Mises和后者有关，也就是
$$
W = \frac{1}{2}\sigma : \varepsilon = \frac{1}{2}(\sigma_v + \sigma_d):(\varepsilon_v + \varepsilon_d) \\
= \frac{1}{2}(\sigma_v : \varepsilon_v + \sigma_v : \varepsilon_d + \sigma_d : \varepsilon_v + \sigma_d : \varepsilon_d)
$$
由于正张量和斜张量的乘积总是零，因此应变能量变为如下式，被分为volumetric 和 deviatoric 部分
$$
W = \frac{1}{2} \sigma_v : \varepsilon_v + \frac{1}{2}\sigma_d : \varepsilon_d
$$
现在我们可以重写deviatoric 能量如下
$$
\large
\varepsilon_d = \frac{1}{2G}\sigma_d \Rightarrow W_{dev} = \frac{1}{4G}\sigma_d : \sigma_d = \frac{1}{4G}\sigma_{rep}^2 \Rightarrow \sigma_{rep} = \sqrt{\sigma_d : \sigma_d}
$$
这就是von mises 应力。1D状态如下
$$
\sigma = \begin{bmatrix}\sigma & 0 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & 0 \end{bmatrix} \Rightarrow \sigma_d = \begin{bmatrix}2\sigma/3 & 0 & 0 \\ 0 & -\sigma/3 & 0 \\ 0 & 0 & -\sigma/3 \end{bmatrix} \Rightarrow \sigma_{rep} = \sqrt{2/3}\sigma
$$
 

有限元一般代求解的问题如下，其中sigma是应力，u是位移，K是刚度矩阵。
$$
\sigma = K u
$$
然而K不能直接求出来，我们现在所知道的只有应力与应变之间的关系
$$
\{\sigma\} = \{\sigma_x \space \sigma_y \space \sigma_z\space \tau_{xy} \space \tau_{yz}\space \sigma_{zx}\} \qquad \{\sigma\} = [E]\{\varepsilon\}
$$
其中E是弹性矩阵
$$
[E] = \begin{bmatrix} \lambda + 2\mu & \lambda & \lambda & 0 & 0 & 0 \\ \lambda & \lambda + 2\mu & \lambda & 0 & 0 & 0 \\ \lambda & \lambda & \lambda + 2\mu & 0 & 0 & 0  \\ 0& 0 & 0 & \mu & 0 & 0 \\ 0& 0 & 0 & 0 & \mu & 0 \\ 0& 0 & 0 & 0 & 0 & \mu\end{bmatrix}
$$


```python
Vol = np.array([[1,1,0],[1,1,0],[0,0,0]]) # lambda
Dia = np.array([[1,0,0],[0,1,0],[0,0,1/2]]) # mu
Dev = Dia - Vol / 3
Elast = np.zeros((9,n_int))
D = np.zeros((450,450))
for i in range(150):
    for j in range(3):
        for k in range(3):
            idx = int(j * 3 + k)
            Elast[idx,i] = 2 * Dev[j,k] * shearMat[i] + Vol[j,k] * bulkMat[i]
```



还有应变与位移之间的关系
$$
\{\varepsilon\} = [D] \{u\} \\
\large
\{\varepsilon\} = \{\varepsilon_x \space \varepsilon_y \space \varepsilon_z\space \gamma_{xy}\space \gamma_{yz}\space \gamma_{zx}\} \qquad \{u\} = \{u \space v \space w\}
$$

$$
[D] = \begin{bmatrix} \partial /\partial x & 0 & 0 \\ 0 & \partial /\partial y & 0 \\ 0 & 0 & \partial /\partial z \\ \partial /\partial y & \partial /\partial x & 0 \\ 0 & \partial /\partial z & \partial /\partial y \\ \partial /\partial z & 0 & \partial /\partial x\end{bmatrix}
$$

而且D矩阵不能直接乘上位移矩阵，而是乘上形函数
$$
[B] = \begin{bmatrix} \partial N/\partial x & 0 & 0 \\ 0 & \partial N/\partial y & 0 \\ 0 & 0 & \partial N/\partial z \\ \partial N/\partial y & \partial N/\partial x & 0 \\ 0 & \partial N/\partial z & \partial N/\partial y \\ \partial N/\partial z & 0 & \partial N/\partial x\end{bmatrix}
$$
二维的如下。有三行是因为有三个顶点，有两列是因为有坐标轴为2，也就是二维的
$$
B = \begin{bmatrix} \partial N_i / \partial x & 0 \\ 0 & \partial N_i / \partial y \\ \partial N_i /\partial y & \partial N_i /\partial  x\end{bmatrix}
$$

```python
# 96 个点
# 150 个元素
B = np.zeros((450,192))
for i in range(150):
    
    idyb = np.zeros((6))
    idyb[0] = elem[0,i]*2 - 2 # 顶点0的x轴位移编号
    idyb[1] = elem[0,i]*2 - 1 # 顶点0的y轴位移编号
    idyb[2] = elem[1,i]*2 - 2
    idyb[3] = elem[1,i]*2 - 1
    idyb[4] = elem[2,i]*2 - 2
    idyb[5] = elem[2,i]*2 - 1
    for j in range(18):
        
        idx = int(i * 3 + j % 3)
        idy = int(idyb[j // 3])
        B[idx,idy] = vB[j,i]
```

而dN/dx并不好求，所以要算dN/dxi和dN/deta
$$
\begin{bmatrix} \partial N_i/\partial \xi  \\ \partial N_i /\partial \eta\end{bmatrix} = \begin{bmatrix} \partial x/\partial \xi  & \partial y /\partial \xi\\ \partial x /\partial \eta & \partial y / \partial \eta\end{bmatrix}\begin{bmatrix} \partial N_i/\partial x  \\ \partial N_i /\partial y\end{bmatrix}
$$

那么反过来
$$
\begin{bmatrix} \partial N_i/\partial x  \\ \partial N_i /\partial y\end{bmatrix} = [J]^{-1}\begin{bmatrix} \partial N_i/\partial \xi  \\ \partial N_i /\partial \eta\end{bmatrix}
$$

```python
for i in range(n_int):
    Dphi1[:,i] = Jinv11[i] * DHatP1[:,0] + Jinv12[i] * DHatP2[:,0]
    Dphi2[:,i] = Jinv21[i] * DHatP1[:,0] + Jinv22[i] * DHatP2[:,0]
```

雅可比的求解稍微需要点技巧。比如有四个点(0,0)和(1,0)和(0,1)和(2,2)，它们的值分别如下
$$
U(0,0) = 0 \quad U(1,0) = 1 \quad U(0,1) = 1 \quad U(2,2) = 4
$$
选用的四个形函数如下
$$
N_1 = \frac{1}{4}(1-\xi)(1- \eta) \quad N_2 = \frac{1}{4}(1+\xi)(1 - \eta) \\ N_3 = \frac{1}{4}(1 -\xi)(1 + \eta) \quad N_4 = \frac{1}{4}(1+\xi)(1 + \eta)
$$


那么它们的导数如下
$$
\frac{\partial N_1}{\partial \xi} = \frac{1}{4}(-1 + \eta) \qquad \frac{\partial N_2}{\partial \xi} = \frac{1}{4}(1 - \eta) \\\frac{\partial N_1}{\partial \xi} = \frac{1}{4}(-1 - \eta) \qquad \frac{\partial N_2}{\partial \xi} = \frac{1}{4}(1 + \eta) \\ \frac{\partial N_1}{\partial \eta} = \frac{1}{4}(\xi - 1) \qquad \frac{\partial N_2}{\partial \eta} = \frac{1}{4}( - 1 - \xi) \\\frac{\partial N_1}{\partial \eta} = \frac{1}{4}(1 - \xi) \qquad \frac{\partial N_2}{\partial \eta} = \frac{1}{4}(1 + \xi)
$$

```matlab
  case 'Q1'
    % - the reference square with coordinates:
    %         [-1,-1], [1,-1], [1,1], [-1,1] 
    % - n_p_s=4, n_q_s=length(xi_1)    
    HatP=[(1-xi_1).*(1-xi_2)/4; (1+xi_1).*(1-xi_2)/4;
          (1+xi_1).*(1+xi_2)/4; (1-xi_1).*(1+xi_2)/4 ];
    DHatP1 = [-(1-xi_2)/4;  (1-xi_2)/4;
               (1+xi_2)/4; -(1+xi_2)/4 ];
    DHatP2 = [-(1-xi_1)/4; -(1+xi_1)/4;
                 (1+xi_1)/4;  (1-xi_1)/4 ];
```

对形函数求导之后的数并不是常数，那么这些导数中的xi和eta中究竟应该是什么值？注意使用形函数的最终目标是积分，因此可以直接使用高斯积分。一维二点高斯积分如下
$$
x_1 = -\frac{1}{\sqrt{3}} \quad x_2 = \frac{1}{\sqrt{3}} \quad w_1 = w_2 = 1
$$
二维二点高斯积分如下
$$

\xi = \begin{bmatrix} -1/\sqrt{3} & -1/\sqrt{3} & 1/\sqrt{3} & 1/\sqrt{3} \end{bmatrix} \\ \eta =  \begin{bmatrix}-1/\sqrt{3} & 1/\sqrt{3} & -1/\sqrt{3} & 1/\sqrt{3}\end{bmatrix}
$$
上式的意思是，假如某个函数为
$$
U = \xi^2 + 5\eta^2
$$
那么积分
$$
\int_{-1}^{1}d\xi \int_{-1}^{1} Ud\eta = ((-\frac{1}{\sqrt{3}})^2 + 5(-\frac{1}{\sqrt{3}})^2) + ((-\frac{1}{\sqrt{3}})^2 + 5(\frac{1}{\sqrt{3}})^2) \\ +((\frac{1}{\sqrt{3}})^2 + 5(-\frac{1}{\sqrt{3}})^2) + ((\frac{1}{\sqrt{3}})^2 + 5(\frac{1}{\sqrt{3}})^2) \\ = 8
$$




Energy Potenial

```matlab
e=(1/2)*U(:)'*K*U(:)-(f_t(:)'+f_V(:)')*U(:); 
```

$$
\sum \frac{1}{2}q^T Kq - \sum q^T f
$$


