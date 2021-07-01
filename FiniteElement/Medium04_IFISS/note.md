https://personalpages.manchester.ac.uk/staff/david.silvester/videolectures.html

配合上面的视频教程
$$
\begin{bmatrix} A & B^T \\ B & 0\end{bmatrix}\begin{bmatrix} u \\ p\end{bmatrix} = \begin{bmatrix} f \\ 0\end{bmatrix} \\
\begin{bmatrix} A & B^T \\0 & - BA^{-1}B^T\end{bmatrix}\begin{bmatrix} u \\ p\end{bmatrix} = \begin{bmatrix} f \\ BA^{-1}f\end{bmatrix}
$$
![image-20210628130514228](D:\图形学书籍\系列流体文章\gif\image-20210628130514228.png)

![image-20210628130801700](D:\图形学书籍\系列流体文章\gif\image-20210628130801700.png)

Unstable modes of the Q1-P0 element Griffiths, David and Silvester, David  
$$
\begin{bmatrix} A & 0 & B_x^t \\ 0 & A & B_y^t \\B_x & B_y & 0\end{bmatrix}\begin{bmatrix} U \\ V \\ P\end{bmatrix} = \begin{bmatrix} F_x \\ F_y \\ G\end{bmatrix}
$$
又
$$
BK^{-1}B^t\bold P = \sigma M \bold P \qquad B^tM^{-1}B\begin{bmatrix}\bold U \\ \bold V \end{bmatrix} = \sigma K \begin{bmatrix} \bold U \\ \bold V \end{bmatrix}
$$

$$
\begin{bmatrix} A & 0 & B_x^t \\ 0 & A & B_y^t \\B_x & B_y & 0\end{bmatrix}\begin{bmatrix} U \\ V \\ P\end{bmatrix} = \lambda \begin{bmatrix} A & 0 & 0 \\ 0 & A & 0 \\ 0 & 0 &M\end{bmatrix}\begin{bmatrix} U \\ V \\ P\end{bmatrix}
$$

Finite Element and Fast Solver
$$
\int_\Omega (\nabla ^2u + f)v = 0 \tag{1.8}
$$
那么
$$
-\int_\Omega v\nabla^2u = \int_\Omega \nabla u \cdot \nabla v - \int_\Omega \nabla \cdot (v \nabla u) = \int_\Omega \nabla u \cdot \nabla v - \int_{\partial \Omega}v\frac{\partial u}{\partial n}
$$
依照neumann边界可以写成
$$
\int_\Omega \nabla u \cdot \nabla v = \int_\Omega vf + \int_{\partial \Omega}vg
$$
![image-20210628135228300](D:\图形学书籍\系列流体文章\gif\image-20210628135228300.png)

![image-20210628135510361](D:\图形学书籍\系列流体文章\gif\image-20210628135510361.png)

![image-20210628135524125](D:\图形学书籍\系列流体文章\gif\image-20210628135524125.png)