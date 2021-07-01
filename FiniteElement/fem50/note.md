状态：完成

难易程度：简单

以下是针对pdf中不懂的地方的笔记

首先，moment reflection如下
$$
\eta_j(x,y) = det\begin{bmatrix} 1 & x & y \\ 1 & x_{j+1} & y_{j+1} \\ 1 & x_{j+1} & y_{j+2}\end{bmatrix} / det\begin{bmatrix} 1 & x_j & y_j \\ 1 & x_{j+1} & y_{j+1} \\ 1 & x_{j+1} & y_{j+2}\end{bmatrix}
$$


但是
$$
G = \begin{bmatrix} 1 & 1 & 1 \\ x_1 & x_2 & x_3 \\ y_1 & y_2 & y_3\end{bmatrix}^{-1}\begin{bmatrix} 0 & 0 \\ 1 & 0 \\ 0 & 1\end{bmatrix} = \\\begin{bmatrix} x_2y_3 - x_3y_2 & y_3 - y_2 & x_3-x_2 \\ x_1y_3 - x_3y_1 & y_3 - y_1 & x_3 - x_1 \\ x_1y_2 - x_2y_1 & y_2 - y_1 & x_2 - x_1\end{bmatrix}\begin{bmatrix} 0 & 0 \\ 1 & 0 \\ 0 & 1\end{bmatrix} / det = \begin{bmatrix} y_3 - y_2 & x_3 - x_2 \\ y_3 - y_1 & x_3 - x_1 \\ y_2 - y_1 & x_2 - x_1\end{bmatrix}
$$


```
G0 = np.array([[1,1,1],[x0,x1,x2],[y0,y1,y2]])
```



