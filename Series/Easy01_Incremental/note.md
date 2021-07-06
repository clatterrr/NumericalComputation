




$$
v <- \frac{(r_{n-1},r_{n-1})}{(s,z)}
$$

```c++
double sigma = dotProduct(_z, _r);
        
        for (int iter = 0; iter < limit; iter++) {
            matrixVectorProduct(_z, _s);
            double alpha = sigma/dotProduct(_z, _s);
```

不过算alpha的时候，分子上式很奇怪？


$$
x_n <- x_{n-1} + v  *s\\
r_n <- r_{n-1} - v * z
$$

```c++
 scaledAdd(_p, _p, _s, alpha);
 scaledAdd(_r, _r, _z, -alpha);
```

公式中的v就是代码alpha，公式中的x就是代码的中的p压力。公式中的s就是代码中的s。
$$
\mu  = \frac{(r_n,r_n)}{(r_{n-1},r_{n-1}}\\
s = r_n + \mu * s
$$


```c++
double sigmaNew = dotProduct(_z, _r);
scaledAdd(_s, _z, _s, sigmaNew/sigma);
sigma = sigmaNew;
```

代码中的_z才是残差，代码中的r是当前残差，z是上一次残差。