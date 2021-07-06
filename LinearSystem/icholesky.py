import numpy as np

def icholesky(a):
    n = a.shape[0]
    L = np.zeros((n,n))
    for k in range(n):
        a[k,k] = np.sqrt(a[k,k])
        for i in range(k+1,n):
            if (a[i,k] !=0):
                a[i,k] = a[i,k]/a[k,k]
        # i,_= a[:,k].nonzero()
        # if len(i) > 0:
        #     a[i,k] = a[i,k]/a[k,k]
        for j in range(k+1,n):
            for i in range(j,n):
                if (a[i,j]!=0):
                    a[i,j] = a[i,j]-a[i,k]*a[j,k]  
            # i,_ = a[j:,j].nonzero()
            # if len(i) > 0: 
            #     a[i,j]  = a[i,j] - a[i,k]*a[j,k]     
        L[k:,k] = a[k:,k]
    return L

def modicholesky(a):
    n = a.shape[0]
    d = np.zeros((n))
    L = np.zeros((n,n))
    d[0] = 1 / a[0,0]
    for i in range(1,n):
        for j in range(i):
            for k in range(i,n):
                a[k,i] = a[k,i] - a[k,i] * a[i,j] * d[j]
        for k in range(i+1,n):
            a[i,i] = a[i,i] + abs(a[k,i])
            a[k,k] = a[k,k] + abs(a[k,i])
            a[k,i] = 0
        d[i] = 1 / a[i,i]
            
    return L

def Modified2(A):
    # A = {L}{d}{L^T}
    n = A.shape[0]
    d = np.zeros((n))
    L = np.zeros((n,n))
    L[0,0] = A[0,0]
    d[0] = 1 / L[0,0]
    for i in range(1,n):
        for j in range(0,i+1):
            lld = A[i,j]
            for k in range(0,j):
                lld -= L[i,k]*L[j,k]*d[k]
            L[i,j] = lld
        d[i] = 1 / L[i,i]
    return L,d

L = np.array([[1,0,0],[2,3,0],[4,5,7]])
A = np.dot(L,np.transpose(L))
# B = icholesky(A)
L,d = Modified2(A)
A2 = np.dot(L,d)
A3 = np.dot(np.transpose(A2),np.transpose((L)))