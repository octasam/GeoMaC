
import numpy as np
from scipy.sparse.linalg import cg
import time

def steepest(A,b,x,tol,kmax):
    r = b - A * x
    res = np.linalg.norm(r)
    k = 0
    while(res > tol and k < kmax):
        alpha = r.T * r / (r.T * A * r)
        x = x + r * alpha
        r = b - A * x
        res = np.linalg.norm(r)
        k += 1
    return x, k

def conjugateGrad(A, b, x=None):
    """
    @brief Solve a linear equation Ax = b with conjugate gradient method.

    @param A: 2d numpy.array of positive semi-definite (symmetric) matrix
    @param b: 1d numpy.array rhs.
    @return x: 1d numpy.array solution.
    """
    n = len(b)
    if not x:
        x = np.ones(n)
    r = np.dot(A, x) - b
    p = - r
    r_k_norm = np.dot(r, r)
    for i in range(2*n):
        Ap = np.dot(A, p)
        alpha = r_k_norm / np.dot(p, Ap)
        x += alpha * p
        r += alpha * Ap
        r_kplus1_norm = np.dot(r, r)
        beta = r_kplus1_norm / r_k_norm
        r_k_norm = r_kplus1_norm
        if r_kplus1_norm < 1e-5:
            print('Itr:', i)
            break
        p = beta * p - r
    return x

if __name__ == '__main__':
    n = 1000
    P = np.random.normal(size=[n, n])
    A = np.dot(P.T, P)
    b = np.ones(n)

    t0 = time.time()
    print('start')
    x0, k = steepest(A, b, b, 0.001, 20)
    t1 = time.time()
    print('k = ', k)
    print(t1 - t0)
    x1 = conjugateGrad(A, b)
    t2 = time.time()
    print(t2 - t1)
    x2 = np.linalg.solve(A, b)
    t3 = time.time()
    print(t3 - t2)
    x3 = cg(A, b)
    t4 = time.time()
    print(t4 - t3)

# print np.dot(A, x) - b
