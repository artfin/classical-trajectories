# file: cython_integrand.pyx

import numpy as np

# use exp from C
from libc.math cimport exp

def f_batch(double[:, ::1] x):
    cdef int i # labels integration point
    cdef int d # labels direction
    cdef int dim = x.shape[1]
    cdef double norm = 1013.2118364296088 ** (dim / 4.)
    cdef double dx2
    cdef double[::1] ans = np.empty(x.shape[0], float)
    
    for i in range(x.shape[0]):
        # integrand for integration point x[i]
        dx2 = 0.0
        for d in range(dim):
            dx2 += (x[i, d] - 0.5) ** 2
        ans[i] = exp(-100. * dx2) * norm
    return ans
