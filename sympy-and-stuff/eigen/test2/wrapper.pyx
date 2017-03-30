import cython

# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

# declare the interface to the C code

cdef extern from 'hamiltonian.h':
    void hamiltonian(double* out, double R, double theta, double pR, double pT, double J, double alpha, double beta)

@cython.boundscheck(False)
@cython.wraparound(False)
def __ham__(np.ndarray[double, ndim = 1, mode = "c"] out, double R, double theta, double pR, double pT, double J, double alpha, double beta):
    """
    Takes a numpy array as input 
    """
    hamiltonian(&out[0], R, theta, pR, pT, J, alpha, beta)
    return None
