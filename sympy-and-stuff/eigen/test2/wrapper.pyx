import cython

# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

# declare the interface to the C code

cdef extern from 'hamiltonian.h':
    void rhs(double* out, double R, double theta, double pR, double pT, double alpha, double beta, double J)
    void hamiltonian(double *out, double R, double theta, double pR, double pT, double alpha, double beta, double J)
    
@cython.boundscheck(False)
@cython.wraparound(False)
def __rhs__(np.ndarray[double, ndim = 1, mode = "c"] out, double R, double theta, double pR, double pT, double alpha, double beta):
    """
        Input: 
            np.array() of length 6 to be filled with right-hand sides
            dynamic variables: R, theta, pR, pT, alpha, beta
        RHS:
            are put in the given np.array()
    """
    J = 1.0 
    rhs(&out[0], R, theta, pR, pT, alpha, beta, J)
    return None

@cython.boundscheck(False)
@cython.wraparound(False)
def __derivatives__(np.ndarray[double, ndim = 1, mode = "c"] out, double R, double theta, double pR, double pT, double alpha, double beta, double J):
    """ 
        Input:
            np.array() of length 6 to be filled with derivatives of hamiltonian by dynamic variables
            dynamic variable: R, theta, pR, pT, alpha, beta, J
    """
    hamiltonian(&out[0], R, theta, pR, pT, alpha, beta, J)
    return None
