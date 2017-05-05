cdef extern from 'hamiltonian.h':
    double hamiltonian(double R, double theta, double pR, double pT, double Jx, double Jy, double Jz)

def __hamiltonian__(R, theta, pR, pT, Jx, Jy, Jz):
    return hamiltonian(R, theta, pR, pT, Jx, Jy, Jz)
