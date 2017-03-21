cdef extern from "ab_initio_potential.h":
    double ab_initio_pot(double R, double Theta)

def potential(R, theta):
    return ab_initio_pot(R, theta)
