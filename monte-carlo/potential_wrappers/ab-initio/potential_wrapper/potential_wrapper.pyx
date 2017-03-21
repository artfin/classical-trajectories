cdef extern from "ab_initio_potential.h":
    double ab_initio_pot(double R, double Theta)

def potential(R, theta):
    return ab_initio_pot(R, theta) /  2.1947 / 10**(5) # from cm^(-1) to hartree
