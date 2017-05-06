cdef extern from 'hamiltonian.h':
    double kinetic_energy(double q1, double q2, double q3, double q4, double p1, double p2, double p3, double p4, double Jx, double Jy, double Jz)

def ke(q1, q2, q3, q4, p1, p2, p3, p4, Jx, Jy, Jz):
    return kinetic_energy(q1, q2, q3, q4, p1, p2, p3, p4, Jx, Jy, Jz)
