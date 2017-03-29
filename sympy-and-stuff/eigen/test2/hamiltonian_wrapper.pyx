cdef extern from "hamiltonian.h":
    double* hamiltonian(double R, double theta, double pR, double pT, double alpha, double beta);

cdef h(double R, double theta, double pR, double pT, double alpha, double beta):
    return hamiltonian(R, theta, pR, pT, alpha, beta)
