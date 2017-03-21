cdef extern from "ar_co2_pes_deriv.h":
    double arco2_pes_derivTheta(double& R, double& Theta)
    double arco2_pes_derivR(double& R, double& Theta)

def derivative_R(double[:] R, double[:] theta):
    return arco2_pes_derivR(&R[0], &theta[0])

def derivative_Theta(double[:] R, double[:] theta):
    return arco2_pes_derivTheta(&R[0], &theta[0])
