cdef extern from "derivatives.h":
    double arco2_pes_derivTheta(double R, double Theta)
    double arco2_pes_derivR(double R, double Theta)

def derivative_R(R, theta):
    return arco2_pes_derivR(R, theta)

def derivative_Theta(R, theta):
    return arco2_pes_derivTheta(R, theta)
