cdef extern from "ar_co2_pes_deriv.h":
    double arco2_pes_derivTheta(double R, double Theta)
    double arco2_pes_derivR(double R, double Theta)

def derivative_R(R, theta):
    return arco2_pes_derivR(R, theta) / 2.1947 / 10**(5) # cm^-1 to hartree

def derivative_Theta(R, theta):
    return arco2_pes_derivTheta(R, theta) / 2.1947 / 10**(5) # cm^-1 to hartree
