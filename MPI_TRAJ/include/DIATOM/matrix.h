#pragma once

#include <cmath>

double ham_value( double R, double pR, double Jx, double Jy );
void hamiltonian( double* out, double R, double pR, bool dip_calc );
void rhs( double* out, double R, double pR, double J );


