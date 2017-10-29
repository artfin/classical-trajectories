#pragma once

#include <cmath>

double ham_value( double R, double pR, double J );
void hamiltonian( double* out, double R, double pR );
double mol_frame_dipole( double R );
void rhs( double* out, double R, double pR, double J );


