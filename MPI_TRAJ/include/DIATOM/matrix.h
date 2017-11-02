#pragma once

#include <cmath>
#include <vector>
#include <iostream>

#include "ar_he_pes.h"
#include "ar_he_pes_derivative.h"
#include "ar_he_dip.h"
#include "ar_he_dip_derivative.h"

#include <Eigen/Dense>

void rhs( double* out, double R, double pR, double theta, double pTheta );
void transform_dipole( std::vector<double> &res, double R, double theta );

