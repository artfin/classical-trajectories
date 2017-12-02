#pragma once

#include <cmath>
#include <vector>
#include <tuple>
#include <iostream>

#include "ar_he_pes.h"
#include "ar_he_pes_derivative.h"
#include "ar_he_dip_buryak_fit.h"

#include <Eigen/Dense>

void rhs( double* out, double R, double pR, double theta, double pTheta );
void transform_dipole( std::vector<double> &res, double R, double theta );
void transform_coordinates( std::tuple<double, double, double> &he_coords, 
							std::tuple<double, double, double> &ar_coords, 
							const double &R, const double &theta );
