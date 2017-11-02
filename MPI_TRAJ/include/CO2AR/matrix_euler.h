#pragma once

#include <cmath>
#include <Eigen/Dense>

#include "psp_pes.h"

using Eigen::Matrix;
using Eigen::Matrix2d;
using Eigen::Matrix3d;

using Eigen::Vector2d;
using Eigen::Vector3d;

void fill_inertia_tensor(Matrix3d &inertia_tensor, double &R, double &theta);
void fill_a_matrix(Matrix2d &a, double &R, double &theta);
void fill_A_matrix(Matrix<double, 3, 2> &A, double &R, double &theta);

void W_matrix( Matrix<double, 3, 3> &W, const double &theta, const double &psi );
void dW_dpsi( Matrix<double, 3, 3> &dW, const double &theta, const double &psi );
void dW_dtheta( Matrix<double, 3, 3> &dW, const double &theta, const double &psi );

void rhs(double* out, double &R, double &Theta, 
					  double pR, double pT, 
					  double phi, double theta, double psi, 
					  double p_phi, double p_theta, double p_psi );


