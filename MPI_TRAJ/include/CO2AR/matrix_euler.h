#pragma once

#include <cmath>
#include <Eigen/Dense>

#include "psp_pes.h"
#include "co2_ar_dipole.h"

using Eigen::Matrix;
using Eigen::Matrix2d;
using Eigen::Matrix3d;

using Eigen::Vector2d;
using Eigen::Vector3d;

void inertia_tensor(Matrix3d &inertia_tensor, const double &R, const double &Theta);
void inertia_tensor_dR(Matrix3d &i_dr, double &R, double &Theta);
void inertia_tensor_dTheta(Matrix3d &i_dt, const double &R, const double &Theta);

void fill_a_matrix(Matrix2d &a, const double &R, const double &theta);
void fill_A_matrix(Matrix<double, 3, 2> &A, const double &R, const double &theta);

void W_matrix( Matrix<double, 3, 3> &W, const double &theta, const double &psi );
void dW_dpsi( Matrix<double, 3, 3> &dW, const double &theta, const double &psi );
void dW_dtheta( Matrix<double, 3, 3> &dW, const double &theta, const double &psi );

void rhs(double* out, const double &R, const double &Theta, 
					  const double &pR, const double &pT, 
					  const double &phi, const double &theta, const double &psi, 
					  const double &p_phi, const double &p_theta, const double &p_psi );

void transform_dipole( std::vector<double> &output, const double &R,
							   					    const double &Theta,
											   		const double &phi,
											   		const double &theta,
											   		const double &psi );


