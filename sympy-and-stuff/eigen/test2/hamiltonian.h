#include <iostream>
#include <math.h>
#include <Eigen/Dense>

using Eigen::Matrix3d;
using Eigen::Matrix2d;
using Eigen::Matrix;

void inertia_tensor_dR(Matrix3d &i_dr, double &R, double &theta);
void inertia_tensor_dtheta(Matrix3d &i_dt, double &R, double &theta);
void inertia_tensor(Matrix3d &inertia_tensor, double &R, double &theta);
void fill_a_matrix(Matrix2d &a, double &R, double &theta);
void fill_A_matrix(Matrix<double, 3, 2> &A, double &R, double &theta);

double* legendre_array(int N, double cosT);
double derivativeTheta(double R, double Theta);
double derivativeR(double R, double Theta);

void hamiltonian(double* out, double R, double theta, double pR, double pT, double alpha, double beta, double J);

void rhs(double* out, double R, double theta, double pR, double pT, double alpha, double beta, double J);


