#include <iostream>
#include <math.h>
#include <Eigen/Dense>

using Eigen::Matrix4d;
using Eigen::Matrix3d;
using Eigen::Matrix;

void fill_inertia_tensor(Matrix3d &inertia_tensor, double &q1, double &q2, double &q3, double &q4);
void fill_a_matrix(Matrix4d &a, double &q1, double &q2, double &q3, double &q4);
void fill_A_matrix(Matrix<double, 3, 4> &A, double &q1, double &q2, double &q3, double &q4);

double kinetic_energy(double q1, double q2, double q3, double q4, double p1, double p2, double p3, double p4, double Jx, double Jy, double Jz); 
