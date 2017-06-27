#include <iostream>
#include <math.h>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

const double m = 1;
const double r0 = 1;
const double I0 = m * r0 * r0;

// matrix 'a' for one-dimensional water is constant
// matrix 'A' is zero-matrix
double a = 0.5 * I0;

void inertia_tensor(Matrix3d &inertia_tensor, double &q)
{
	double cos_q = cos(q);

	inertia_tensor(0, 0) = I0 * (1 - cos_q);
	inertia_tensor(0, 1) = inertia_tensor(0, 2) = 0;

	inertia_tensor(1, 1) = 2 * I0;
	inertia_tensor(1, 0) = inertia_tensor(1, 2) = 0;

	inertia_tensor(2, 2) = I0 * (1 + cos_q);
	inertia_tensor(2, 0) = inertia_tensor(2, 1) = 0;
}

void inertia_tensor_dq(Matrix3d &i_dq, double &q)
{
	double sin_q = sin(q);

	i_dq(0, 0) = I0 * sin_q;
	i_dq(0, 1) = i_dq(0, 2) = 0;

	i_dq(1, 0) = i_dq(1, 1) = i_dq(1, 2) = 0;

	i_dq(2, 2) = - i_dq(0, 0);
	i_dq(2, 0) = i_dq(2, 1) = 0;
}

void hamiltonian(double *out, double q, double p_q, double alpha, double beta, double J)
{
	Vector3d j_vector(J * cos(alpha) * sin(beta),
			  J * sin(alpha) * sin(beta),
			  J * cos(beta));

	Vector3d j_vector_dalpha( - J * sin(alpha) * sin(beta),
				    J * cos(alpha) * sin(beta),
				    0);
	Vector3d j_vector_dbeta( J * cos(alpha) * cos(beta),
				 J * sin(alpha) * cos(beta),
				-J * sin(beta));

	// initializing and filling inertia tensor 
	Matrix<double, 3, 3> I;
	Matrix3d I_dq;
	inertia_tensor(I, q);
	inertia_tensor_dq(I_dq, q);

	Matrix3d G11 = I.inverse();
	Matrix3d G11_dq = - G11 * I_dq * G11;

	// calculating derivative dH/dq
	double h_dq = 0.5 * j_vector.transpose() * G11_dq * j_vector;
	
	// calculating derivative dH/dp
	double h_dp = 2 * p_q / I0;

	// derivatives dH/dalpha, dH/dbeta
	double h_dalpha = j_vector.transpose() * G11 * j_vector_dalpha;
	double h_dbeta = j_vector.transpose() * G11 * j_vector_dbeta;

	out[0] = h_dq;
	out[1] = h_dp;
	out[2] = h_dalpha;
	out[3] = h_dbeta;
}

int main()
{
	double *out = new double[4];
	hamiltonian(out, 1.2, 3.0, 0.55, 0.1, 30);

	for (int i = 0; i < 4; i++) {
		cout << out[i] << endl;
	}

	return 0;

}
