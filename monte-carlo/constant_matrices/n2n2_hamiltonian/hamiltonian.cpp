#include <iostream>
#include <math.h>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

const double dalton = 1.660539040 * pow(10, -27); // dalton to kg
const double amu = 9.1093835 * pow(10, -31); // amu to  kg 

const double mu1 = 7.0035 * dalton / amu; // amu
const double mu2 = mu1; // amu
const double mu3 = 14.007 * dalton / amu; // amu
const double l1 = 2.0744; // bohrs
const double l2 = l1; // bohrs

const double l1_sq = l1 * l1;
const double l2_sq = l2 * l2;

void fill_inertia_tensor(Matrix3d &inertia_tensor, double &q1, double &q2, double &q3, double &q4)
{
    double sin_q1 = sin(q1);
    double cos_q1 = cos(q1);
    double sin_q1_sq = sin_q1 * sin_q1;
    double cos_q1_sq = cos_q1 * cos_q1;

    double sin_q2 = sin(q2);
    double cos_q2 = cos(q2);
    double sin_q2_sq = sin_q2 * sin_q2;
    double cos_q2_sq = cos_q2 * cos_q2;

    double sin_q3 = sin(q3);
    double cos_q3 = cos(q3);    
    double sin_q3_sq = sin_q3 * sin_q3;
    double cos_q3_sq = cos_q3 * cos_q3;

    double q4_sq = q4 * q4;

    inertia_tensor(0, 0) = mu1 * l1_sq * cos_q1_sq + mu2 * l2_sq * (sin_q2_sq  * sin_q3_sq + cos_q3_sq) + mu3 * q4_sq; 
    inertia_tensor(1, 0) = - mu2 * l2_sq * sin_q2 * cos_q2 * sin_q3_sq;
    inertia_tensor(2, 0) = - mu1 * l1_sq * sin_q1 * cos_q1 - mu2 * l2_sq * cos_q2 * sin_q3 * cos_q3;

    inertia_tensor(0, 1) = inertia_tensor(1, 0);
    inertia_tensor(1, 1) = mu1 * l1_sq + mu2 * l2_sq * (cos_q2_sq * sin_q3_sq + cos_q3_sq) + mu3 * q4_sq;
    inertia_tensor(2, 1) = - mu2 * l2_sq * sin_q2 * sin_q3 * cos_q3;
    
    inertia_tensor(0, 2) = inertia_tensor(2, 0);
    inertia_tensor(1, 2) = inertia_tensor(2, 1);
    inertia_tensor(2, 2) = mu1 * l1_sq * sin_q1_sq + mu2 * l2_sq * sin_q3_sq;
}

void fill_a_matrix(Matrix4d &a, double &q1, double &q2, double &q3, double &q4)
{
    double sin_q3 = sin(q3);

    a(0, 0) = mu1 * l1_sq;
    a(1, 1) = mu2 * l2_sq * sin_q3 * sin_q3;
    a(2, 2) = mu2 * l2_sq;
    a(3, 3) = mu3;

    a(0, 1) = a(0, 2) = a(0, 3) = 0;
    a(1, 0) = a(1, 2) = a(1, 3) = 0;
    a(2, 0) = a(2, 1) = a(2, 3) = 0;
    a(3, 0) = a(3, 1) = a(3, 2) = 0;
}

void fill_A_matrix(Matrix<double, 3, 4> &A, double &q1, double &q2, double &q3, double &q4)
{
    double sin_q2 = sin(q2);
    double cos_q2 = cos(q2);

    double sin_q3 = sin(q3);
    double cos_q3 = cos(q3);

    A(0, 0) = 0;
    A(1, 0) = mu1 * l1_sq;
    A(2, 0) = 0;

    A(0, 1) = - mu2 * l2_sq * cos_q2 * sin_q3 * cos_q3;
    A(1, 1) = - mu2 * l2_sq * sin_q2 * sin_q3 * cos_q3;
    A(2, 1) = mu2 * l2_sq * sin_q3 * sin_q3;

    A(0, 2) = - mu2 * l2_sq * sin_q2;
    A(1, 2) = mu2 * l2_sq * cos_q2;
    A(2, 2) = 0;

    A(0, 3) = 0;
    A(1, 3) = 0;
    A(2, 3) = 0;
}

double kinetic_energy(double q1, double q2, double q3, double q4, double p1, double p2, double p3, double p4, double Jx, double Jy, double Jz)
{
	Vector3d j_vector(Jx, Jy, Jz);
	Vector4d p_vector(p1, p2, p3, p4);

	Matrix<double, 3, 3> I;
	Matrix<double, 4, 4> a;
	Matrix<double, 3, 4> A;

	fill_inertia_tensor(I, q1, q2, q3, q4);
	fill_a_matrix(a, q1, q2, q3, q4);
	fill_A_matrix(A, q1, q2, q3, q4);

	Matrix<double, 3, 3> I_inv = I.inverse();
	Matrix<double, 4, 4> a_inv = a.inverse();

	Matrix<double, 3, 3> G11;
	Matrix<double, 4, 4> G22;
	Matrix<double, 3, 4> G12;

	Matrix<double, 3, 3> t1;
	Matrix<double, 4, 4> t2;

	t1 = I;
	t1.noalias() -= A * a_inv * A.transpose();
	G11 = t1.inverse();

	t2 = a;
	t2.noalias() -= A.transpose() * I_inv * A;
	G22 = t2.inverse();

	G12.noalias() = - G11 * A * a.inverse();

	//cout << "p vector: " << endl << p_vector << endl;
	//cout << "a: " << endl << a << endl;
	//cout << "G22: " << endl << G22 << endl;

	double ang_term = 0.5 * j_vector.transpose() * G11 * j_vector;
	double kin_term = 0.5 * p_vector.transpose() * G22 * p_vector;
	double cor_term = j_vector.transpose() * G12 * p_vector;

	//cout << "angular term: " << ang_term << endl;
	//cout << "kin term: " << kin_term << endl;
	//cout << "cor term: " << cor_term << endl;

	return ang_term + kin_term + cor_term; 
}

//int main()
//{
	//double q1 = 1.95;
	//double q2 = 5.28;
	//double q3 = 0.14;
	//double q4 = 4.44;

	//double p1 = -24.04;
	//double p2 = 49.58;
	//double p3 = 36.22;
	//double p4 = 1.48;

	//double Jx = -4.20;
	//double Jy = -22.97;
	//double Jz = -2.96;

	//double h = kinetic_energy(q1, q2, q3, q4, p1, p2, p3, p4, Jx, Jy, Jz);
	//cout << "hamiltonian value: " << h << endl;

	//return 0;
//}






