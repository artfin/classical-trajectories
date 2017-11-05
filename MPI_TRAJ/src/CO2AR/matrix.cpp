#include "matrix.h"
#include "co2_ar_dipole.h"
#include "psp_pes.h"

using namespace Eigen;
using namespace std;

const double mu1 = 14579.0;
const double mu2 = 38183.0;
const double l = 4.398; 

void inertia_tensor_dR(Matrix3d &i_dr, double &R, double &Theta)
{
    i_dr(0, 0) = i_dr(1, 1) = 2 * mu2 * R;
    
    i_dr(1, 0) = i_dr(2, 0) = 0;    
    i_dr(0, 2) = i_dr(1, 2) = i_dr(2, 2) = 0;
    i_dr(0, 1) = i_dr(2, 1) = 0;
}

void inertia_tensor_dTheta(Matrix3d &i_dt, double &R, double &Theta)
{
    double l2 = l * l;
    double sin_t = sin(Theta);
    double cos_t = cos(Theta);

    i_dt(0, 0) = - 2 * mu1 * l2 * sin_t * cos_t;
    i_dt(1, 0) = 0;
    i_dt(2, 0) = i_dt(0, 2) = - mu1 * l2 * (cos_t * cos_t - sin_t * sin_t);

    i_dt(1, 2) = 0;
    i_dt(2, 2) = - i_dt(0, 0); 

    i_dt(0, 1) = i_dt(1, 1) = i_dt(2, 1) = 0;
}

void inertia_tensor(Matrix3d &inertia_tensor, double &R, double &Theta)
{
    double sin_t = sin(Theta);
    double cos_t = cos(Theta);
    double l2 = l * l;
    double R2 = R * R;

    inertia_tensor(0, 0) = mu1 * l2 * cos_t * cos_t + mu2 * R2; 
    inertia_tensor(1, 0) = 0;
    inertia_tensor(2, 0) = -mu1 * l2 * sin_t * cos_t;

    inertia_tensor(0, 2) = inertia_tensor(2, 0);
    inertia_tensor(1, 2) = 0;
    inertia_tensor(2, 2) = mu1 * l2 * sin_t * sin_t;

    inertia_tensor(0, 1) = 0;
    inertia_tensor(1, 1) = inertia_tensor(0, 0) + inertia_tensor(2, 2);
    inertia_tensor(2, 1) = 0;
}

void fill_a_matrix(Matrix2d &a, double &R, double &Theta)
{
    a(0, 0) = mu2;
    a(0, 1) = 0;
    a(1, 0) = 0;
    a(1, 1) = mu1 * l * l;
}

void fill_A_matrix(Matrix<double, 3, 2> &A, double &R, double &Theta)
{
    A(0, 0) = A(0, 1) = 0;
    A(1, 0) = 0;
    A(1, 1) = mu1 * l * l;
    A(2, 0) = A(2, 1) = 0;
}

double ham_value(double R, double theta, double pR, double pT, double Jx, double Jy, double Jz)
{
	Vector3d j_vector(Jx, Jy, Jz);
	Vector2d p_vector(pR, pT);

	Matrix<double, 3, 3> I;
	Matrix<double, 2, 2> a;
	Matrix<double, 3, 2> A;

	inertia_tensor(I, R, theta);
	fill_a_matrix(a, R, theta);
	fill_A_matrix(A, R, theta);

	Matrix<double, 3, 3> I_inv = I.inverse();
	Matrix<double, 2, 2> a_inv = a.inverse();

	Matrix<double, 3, 3> G11;
	Matrix<double, 2, 2> G22;
	Matrix<double, 3, 2> G12;

	Matrix<double, 3, 3> t1;
	Matrix<double, 2, 2> t2;

	t1 = I;
	t1.noalias() -= A * a_inv * A.transpose();
	G11 = t1.inverse();

	t2 = a;
	t2.noalias() -= A.transpose() * I_inv * A;
	G22 = t2.inverse();

	G12.noalias() = - G11 * A * a.inverse();

	double ang_term = 0.5 * j_vector.transpose() * G11 * j_vector;
	double kin_term = 0.5 * p_vector.transpose() * G22 * p_vector;
	double cor_term = j_vector.transpose() * G12 * p_vector;

	return ang_term + kin_term + cor_term + psp_pes(R, theta); 
}


void hamiltonian(double* out, double R, double Theta, double pR, double pT, double phi, double theta, double J, bool dip_calc)
{
    //declaring angular momentum and its derivatives
    Vector3d j_vector(J * cos(phi) * sin(theta),
                      J * sin(phi) * sin(theta),
                      J * cos(theta));
   
    Vector3d j_vector_dphi(- J * sin(phi) * sin(theta),
                             J * cos(phi) * sin(theta),
                             0);
    Vector3d j_vector_dtheta(J * cos(phi) * cos(theta),
                             J * sin(phi) * cos(theta),
                            -J * sin(theta));
    
    Vector2d p_vector(pR, pT);

    // declaring inertia tensor and it's derivatives
    Matrix<double, 3, 3> I; 
    Matrix<double, 3, 3> I_dr;
    Matrix<double, 3, 3> I_dTheta;
    
    // declaring matrices a, A
    Matrix<double, 2, 2> a;
    Matrix<double, 3, 2> A;

    // filling matrices 
    inertia_tensor(I, R, Theta); 
    inertia_tensor_dR(I_dr, R, Theta);
    inertia_tensor_dTheta(I_dTheta, R, Theta);
    fill_a_matrix(a, R, Theta);
    fill_A_matrix(A, R, Theta);

    // supplementary matrices
    Matrix<double, 3, 3> I_inv = I.inverse();
    Matrix<double, 2, 2> a_inv = a.inverse();

    Matrix<double, 3, 3> G11;
    Matrix<double, 2, 2> G22;
    Matrix<double, 3, 2> G12;
   
    Matrix<double, 3, 3> t1;
    Matrix<double, 2, 2> t2;
    
    t1 = I;
    t1.noalias() -= A * a_inv * A.transpose();
    G11 = t1.inverse();
  
    t2 = a;
    t2.noalias() -= A.transpose() * I_inv * A;
    G22 = t2.inverse();

    G12.noalias() = - G11 * A * a.inverse();
    
    Matrix<double, 3, 3> G11_dr;
    Matrix<double, 3, 3> G11_dTheta;
    Matrix<double, 2, 2> G22_dr;
    Matrix<double, 2, 2> G22_dTheta;
    Matrix<double, 3, 2> G12_dr;
    Matrix<double, 3, 2> G12_dTheta;

    
    // (vector of) derivatives dH/dp
    Vector2d h_dp = G22 * p_vector + G12.transpose() * j_vector;
    
    // usual derivative calculation
    if ( dip_calc == false )
    {
    	// derivatives dH/dr and dH/dTheta
    	G11_dr = - G11 * I_dr * G11;
    	G11_dTheta = - G11 * I_dTheta * G11;
   
    	G22_dr = - G22 * A.transpose() * I_inv * I_dr * I_inv * A * G22; 
    	G22_dTheta = - G22 * A.transpose() * I_inv * I_dTheta * I_inv * A * G22;

    	G12_dr = - G11_dr * A * a_inv;
    	G12_dTheta = -G11_dTheta * A * a_inv;

    	double ang_term, kin_term, cor_term;

    	ang_term = 0.5 * j_vector.transpose() * G11_dr * j_vector;
    	kin_term = 0.5 * p_vector.transpose() * G22_dr * p_vector;
    	cor_term = j_vector.transpose() * G12_dr * p_vector;
    	double h_dr = ang_term + kin_term + cor_term;
    
   		ang_term = 0.5 * j_vector.transpose() * G11_dTheta * j_vector;
		kin_term = 0.5 * p_vector.transpose() * G22_dTheta * p_vector;
		cor_term = j_vector.transpose() * G12_dTheta * p_vector;
		double h_dTheta = ang_term + kin_term + cor_term;
    	
		// derivatives dH/dphi and dH/dtheta
		ang_term = j_vector.transpose() * G11 * j_vector_dphi;
		cor_term = j_vector_dphi.transpose() * G12 * p_vector;
		double h_dphi = ang_term + cor_term;

		ang_term = j_vector.transpose() * G11 * j_vector_dtheta;
		cor_term = j_vector_dtheta.transpose() * G12 * p_vector;
		double h_dtheta = ang_term + cor_term;
	    
		out[0] = h_dr; 
		out[1] = h_dTheta;
		out[2] = h_dp(0);
		out[3] = h_dp(1);
		out[4] = h_dphi;
		out[5] = h_dtheta;
    }
    	
    // calculating derivatives of dipole moment in laboratory frame
    else
    {
    	Vector3d omega = G11 * j_vector + G12 * p_vector;
	
		// ! R, Theta -- internal coordinates ! 
		// values of dipole components
		double dipole_x = dipx(R, Theta);
		double dipole_y = 0;
		double dipole_z = dipz(R, Theta);

		// values of dipole derivatives
		double ddipolex_dR = ddipxdR(R, Theta);
		double ddipoley_dR = 0; 
		double ddipolez_dR = ddipzdR(R, Theta);

		double ddipolex_dTheta = ddipxdTheta(R, Theta);
		double ddipoley_dTheta = 0; 
		double ddipolez_dTheta = ddipzdTheta(R, Theta);

		// derivatives of dipole in molecular frame
		// h_dp(0) = dH/dpR, \dot{R} = dH/dpR
		// h_dp(1) = dH/dpTheta, \dot{\theta} = dH/dpTheta
		double ddipolex_dt = ddipolex_dR * h_dp(0) + ddipolex_dTheta * h_dp(1);
		double ddipoley_dt = ddipoley_dR * h_dp(0) + ddipoley_dTheta * h_dp(1);
		double ddipolez_dt = ddipolez_dR * h_dp(0) + ddipolez_dTheta * h_dp(1);

		// derivatives of dipole in laboratory frame
		double ddipolex_dt_lab = ddipolex_dt + omega(1) * dipole_z - omega(2) * dipole_y;
		double ddipoley_dt_lab = ddipoley_dt + omega(2) * dipole_x - omega(0) * dipole_z;
		double ddipolez_dt_lab = ddipolez_dt + omega(0) * dipole_y - omega(1) * dipole_x;
   	
		out[0] = ddipolex_dt_lab;
		out[1] = ddipoley_dt_lab;
		out[2] = ddipolez_dt_lab;
    }
}

void rhs(double* out, double R, double Theta, double pR, double pT, double phi, double theta, double J)
// input:
//      out -- prepared array to fill the right-hand sides of ODEs
{
    double Jsint = J * sin(theta);

    double* derivatives = new double[6];
    hamiltonian(derivatives, R, Theta, pR, pT, phi, theta, J, false);
    
    out[0] = derivatives[2]; // /dot(R) = dH/dpR
    out[1] = derivatives[3]; // /dot(Theta) = dH/dpT
    out[2] = - derivatives[0] - dpsp_pesdR(R,Theta);  // /dot(pR) = - dH/dR = - dT/dR - dU/dR (to hartrees from cm^-1)
    out[3] = - derivatives[1] - dpsp_pesdTheta(R, Theta); // /dot(pT) = - dH/dTheta = -dT/dTheta - dU/dTheta (to hartrees from cm^-1)
    out[4] = 1 / Jsint * derivatives[5]; // /dot(varphi) = 1 / J / sin(teta) * dH/dTheta
    out[5] = - 1 / Jsint * derivatives[4]; // /dot(Theta) = - 1 / J / sin(Theta) * dH/dvarphi

    delete[] derivatives;
}
