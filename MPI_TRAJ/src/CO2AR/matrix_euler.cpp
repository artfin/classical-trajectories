#include "matrix_euler.h"

#include <iostream>

const double mu1 = 14579.0;
const double mu2 = 38183.0;
const double l = 4.398; 

void inertia_tensor_dR(Matrix3d &i_dr, const double &R, const double &Theta)
{
    i_dr(0, 0) = i_dr(1, 1) = 2 * mu2 * R;
    
    i_dr(1, 0) = i_dr(2, 0) = 0;    
    i_dr(0, 2) = i_dr(1, 2) = i_dr(2, 2) = 0;
    i_dr(0, 1) = i_dr(2, 1) = 0;
}

void inertia_tensor_dTheta(Matrix3d &i_dt, const double &R, const double &Theta)
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

void inertia_tensor(Matrix3d &inertia_tensor, const double &R, const double &Theta)
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

void fill_a_matrix(Matrix2d &a, const double &R, const double &Theta)
{
    a(0, 0) = mu2;
    a(0, 1) = 0;
    a(1, 0) = 0;
    a(1, 1) = mu1 * l * l;
}

void fill_A_matrix(Matrix<double, 3, 2> &A, const double &R, const double &Theta)
{
    A(0, 0) = A(0, 1) = 0;
    A(1, 0) = 0;
    A(1, 1) = mu1 * l * l;
    A(2, 0) = A(2, 1) = 0;
}

// J = (V^{-1})^\top pe = W pe
void W_matrix( Matrix<double, 3, 3> &W, const double &theta, const double &psi )
{
	double sin_psi = sin( psi );
	double cos_psi = cos( psi );

	double sin_theta = sin( theta );
	double cos_theta = cos( theta );
	double ctg_theta = cos_theta / sin_theta;

	W(0, 0) = sin_psi / sin_theta;
	W(0, 1) = cos_psi;
	W(0, 2) = - sin_psi * ctg_theta;
	
	W(1, 0) = cos_psi / sin_theta;
	W(1, 1) = - sin_psi;
	W(1, 2) = - cos_psi * ctg_theta;

	W(2, 0) = 0;
	W(2, 1) = 0;
	W(2, 2) = 1;
}

// \frac{\partial W}{\partial \psi}
void dW_dpsi( Matrix<double, 3, 3> &dW, const double &theta, const double &psi )
{
	double cos_psi = cos( psi );
	double sin_psi = sin( psi );

	double sin_theta = sin( theta );
	double cos_theta = cos( theta );
	double ctg_theta = cos_theta / sin_theta;

	dW(0, 0) = cos_psi / sin_theta;
   	dW(0, 1) = - sin_psi;
	dW(0, 2) = - cos_psi * ctg_theta;

	dW(1, 0) = - sin_psi / sin_theta;
	dW(1, 1) = - cos_psi;
	dW(1, 2) = sin_psi * ctg_theta;

	dW(2, 0) = 0;
	dW(2, 1) = 0;
	dW(2, 2) = 0;
}

// \frac{\partial W}{\partial \theta}
void dW_dtheta( Matrix<double, 3, 3> &dW, const double &theta, const double &psi )
{
	double sin_psi = sin( psi );
	double cos_psi = cos( psi );

	double sin_theta = sin( theta );
	double cos_theta = cos( theta );

	dW(0, 0) = - sin_psi * cos_theta / pow(sin_theta, 2);
	dW(0, 1) = 0;
	dW(0, 2) = sin_psi / pow(sin_theta, 2);

	dW(1, 0) = - cos_psi * cos_theta / pow(sin_theta, 2);
	dW(1, 1) = 0;
	dW(1, 2) = cos_psi / pow(sin_theta, 2);

	dW(2, 0) = 0;
	dW(2, 1) = 0;
	dW(2, 2) = 0;
}

void transform_dipole( std::vector<double> &output, const double &R,
								   			  	    const double &Theta,
											   	    const double &phi, 
											        const double &theta,
											        const double &psi )
{
	Matrix<double, 3, 3> S;
	
	double sin_phi = sin( phi );
	double cos_phi = cos( phi );

	double sin_theta = sin( theta );
	double cos_theta = cos( theta );

	double sin_psi = sin( psi );
	double cos_psi = cos( psi );

	S(0, 0) = cos_psi * cos_phi - cos_theta * sin_phi * sin_psi;
	S(0, 1) = - sin_psi * cos_phi - cos_theta * sin_phi * cos_psi;
	S(0, 2) = sin_theta * sin_phi;

	S(1, 0) = cos_psi * sin_phi + cos_theta * cos_phi * sin_psi;
	S(1, 1) = - sin_psi * sin_phi + cos_theta * cos_phi * cos_psi;
	S(1, 2) = - sin_theta * cos_phi;

	S(2, 0) = sin_theta * sin_psi;
	S(2, 1) = sin_theta * cos_psi;
	S(2, 2) = cos_theta;

	// vector of dipole in molecular frame
	Vector3d mol_dipole;

	mol_dipole(0) = dipx( R, Theta );
	mol_dipole(1) = 0;
	mol_dipole(2) = dipz( R, Theta );

	// vector of dipole in laboratory frame
	Vector3d lab_dipole = S * mol_dipole;

	output[0] = lab_dipole(0);
	output[1] = lab_dipole(1);
	output[2] = lab_dipole(2);
}

void rhs(double* out, const double &R, const double &Theta, 
					  const double &pR, const double &pT, 
					  const double &phi, const double &theta, const double &psi, 
					  const double &p_phi, const double &p_theta, const double &p_psi )
// input:
//     R, Theta, 
// 	   pR, pTheta, 
// 	   phi, theta, psi, 
// 	   p_phi, p_theta, p_psi
// output:
// 	   dH_dR, dH_dTheta
// 	   dH_dp(vector2d): dH_dp(0), dH_dp(1)
//	   dH_dpe(vector3d): dH_dpe(0), dH_dpe(1), dH_dpe(2)

{
	//std::cout << "inside rhs" << std::endl;
	//std::cout << "R: " << R << std::endl;
	//std::cout << "Theta: " << Theta << std::endl;

	// vector of euler impulses
  	Vector3d pe ( p_phi, p_theta, p_psi );

	// vector of impulses
    Vector2d p_vector(pR, pT);
	
	// ##################################################################
    // filling matrices I, a, A and derivatives of I
	Matrix<double, 3, 3> I; 
    inertia_tensor(I, R, Theta); 
	//std::cout << "I: " << I << std::endl;

	Matrix<double, 3, 3> I_dR;
    inertia_tensor_dR(I_dR, R, Theta);
    
	Matrix<double, 3, 3> I_dTheta;
    inertia_tensor_dTheta(I_dTheta, R, Theta);
    
    Matrix<double, 2, 2> a;
    fill_a_matrix(a, R, Theta);
    
	Matrix<double, 3, 2> A;
    fill_A_matrix(A, R, Theta);

	Matrix<double, 3, 3> W;
	W_matrix( W, theta, psi );
	// ##################################################################


	// ##################################################################
    // filling matrices G11, G22, G12
	Matrix<double, 3, 3> I_inv = I.inverse();
    Matrix<double, 2, 2> a_inv = a.inverse();

    Matrix<double, 3, 3> G11;
    Matrix<double, 2, 2> G22;
    Matrix<double, 3, 2> G12;

	//std::cout << "a: " << a << std::endl;
	//std::cout << "a_inv: " << a_inv << std::endl;
	//std::cout << "A: " << A << std::endl;

	Matrix<double, 3, 3> t1 = I;
    t1.noalias() -= A * a_inv * A.transpose();
    G11 = t1.inverse();
	//std::cout << "G11: " << G11 << std::endl;
 
    Matrix<double, 2, 2> t2 = a;	
	t2.noalias() -= A.transpose() * I_inv * A;
    G22 = t2.inverse();
	//std::cout << "G22: " << G22 << std::endl;

    G12.noalias() = - G11 * A * a.inverse();
	//std::cout << "G12: " << G12 << std::endl;
	// ##################################################################
	

	// ##################################################################
    // derivatives \frac{\partial \mH}{\partial p} 
	Vector2d dH_dp = G22 * p_vector + G12.transpose() * W * pe; 
	//std::cout << "dH_dp: " << dH_dp << std::endl;
	// ##################################################################


	// auxiliary variables
   	double ang_term, kin_term, cor_term;
 
	// ##################################################################
   	// derivative \frac{\partial \mH}{\partial R}
    Matrix<double, 3, 3> G11_dR;
	Matrix<double, 2, 2> G22_dR;
	Matrix<double, 3, 2> G12_dR;

   	G11_dR = - G11 * I_dR * G11;
   	G12_dR = - G11_dR * A * a_inv;
   	G22_dR = - G22 * A.transpose() * I_inv * I_dR * I_inv * A * G22; 
   
	//std::cout << "G11_dR: " << G11_dR << std::endl;
	//std::cout << "G12_dR: " << G12_dR << std::endl;
	//std::cout << "G22_dR: " << G22_dR << std::endl;	

	//std::cout << "W: " << W << std::endl;
	//std::cout << "pe: " << pe << std::endl;

   	ang_term = 0.5 * pe.transpose() * W.transpose() * G11_dR * W * pe;
	//std::cout << "ang_term: " << ang_term << std::endl;

   	kin_term = 0.5 * p_vector.transpose() * G22_dR * p_vector;
	//std::cout << "kin_term: " << kin_term << std::endl;

	cor_term = pe.transpose() * W.transpose() * G12_dR * p_vector;
	//std::cout << "cor_term: " << cor_term << std::endl;

	double dH_dR = ang_term + kin_term + cor_term;
	//std::cout << "dH_dR: " << dH_dR << std::endl;
	// ##################################################################


	// ##################################################################
   	// derivative \frac{\partial \mH}{\partial Theta}
	Matrix<double, 3, 3> G11_dTheta;
    Matrix<double, 2, 2> G22_dTheta;
    Matrix<double, 3, 2> G12_dTheta;

	G11_dTheta = - G11 * I_dTheta * G11;
	//std::cout << "G11_dTheta: " << G11_dTheta << std::endl;

	G22_dTheta = - G22 * A.transpose() * I_inv * I_dTheta * I_inv * A * G22;
	//std::cout << "G22_dTheta: " << G22_dTheta << std::endl;

	G12_dTheta = -G11_dTheta * A * a_inv;
	//std::cout << "G12_dTheta: " << G12_dTheta << std::endl;

   	ang_term = 0.5 * pe.transpose() * W.transpose() * G11_dTheta * W * pe;
	//std::cout << "W * pe: " << W * pe << std::endl;
	//std::cout << "ang_term: " << ang_term << std::endl;

	kin_term = 0.5 * p_vector.transpose() * G22_dTheta * p_vector;
	//std::cout << "kin_term: " << kin_term << std::endl;

	cor_term = pe.transpose() * W.transpose() * G12_dTheta * p_vector;
	//std::cout << "cor_term: " << cor_term << std::endl;

	double dH_dTheta = ang_term + kin_term + cor_term;	
	//std::cout << "dH_dTheta: " << dH_dTheta << std::endl;
	// ##################################################################

	// auxiliary variables
	double ang_term1, ang_term2;
	// ##################################################################
	// derivative /frac{\partial \mH}{\partial \theta} (euler angle)
	Matrix<double, 3, 3> W_dtheta;
   	dW_dtheta( W_dtheta, theta, psi );

	ang_term1 = 0.5 * pe.transpose() * W_dtheta.transpose() * G11 * W * pe;
	ang_term2 = 0.5 * pe.transpose() * W.transpose() * G11 * W_dtheta * pe;
	cor_term = pe.transpose() * W_dtheta.transpose() * G12 * p_vector;

	double dH_dtheta = ang_term1 + ang_term2 + cor_term;
	//std::cout << "dH_dtheta: " << dH_dtheta << std::endl;	
	// ##################################################################

	// ##################################################################
	// derivative \frac{\partial \mH}{\partial \psi} (euler angle)
	Matrix<double, 3, 3> W_dpsi;
	dW_dpsi( W_dpsi, theta, psi );

	ang_term1 = 0.5 * pe.transpose() * W_dpsi.transpose() * G11 * W * pe;
	ang_term2 = 0.5 * pe.transpose() * W.transpose() * G11 * W_dpsi * pe;
	cor_term = pe.transpose() * W_dpsi.transpose() * G12 * p_vector;

	double dH_dpsi = ang_term1 + ang_term2 + cor_term;
	//std::cout << "dH_dpsi: " << dH_dpsi << std::endl;
	// ##################################################################

	// ##################################################################
	// derivative \frac{\partial \mH}{\partial \phi} (euler angle)
	double dH_dphi = 0;
	// ##################################################################

	// ##################################################################
	// derivatives \frac{\partial \mH}{\partial pe} (euler impulses)
	
	Vector3d dH_dpe = W.transpose() * G11 * W * pe + W.transpose() * G12 * p_vector;
	//std::cout << "dH_dpe: " << dH_dpe << std::endl;
	// ##################################################################
  
	out[0] = dH_dp(0); 
	out[1] = dH_dp(1); 
	out[2] = - dH_dR - dpsp_pesdR( R, Theta ); 
	out[3] = - dH_dTheta - dpsp_pesdTheta( R, Theta );

	out[4] = dH_dpe(0);
	out[5] = dH_dpe(1);
   	out[6] = dH_dpe(2);

	out[7] = - dH_dphi; 
	out[8] = - dH_dtheta;
	out[9] = - dH_dpsi;
}
