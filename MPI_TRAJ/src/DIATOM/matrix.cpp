#include "matrix.h"

#include <iostream>

const double HE_MASS = 4.00260325413;
const double AR_MASS = 39.9623831237; 
const double PROTON_TO_ELECTRON_RATIO = 1836.15267389; 

const double MU = HE_MASS * AR_MASS / ( HE_MASS + AR_MASS ) * PROTON_TO_ELECTRON_RATIO; 

void transform_dipole( std::vector<double> &res, double R, double theta )
{
	// constructing simple S matrix
	Eigen::Matrix<double, 3, 3> S;
	S(0, 0) = 1.0;
	S(0, 1) = 0.0;
	S(0, 2) = 0.0;

	S(1, 0) = 0.0;
	S(1, 1) = cos( theta );
	S(1, 2) = -sin( theta );

	S(2, 0) = 0.0;
	S(2, 1) = sin(theta);
	S(2, 2) = cos( theta );

	double dipz = ar_he_dip( R );
	Eigen::Vector3d mol_dipole( 0, 0, dipz );
	Eigen::Vector3d lab_dipole = S * mol_dipole;

	res[0] = lab_dipole[0];
	res[1] = lab_dipole[1];
	res[2] = lab_dipole[2];
}

void rhs( double* out, double R, double pR, double theta, double pTheta )
{
	out[0] = pR / MU; // dot{R} 
	out[1] = pow(pTheta, 2) / MU / pow(R, 3) - ar_he_pot_derivative( R ); // dot{pR}
	out[2] = pTheta / MU / pow(R, 2); // dot{theta}
	out[3] = 0; // dot[pTheta}
}
