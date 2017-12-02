#include "matrix.h"


const double DALTON_UNIT = 1.660539040 * 1e-27;
const double AMU = 9.10938356 * 1e-31; 

const double HE_MASS = 4.00260325413;
const double AR_MASS = 39.9623831237; 
const double MU_SI = HE_MASS * AR_MASS / ( HE_MASS + AR_MASS ) * DALTON_UNIT; 
const double MU = MU_SI / AMU;

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

	double dipz = ar_he_dip_buryak_fit( R );
	Eigen::Vector3d mol_dipole( 0, 0, dipz );
	Eigen::Vector3d lab_dipole = S * mol_dipole;

	res[0] = lab_dipole[0];
	res[1] = lab_dipole[1];
	res[2] = lab_dipole[2];
}

void transform_coordinates( std::tuple<double, double, double> &he_coords, 
							std::tuple<double, double, double> &ar_coords, 
							const double &R, const double &theta )
{
	Eigen::Vector3d he_coords_mol( 0.0, 0.0, -AR_MASS / (HE_MASS + AR_MASS) * R );
	Eigen::Vector3d ar_coords_mol( 0.0, 0.0, HE_MASS / (HE_MASS + AR_MASS) * R );

	Eigen::Matrix<double, 3, 3> S;
	S(0, 0) = 1.0;
	S(0, 1) = 0.0;
	S(0, 2) = 0.0;

	S(1, 0) = 0.0;
	S(1, 1) = cos( theta );
	S(1, 2) = -sin( theta );

	S(2, 0) = 0.0;
	S(2, 1) = sin( theta );
	S(2, 2) = cos( theta );

	Eigen::Vector3d he_coords_lab = S * he_coords_mol;
	Eigen::Vector3d ar_coords_lab = S * ar_coords_mol;

	std::get<0>( he_coords ) = he_coords_lab[0];
	std::get<1>( he_coords ) = he_coords_lab[1];
	std::get<2>( he_coords ) = he_coords_lab[2];

	std::get<0>( ar_coords ) = ar_coords_lab[0];
	std::get<1>( ar_coords ) = ar_coords_lab[1];
	std::get<2>( ar_coords ) = ar_coords_lab[2];
}


void rhs( double* out, double R, double pR, double theta, double pTheta )
{
	out[0] = pR / MU; // dot{R} 
	out[1] = pow(pTheta, 2) / MU / pow(R, 3) - ar_he_pot_derivative( R ); // dot{pR}
	out[2] = pTheta / MU / pow(R, 2); // dot{theta}
	out[3] = 0; // dot[pTheta}
}
