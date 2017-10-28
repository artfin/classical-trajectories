#include "matrix.h"

#include <iostream>

#include "ar_he_pes.h"
#include "ar_he_pes_derivative.h"
#include "ar_he_dip.h"

const double HE_MASS = 4.00260325413;
const double AR_MASS = 39.9623831237; 
const double PROTON_TO_ELECTRON_RATIO = 1836.15267389; 

const double MU = HE_MASS * AR_MASS / ( HE_MASS + AR_MASS ) * PROTON_TO_ELECTRON_RATIO; 

double ham_value( double R, double pR, double J )
{
	double kin_term = pow(pR, 2) / (2 * MU);
	double ang_term = pow(J, 2) / (2 * MU * pow(R, 2));

	return ang_term + kin_term + ar_he_pot(R); 
}


void hamiltonian( double* out, double R, double pR, double J, bool dip_calc )
{
    if ( dip_calc == false )
    {
		out[0] = - pow(J, 2) / (MU * pow(R, 3) );
		out[1] = pR / MU;
    }
    	
    else
    {
		out[0] = ar_he_dip( R );
    }
}

void rhs( double* out, double R, double pR, double J )
{
    double* derivatives = new double[2];
    hamiltonian(derivatives, R, pR, J, false);

	std::cout << "inside rhs" << std::endl;

	std::cout << "derivatives[0]: " << derivatives[0] << std::endl;
	std::cout << "derivatives[1]: " << derivatives[1] << std::endl;

    out[0] = derivatives[1]; // /dot(R) = dH/dpR
	std::cout << "out[0]: " << out[0] << std::endl;

	out[1] = - derivatives[0] - ar_he_pot_derivative(R);  // /dot(pR) = - dH/dR = - dT/dR - dU/dR (to hartrees from cm^-1)
	std::cout << "out[1]: " << out[1] << std::endl;

    delete [] derivatives;

	std::cout << "finished rhs" << std::endl;
}
