#include "ar_he_pes.h"
#include "ar_he_pes_derivative.h"
#include "ar_he_dip.h"

#include <iostream>
#include <vector>
#include <string>

using std::cout;
using std::endl;

using std::vector;
using std::string;

void test_dipole( void )
{
	double r_step = 0.5;
	double r0 = 3.0;

	double r, dip;

	for ( int i = 0; i < 20; i++ )
	{
		r = r0 + i * r_step;
		dip = ar_he_dip( r ); 
		cout << "r: " << r << "; dip: " << dip << endl;
	}
}

void test_potential( void )
{
	double r_step = 0.25;
	double r0 = 4.25;

	double r, pot;

	for ( int i = 0; i < 100; i++ )
	{
		r = r0 + i * r_step;
		pot = ar_he_pot( r );
		cout << "r: " << r << "; pot: " << pot << " E_h" << endl;
	}
}

inline void dampfun( const double x, double* D, const int n )
{
	D[0] = 1.0;
	double exp_mx = exp( -x );
	double xk = 1.0;

	for ( int k = 1; k <= n; k++ )
	{
		xk *= x / k;
		D[k] = D[k - 1] + xk;
		D[k - 1] = 1.0 - exp_mx * D[k - 1];
	}
	
	D[n] = 1 - exp_mx * D[n];
}

int main( int argc, char* argv[] )
{
	if ( argc != 2 )
	{
		cout << "USAGE: ./.. (string) type of test" << endl;
		cout << "Types: potential dipole" << endl;
		exit( 1 );
	}

	string type = argv[1];
	cout << "type: " << type << endl;

	if ( type == "dipole" )
	{
		test_dipole();
	}

	if ( type == "potential" )
	{
		test_potential();
	}

	return 0;
}

