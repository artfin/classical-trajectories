#include <iostream>
#include <random>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

//a Mersenne Twister pseudo−random generator of 32−bit numbers with a state size of 19937 bits 
random_device rd;
mt19937 eng( rd() );
uniform_real_distribution<double> distr( 0.0, 1.0);

void rotationMatrix( Matrix3d &m )
{
	double theta = 2 * M_PI * distr( eng ); // a rotation about the pole
	double phi = 2 * M_PI * distr( eng ); // a direction to deflect the pole 
	double z = distr( eng ); // the amount of pole deflection	 
	double sz = sqrt( z );

	// a vector to perform the reflection
	Vector3d v ( cos(phi) * sz, sin(phi) * sz, sqrt(1 - z) );
	// the Householder matrix
	Matrix3d s = 2 * v * v.transpose() - Matrix<double, 3, 3>::Identity();

	Matrix3d r;
   	r << cos(theta), sin(theta), 0,
		-sin(theta), cos(theta), 0,
		0, 0, 1;
	
	m = s * r;
}

int main( int argc, char* argv[] )
{
	int n = atoi( argv[1] );
	
	// initial vector
	Vector3d v ( 0.0, 0.0, 1.0 );
	// resulting vector
	Vector3d r;

	for ( int i = 0; i < n; i++ )
	{
		// filling rotation matrix
    	Matrix3d m;
    	rotationMatrix( m );

		// performing a random rotation of OZ-vector
		r = m * v;
		
		// displaying the components of resulting vector
		cout << r(0) << " " << r(1) << " " << r(2) << endl;
	}

    return 0;
}
