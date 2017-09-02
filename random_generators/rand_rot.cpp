#include <iostream>
#include <random>

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

random_device rd;
mt19937 eng( rd() );
uniform_real_distribution<double> distr( 0.0, 1.0);

void RMatrix( Matrix3d &m )
{
	// pick a rotation about the pole
    double theta = 2 * M_PI * distr( eng );
	// pick a direction to deflect the pole
	double phi = 2 * M_PI * distr( eng );
	// pick the amount of pole deflection
	double z = distr( eng );
	double sz = sqrt( z );

	Vector3d v ( cos(phi) * sz, sin(phi) * sz, sqrt(1-z) );
	Matrix3d s = 2 * v * v.transpose() - Matrix<double, 3, 3>::Identity();

	Matrix3d r;
	r(0, 0) = cos(theta);
	r(0, 1) = sin(theta);
	r(0, 2) = 0;
	r(1, 0) = -sin(theta);
	r(1, 1) = cos(theta);
	r(1, 2) = 0;
	r(2, 0) = 0;
	r(2, 1) = 0;
	r(2, 2) = 1;
	
	m = s * r;
}

int main( int argc, char* argv[] )
{
	int n = atoi( argv[1] );

	Vector3d v (0.0, 0.0, 1.0 );
	Vector3d r;

	for ( int i = 0; i < n; i++ )
	{
    	Matrix3d m;
    	RMatrix( m );

		r = m * v;

		cout << r(0) << " " << r(1) << " " << r(2) << endl;
	}

    return 0;
}
