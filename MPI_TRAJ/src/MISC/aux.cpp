#include "aux.h"

// similar to python's linspace
std::vector<double> linspace( const double min, const double max, const int npoints )
{
	const double step = ( max - min ) / ( npoints - 1 );
	
	std::vector<double> res;
	for ( double temp = min; temp <= max; temp += step )
	{
		res.push_back( temp );
	}

	return res;
}

// multiply vector by given factor
void multiply_vector( std::vector<double> &v, const double factor )
{
	for ( int i = 0; i < v.size(); i++ )
	{
		v[i] *= factor;
	}
}

// copy from vector to fftw_complex*
void copy_to( std::vector<double> &v, fftw_complex* arr )
{
	for ( int i = 0; i < v.size(); i++ )
	{
		arr[i][0] = v[i];
		arr[i][1] = 0;
	}
}


