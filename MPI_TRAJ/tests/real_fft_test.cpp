#include "fft.h"

#include <cmath>
#include <vector>
#include <string>

#include <iostream>

using std::cout;
using std::endl;

using std::vector;
using std::string;

void sample_sin( vector<double> &input, const double sampling_time, const int N )
{
	double x0 = 0.0;
	double x;	

	for ( int i = 0; i < N; i++ )
	{
		x = x0 + i * sampling_time; 
		
		input.push_back( 7 * sin( 2 * M_PI * x ) + 3 * sin( 2 * M_PI * 5 * x ) + 5 * sin(2 * M_PI * x * 0.5) + sin( 2 * M_PI * x * 0.1 ) ); 
	}
}

int main( int argc, char* argv[] )
{
	vector<double> input;

	int N = 100;
	double sampling_time = 0.1;

	int freq_points = std::ceil( (N + 1) / 2.0 );

	sample_sin( input, sampling_time, 30 );

	Fourier fourier( N );

	copy_to( input, fourier.inx );

	fourier.do_fourier();
	
	//for ( int i = 0; i < N; i++ )
	//{
		//cout << (*fourier.outx)[i] << endl;
	//}

	return 0;
}
