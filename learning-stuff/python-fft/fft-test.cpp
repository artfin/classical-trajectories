#include <fftw3.h>
#include <math.h>
#include <iostream>

using namespace std;

#define REAL_PART 0
#define IMAG_PART 1

int main() 
{
	const int N = 10;
	const int NC = ( N / 2 ) + 1;

	fftw_plan plan_forward;

	double *in = (double*) fftw_malloc( sizeof( double ) * N );	
	fftw_complex *out = (fftw_complex*) fftw_malloc( sizeof( fftw_complex ) * NC );

	double power_spectrum[(N + 1)/2];

	for ( int i = 0; i < N; i++ )
	{
		in[i] =  sin ( 2.0 * M_PI * i / N );
		cout << "in[" << i << "] = " << in[i] << endl;
	}	

	plan_forward = fftw_plan_dft_r2c_1d( N, in, out, FFTW_BACKWARD && FFTW_ESTIMATE); 

	fftw_execute( plan_forward );

	cout << "plan executed" << endl;

	for ( int i = 0; i < NC; i++ )
	{
		cout << "out[" << i << "] = " << out[i][0] << " + i * " << out[i][1] << endl;
	}

	power_spectrum[0] = out[0][REAL_PART] * out[0][IMAG_PART];
	for ( int k = 1; k < (N + 1) / 2; k++ )
	{
			power_spectrum[k] = out[k][REAL_PART] * out[k][REAL_PART] + out[k][IMAG_PART] * out[k][IMAG_PART] + 
								out[N - k][REAL_PART] * out[N - k][REAL_PART] + out[N -k][IMAG_PART] * out[N - k][IMAG_PART];
	}
	if ( N % 2 == 0 )
	{
			power_spectrum[N/2] = out[N/2][REAL_PART] * out[N/2][REAL_PART] + out[N/2][IMAG_PART] * out[N/2][IMAG_PART]; // nyquist frequency
	}
	for ( int i = 0; i < (N + 1) / 2; i++ )
	{
		cout << "power_spectrum[" << i << "] = " << power_spectrum[i] << endl;
	}
	
	fftw_destroy_plan( plan_forward );
	fftw_cleanup ();

	return 0;
}


