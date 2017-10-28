#include "fft.h"

std::vector<double> fft( std::vector<double> signal )
{
	int npoints = signal.size();

	// td == time domain
	fftw_complex *signal_td = (fftw_complex*) fftw_malloc( sizeof(fftw_complex) * npoints );
	// fd == frequency domain
	fftw_complex *signal_fd = (fftw_complex*) fftw_malloc( sizeof(fftw_complex) * npoints );

	// copying from vector<double> to fftw_complex*
	copy_to( signal, signal_td );

	// +1 -- backward Fourier transform
	fftw_plan plan = fftw_plan_dft_1d( npoints, signal_td, signal_fd, +1, FFTW_ESTIMATE );

	fftw_execute( plan );

	fftw_destroy_plan( plan );

	double power;
	std::vector<double> ints;

	for ( int i = 0; i < npoints / 2; i++ )
	{
		power = signal_fd[i][REALPART] * signal_fd[i][REALPART] +
				signal_fd[i][IMAGPART] * signal_fd[i][IMAGPART];
		
		ints.push_back( power );
	}

	fftw_free( signal_td );
	fftw_free( signal_fd );

	return ints;
}


