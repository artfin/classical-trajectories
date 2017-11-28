#include "fft.h"

void Fourier::zero_out_input( void )
{
	for ( int i = 0; i < this->MaxTrajectoryLength; i++ )
	{
		this->inx[i] = 0.0;
		this->iny[i] = 0.0;
		this->inz[i] = 0.0;
	}
}

Fourier::Fourier( int MaxTrajectoryLength )
{
	this->MaxTrajectoryLength = MaxTrajectoryLength;
	
	int outlen = MaxTrajectoryLength / 2 + 1;

	this->inx = (double*) fftw_malloc( sizeof(double) * MaxTrajectoryLength );
	this->iny = (double*) fftw_malloc( sizeof(double) * MaxTrajectoryLength );
	this->inz = (double*) fftw_malloc( sizeof(double) * MaxTrajectoryLength );

	this->outx = (fftw_complex*) fftw_malloc( sizeof(fftw_complex) * outlen ); 
	this->outy = (fftw_complex*) fftw_malloc( sizeof(fftw_complex) * outlen );
	this->outz = (fftw_complex*) fftw_malloc( sizeof(fftw_complex) * outlen ); 

	this->px = fftw_plan_dft_r2c_1d( MaxTrajectoryLength, this->inx, this->outx, FFTW_ESTIMATE );
	this->py = fftw_plan_dft_r2c_1d( MaxTrajectoryLength, this->iny, this->outy, FFTW_ESTIMATE );
	this->pz = fftw_plan_dft_r2c_1d( MaxTrajectoryLength, this->inz, this->outz, FFTW_ESTIMATE );

	zero_out_input( );	
}

Fourier::~Fourier()
{
	//std::cout << "Doing fourier cleanup" << std::endl;

	fftw_destroy_plan( this->px );
	fftw_destroy_plan( this->py );
	fftw_destroy_plan( this->pz );

	//std::cout << "Destroyed plans" << std::endl;

	fftw_free( this->inx );
	fftw_free( this->iny );
	fftw_free( this->inz );

	fftw_free( this->outx );
	fftw_free( this->outy );
	fftw_free( this->outz );
}

void Fourier::do_fourier( void )
{
	//for ( int i = 0; i < 10; i++ )
	//{
		//if ( inz[i] != 0 )
		//{
			//std::cout << "inz[" << i << "] = " << inz[i] << std::endl;  
		//}
	//}
	
	fftw_execute( this->px );
	fftw_execute( this->py );
	fftw_execute( this->pz );
}


//void fft_positive( std::vector<double> &signal )
//{
	//int npoints = signal.size();
	//int freq_points = ceil( npoints / 2.0 );

	//double* in = (double*) fftw_malloc( sizeof(double) * npoints );
	//fftw_complex* out = (fftw_complex*) fftw_malloc( sizeof(fftw_complex) * freq_points );

	//fftw_plan p = fftw_plan_dft_r2c_1d( npoints, in, out, FFTW_ESTIMATE );

	//for ( int i = 0; i < npoints; i++ )
	//{
		//in[i] = signal[i];
	//}

	//// performing fft
	//fftw_execute( p );

	//std::ofstream output = std::fopen( "spectrum.txt", "w" );

	//double realVal;
	//double imagVal;
	//double powVal;
	//double absVal;

	//for ( int i = 0; i < freq_points; i++ )
	//{
		//realVal = out[i][0] / npoints;
		//imagVal = out[i][1] / npoints;

		//powVal = 2 * ( realVal * realVal + imagVal * imagVal );
		//absVal = sqrt( powVal / 2.0 );

		//if ( i == 0 )
		//{
			//powVal /= 2;
		//}

		//output << realVal << " " <<
				  //imagVal << " " <<
				  //powVal << " " <<
				  //absVal << std::endl;
	//}

	//fclose( output );

	//fftw_destroy_plan( p );
	//fftw_free( in );
	//fftw_free( out );
//}

std::vector<double> fft_one_side( std::vector<double> signal )
{
	int npoints = signal.size();
	int freq_points = npoints / 2 + 1;

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

	for ( int i = 0; i < freq_points; i++ )
	{
		power = sqrt( 
					  signal_fd[i][REALPART] * signal_fd[i][REALPART] +
					  signal_fd[i][IMAGPART] * signal_fd[i][IMAGPART]
					);
		ints.push_back( power );
	}

	fftw_free( signal_td );
	fftw_free( signal_fd );

	return ints;
}

std::vector<double> fft_two_side( std::vector<double> &signal )
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

	for ( int i = 0; i < npoints; i++ )
	{
		power = sqrt( 
						signal_fd[i][REALPART] * signal_fd[i][REALPART] +
						signal_fd[i][IMAGPART] * signal_fd[i][IMAGPART]
					);
		ints.push_back( power );
	}

	fftw_free( signal_td );
	fftw_free( signal_fd );

	return ints;
}

void fftfreq( std::vector<double> &freqs, const int len, const double d )
{
	if ( len % 2 == 0 )
	{
		int i;
		for ( i = 0; i < ( len / 2 ); i++ )
		{
			freqs.push_back( i / ( len * d ) );
		}
		for ( i = - len / 2; i < 0; i++ )
		{
			freqs.push_back( i / ( len * d ) );
		}
	}
	else
	{
		int i;
		for ( i = 0; i < ( len + 1 ) / 2; i++ )
		{
			freqs.push_back( i / ( len * d ) );
		}
		for ( i = - ( len - 1 ) / 2; i < 0; i++ )
		{
			freqs.push_back( i / ( len * d ) );
		}
	}
}

void fftshift( std::vector<double> &arr, const int len )
{
	int center = ( len + 1 ) / 2.0;
	std::rotate( arr.begin(), arr.begin() + center, arr.end() );
}

void ifftshift( std::vector<double> &arr, const int len )
{
	int center = len - ( len + 1 ) / 2.0;
	std::rotate( arr.begin(), arr.begin() + center, arr.end() );
}
