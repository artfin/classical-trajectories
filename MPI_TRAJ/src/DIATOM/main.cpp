#include <mpi.h>

// matrix multiplication
#include "matrix.h"

// should be included BEFORE Gear header files
// due to namespace overlap
#include <vector>

// to perform double to string conversion
#include <sstream>
#include <string>

#include <fstream>

#include <iostream>
#include <stdlib.h>
#include <stdio.h>

// Gear header files
#include "basis.h"
#include "vmblock.h"
#include "gear.h"
#include "t_dgls.h"

// library for FFT
#include <fftw3.h>

// gsl histogram
#include <gsl/gsl_histogram.h>

// macros for real and imaginary parts
#define REALPART 0
#define IMAGPART 1

const int EXIT_TAG = 42; // exiting tag
const int ICPERTRAJ = 8; // number of initial conditions per trajectory

// cm^-1 to Hz
//const double CMTOHZ = 2.99793 * pow(10, 10);

// hartree to joules
const double HTOJ = 4.35974417 * pow(10, -18);
// boltzmann constant
const long double BOLTZCONST = 1.38064852 * pow(10, -23);

const double ATU = 2.418884326505 * pow( 10, -17 );
const double HZTOCM = 3.335631 * pow( 10, -11 );

const double Temperature = 298; // K

const int NBINS = 306;
const double LBOUND = 0.0;
const double RBOUND = 306.0;

using namespace std;

void syst (REAL t, REAL *y, REAL *f)
{
  (void)(t); // avoid unused parameter warning 

  double *out = new double[6];
  rhs(out, y[0], y[1], y[2], y[3], y[4], y[5], y[6]);
  // R  Theta pR pT phi theta J

  f[0] = out[0]; // dR/dt  
  f[1] = out[1]; // d(Theta)/dt
  f[2] = out[2]; // d(pR)/dt
  f[3] = out[3]; // d(pT)/dt
  f[4] = out[4]; // d(phi)/dt
  f[5] = out[5]; // d(theta)/dt

  delete [] out;
}

void save_histogram( gsl_histogram *histogram )
{
	ofstream file( "test" );

	double lower_bound, higher_bound, bin_content;
	for ( int counter = 0; counter < NBINS; counter++ )
	{
		gsl_histogram_get_range( histogram, counter, &lower_bound, &higher_bound );
		bin_content = gsl_histogram_get( histogram, counter );

		file << lower_bound << " " << higher_bound << " " << bin_content << endl;
	}
	
	file.close();	
}

// similar to python's linspace
vector<double> linspace( double min, double max, double npoints )
{
	double step = ( max - min ) / ( npoints - 1 );
	
	vector<double> res;
	for ( double temp = min; temp <= max; temp += step )
	{
		res.push_back( temp );
	}

	return res;
}

// multiply vector by given factor
void multiply_vector( vector<double> &v, double factor )
{
	for ( int i = 0; i < v.size(); i++ )
	{
		v[i] *= factor;
	}
}

void master_code( int world_size )
{
	MPI_Status status;
	int source;

	//FILE* inputfile = fopen("input/test", "r");
	FILE* inputfile = fopen("input/ics.txt", "r" );

	// counter of calculated trajectories
	int NTRAJ = 0;

	double *ics = new double[ICPERTRAJ];

	// result of reading values from file
	int scanfResult;

	// sending first message to slaves
	for ( int i = 1; i < world_size; i++ ) 
	{
		scanfResult = fwscanf(inputfile, L"%lf %lf %lf %lf %lf %lf %lf %lf\n", &ics[0], &ics[1], &ics[2], &ics[3], &ics[4], &ics[5], &ics[6], &ics[7]);
		MPI_Send(&ics[0], ICPERTRAJ, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);

		NTRAJ++;
	}
	
	// number of alive processes
	int alive = world_size - 1;

	// prepare bins
	gsl_histogram *histogram = gsl_histogram_alloc( NBINS );
	gsl_histogram_set_ranges_uniform( histogram, LBOUND, RBOUND );

	while ( true )
	{	
		// exit when all slaves are killed
		if ( alive == 0 )
		{
			break;
		}
	
		if ( NTRAJ % 100 == 0 )
		{
			cout << ">> Saving histogram... " << endl;
			save_histogram( histogram );
		}

		// receiving message from any of slaves	
		int package_size;
		MPI_Recv( &package_size, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
		source = status.MPI_SOURCE;
	
		//cout << "Master received a package size = " << package_size << " from process " << source << endl;

		vector<double> freqs_package( package_size );
		vector<double> intensities_package( package_size );

		if ( package_size != 0 )
		{
			MPI_Recv( &freqs_package[0], package_size, MPI_DOUBLE, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
			//cout << "Master received frequency package from process " << source << endl;
			
			MPI_Recv( &intensities_package[0], package_size, MPI_DOUBLE, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
			//cout << "Master received intensities_package from process " << source << endl;
		}
		
		// applying immediate binning of values
		for ( int i = 0; i < freqs_package.size(); i++ )
		{
			// increases the value of appropriate bin by intensity
			gsl_histogram_accumulate( histogram, freqs_package[i], intensities_package[i] );
		}

		// reading another line from file
		scanfResult = fwscanf(inputfile, L"%lf %lf %lf %lf %lf %lf %lf %lf\n", &ics[0], &ics[1], &ics[2], &ics[3], &ics[4], &ics[5], &ics[6], &ics[7]);	
		// if it's not ended yet then sending new chunk of work to slave
		if ( scanfResult != -1)
   		{
			//cout << "Master sends new initial conditions to process " << source << endl;
			MPI_Send(&ics[0], ICPERTRAJ, MPI_DOUBLE, source, 0, MPI_COMM_WORLD);
			NTRAJ++;
		}
		// work is done, sending a killing message
		else
		{
			MPI_Send(&ics[0], ICPERTRAJ, MPI_DOUBLE, source, EXIT_TAG, MPI_COMM_WORLD);
			alive--;
		}
	}

	gsl_histogram_free( histogram );

	fftw_cleanup();

	delete [] ics;	

	fclose(inputfile);
}

void copy_to( vector<double> &v, fftw_complex* arr )
{
	for ( int i = 0; i < v.size(); i++ )
	{
		arr[i][0] = v[i];
		arr[i][1] = 0;
	}
}

vector<double> fft( vector<double> signal )
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
	vector<double> ints;

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

void show( string name, vector<double> v )
{
	cout << "#########################" << endl;
	
	for ( int i = 0; i < v.size(); i++ )
	{
		cout << name << "[" << i << "] = " << v[i] << endl;
	}
	
	cout << "#########################" << endl;
}

void show( string name1, string name2, vector<double> v1, vector<double> v2 )
{
	cout << "#########################" << endl;

	for ( int i = 0; i < v1.size(); i++ )
	{
		cout << name1 << "[" << i << "] = " << v1[i] << "; " <<
				name2 << "[" << i << "] = " << v2[i] << endl;
	}	
	
	cout << "#########################" << endl;
}

void slave_code( int world_rank )
{
	REAL epsabs;    //  absolute error bound
	REAL epsrel;    //  relative error bound    
	REAL t0;        // left edge of integration interval
	REAL *y0;       // [0..n-1]-vector: initial value, approxim. 
	REAL h;         // initial, final step size
	REAL xend;      // right edge of integration interval 

	long fmax;      // maximal number of calls of right side in gear4()
	long aufrufe;   // actual number of function calls
	int  N;         // number of DEs in system
	int  fehler;    // error code from umleiten(), gear4()
	int  i;         // loop counter
								 
	void *vmblock;  // List of dynamically allocated vectors
	
	// array to store initial conditions
	double *ics = new double [ICPERTRAJ];

	// auxiliary variable to store status of message
	MPI_Status status;
		
	while ( true )
	{
		MPI_Recv(&ics[0], ICPERTRAJ, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		//cout << ">> Process " << world_rank << " received new initial conditions." << endl; 

		if ( status.MPI_TAG == EXIT_TAG )
		{
			break;
		}

		N = 6;
		vmblock = vminit();
		y0 = (REAL*) vmalloc(vmblock, VEKTOR, N, 0);
			
		// out of memory?
  		if ( !vmcomplete(vmblock) )
		{ 
   			printf("mgear: out of memory.\n");
   			return;
 		}

		// according to Ivanov:
		// sampling time = 50 fs

		// in atomic time units
		const double sampling_time = 2250.0; 

  		epsabs = 1E-13;
  		epsrel = 1E-13;
  
  		t0 = 0.0;
			
  		h = 0.1;         		// initial step size
  		xend = sampling_time;   // initial right bound of integration
  		fmax = 1e8;  	 		// maximal number of calls 
			
  		double *dipole = new double [3];
			
  		vector<double> ddipx;
  		vector<double> ddipy;
  		vector<double> ddipz;
		
		// r, theta, pr, ptheta, phi, theta, j
		for ( int i = 0; i < 7; i++ )
		{
			y0[i] = ics[i + 1];
		}
	
		// calculating initial weight of trajectory
		// j == ics[7]; theta == ics[6]; phi == ics[5]
		double jx = ics[7] * sin(ics[6]) * cos(ics[5]);
		double jy = ics[7] * sin(ics[6]) * sin(ics[5]);
		double jz = ics[7] * cos(ics[6]);
		double h0 = ham_value( ics[1], ics[2], ics[3], ics[4], jx, jy, jz);
		double exp_hkt = exp( - h0 * HTOJ / ( BOLTZCONST * Temperature ));
		int counter = 0;
		double end_value = ics[1] + 0.1;

		while ( y0[0] < end_value ) 
		{
     		fehler = gear4(&t0, xend, N, syst, y0, epsabs, epsrel, &h, fmax, &aufrufe);
     			
     		if ( fehler != 0 ) 
			{
     			printf(" Gear4: error n° %d\n", 10 + fehler);
				break;
     		}
     
     		hamiltonian(dipole, y0[0], y0[1], y0[2], y0[3], y0[4], y0[5], y0[6], true);
     		
			// collecting derivatives of dipole in laboratory frame
			ddipx.push_back( dipole[0] );
     		ddipy.push_back( dipole[1] );
     		ddipz.push_back( dipole[2] );
				
     		xend = sampling_time * (counter + 2);
     		aufrufe = 0;  // actual number of calls
  				
			counter++;
		}
	
		// freeing dipole memory
		delete [] dipole;

		// length of dipole vector = number of samples
 	    int npoints = ddipx.size(); 
		int freqs_size = npoints / 2.0;		

		if ( freqs_size != 0 )
		{
			// given the number of points and sampling time we can calculate freqs vector
			vector<double> freqs = linspace( 0.0, 1.0 / ( 2.0 * sampling_time ), freqs_size );
		
			// deleting first frequency == 0
			freqs.erase( freqs.begin() );

			// due to 2pi inside Fourier transofrm
			multiply_vector( freqs, 2 * M_PI ); 
			// transforming reverse atomic time units to cm^-1
			multiply_vector( freqs, HZTOCM / ATU );

			cout << ">> Processing " << ics[0] << " trajectory. npoints = " << npoints << endl;

			vector<double> ddipx_out = fft( ddipx );
	   		vector<double> ddipy_out = fft( ddipy );
			vector<double> ddipz_out = fft( ddipz );	
	
			// auxiliary variables to store interim variables
			double power;
			vector<double> intensities; 

	  		for ( int i = 0; i < freqs_size; i++ )
	  		{
				power = ddipx_out[i] + ddipy_out[i] + ddipz_out[i];
				power = power / pow( freqs[i], 2 ) * exp_hkt ;	
				
				intensities.push_back( power );	
			}
		
			show( "freqs", "ints", freqs, intensities );

			MPI_Send( &freqs_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD );
			//cout << ">> Process " << world_rank << " sends package size = " << freqs_size << " to root." << endl;

			MPI_Send( &freqs[0], freqs_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
			//cout << ">> Process " << world_rank << " sends frequency package to root" << endl;

			MPI_Send( &intensities[0], freqs_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
			//cout << ">> Process " << world_rank << " sends intensities package to root" << endl;
		}
		else
		{
			MPI_Send( &freqs_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD );
			//cout << ">> Process " << world_rank << " sends package size = " << freqs_size << " to root." << endl;
		}
	}
}

int main( int argc, char* argv[] )
{
	// Initialize the MPI environment
	MPI_Init( &argc, &argv );

	// getting id of the current process
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	// getting number of running processes
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size); 
	
	if ( world_rank == 0 ) 
	{
		master_code( world_size );
	}	
	else
	{
		slave_code( world_rank );
	}

	MPI_Finalize();
	
	return 0;
}
