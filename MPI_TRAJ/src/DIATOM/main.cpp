#include <mpi.h>

#include <iostream>

// matrix multiplication
#include "matrix.h"

// should be included BEFORE Gear header files
// due to namespace overlap
#include <vector>

// miscellaneous functions
#include "fft.h"
#include "gsl_custom.h"
#include "aux.h"

// physical constants
#include "constants.h"

// Gear header files
#include "basis.h"
#include "vmblock.h"
#include "gear.h"
#include "t_dgls.h"

// library for FFT
#include <fftw3.h>

// gsl histogram
#include <gsl/gsl_histogram.h>

// ############################################
// Exit tag for killing slave
const int EXIT_TAG = 42;
// ############################################

// ############################################
// Number of initial conditions per trajectory
const int ICPERTRAJ = 5; 
// ############################################

// ############################################
const double Temperature = 298; 
// ############################################

// ############################################
// GSL Histogram parameters (spectrum boundaries)
const int NBINS = 50;
const double LBOUND = 0.0;
const double RBOUND = 400.0;
// ############################################

const double HE_MASS = 4.00260325413;
const double AR_MASS = 39.9623831237; 
const double PROTON_TO_ELECTRON_RATIO = 1836.15267389; 

const double MU = HE_MASS * AR_MASS / ( HE_MASS + AR_MASS ) * PROTON_TO_ELECTRON_RATIO; 

const double B_STEP = 0.100;
const double V_STEP = 1.0;

using namespace std;

void show( string name, vector<double> v, double sampling_time )
{
	cout << "#########################" << endl;
	
	for ( int i = 0; i < v.size(); i++ )
	{
		cout << "t = " << sampling_time * i << "; " << name << "[" << i << "] = " << v[i] << endl;
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

void syst (REAL t, REAL *y, REAL *f)
{
  	(void)(t); // avoid unused parameter warning 

	double *out = new double[4];

	rhs( out, y[0], y[1], y[2], y[3] );

	f[0] = out[0]; // \dot{R} 
	f[1] = out[1]; // \dot{p_R}
	f[2] = out[2]; // \dot{\theta}
	f[3] = out[3]; // \dot{p_\theta}

	delete [] out;
}

void master_code( int world_size )
{
	MPI_Status status;
	int source;

	//FILE* inputfile = fopen("input/DIATOM/weights_test", "r" );
	FILE* inputfile = fopen("input/DIATOM/test2", "r" );

	string spectrum_filename = "weights_test";

	// counter of calculated trajectories
	int NTRAJ = 0;

	double *ics = new double[ICPERTRAJ];

	// result of reading values from file
	int scanfResult;

	// sending first message to slaves
	for ( int i = 1; i < world_size; i++ ) 
	{
		scanfResult = fwscanf( inputfile, L"%lf %lf %lf %lf\n", &ics[0], &ics[1], &ics[2], &ics[3] );
		cout << "Read initial condition: " << ics[0] << " " 
										   << ics[1] << " "
										   << ics[2] << " "
										   << ics[3] << endl;
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
			save_histogram( histogram, NBINS, spectrum_filename );
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
			gsl_histogram_accumulate( histogram, freqs_package[i], intensities_package[i] * B_STEP * V_STEP );
		}

		// reading another line from file
		scanfResult = fwscanf(inputfile, L"%lf %lf %lf %lf\n", &ics[0], &ics[1], &ics[2], &ics[3], &ics[4] );	
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

void copy_initial_conditions( double* ics, REAL* y0, const int length )
{
	for ( int i = 0; i < length; i++ )
	{
		y0[i] = ics[i];
	}
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

		cout << "Received conditions: " << endl;
		cout << "ics[0] = " << ics[0] << endl;
		cout << "ics[1] = " << ics[1] << endl;
		cout << "ics[2] = " << ics[2] << endl;
		cout << "ics[3] = " << ics[3] << endl;
		cout << "########################" << endl;

		int trajectory_number = ics[0];
		ics[0] = ics[1];

		double b = ics[2];
		double v0 = ics[3];
		cout << "b = " << b << endl;
		cout << "v0 = " << v0 << endl;

		double b_weight = - 2 * b * v0;
		double v_weight = pow(v0, 2) * exp( - MU * pow(v0, 2) * constants::HTOJ / (2 * constants::BOLTZCONST * Temperature));


		ics[1] = MU * v0; // pR
		ics[2] = 0; // theta 
		ics[3] = abs( MU * v0 * b ); // p_\theta

		cout << "ics[0](R) = " << ics[0] << endl;
		cout << "ics[1](pR) = " << ics[1] << endl;
		cout << "ics[2](theta) = " << ics[2] << endl;
		cout << "ics[3](pTheta = J) = " << ics[3] << endl;

		N = 4;
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
		const double sampling_time = 5; 

		epsabs = 1E-13;
		epsrel = 1E-13;

		t0 = 0.0;

		h = 0.1;         		// initial step size
		xend = sampling_time;   // initial right bound of integration
		fmax = 1e8;  	 		// maximal number of calls 

		// calculating initial weight of trajectory
		int counter = 0;
		double end_value = ics[0] + 0.1;

		copy_initial_conditions( ics, y0, N );
		cout << "y0[0] = " << y0[0] << endl;
		cout << "y0[1] = " << y0[1] << endl;
		cout << "y0[2] = " << y0[2] << endl;
		cout << "y0[3] = " << y0[3] << endl;

		// dipole moment in laboratory frame
		vector<double> temp( 3 );

		vector<double> dipx;
		vector<double> dipy;
		vector<double> dipz;

		while ( y0[0] < end_value ) 
		{
			fehler = gear4(&t0, xend, N, syst, y0, epsabs, epsrel, &h, fmax, &aufrufe);

			if ( fehler != 0 ) 
			{
				printf(" Gear4: error n° %d\n", 10 + fehler);
				break;
			}

			cout << "t0: " << t0 << "; r: " << y0[0] << endl;
			
			// y0[0] -- R
			// y0[1] -- PR
			// y0[2] -- \theta
			// y0[3] -- p_\theta
			transform_dipole( temp, y0[0], y0[2] );

			dipx.push_back( temp[0] );
			dipy.push_back( temp[1] );
			dipz.push_back( temp[2] );

			xend = sampling_time * (counter + 2);
			aufrufe = 0;  // actual number of calls

			counter++;
		}
			
		// length of dipole vector = number of samples
		int npoints = dipz.size(); 
		int freqs_size = npoints / 2.0;		

		if ( freqs_size != 0 )
		{
			// given the number of points and sampling time we can calculate freqs vector
			vector<double> freqs = linspace( 0.0, 1.0 / ( 2.0 * sampling_time ), freqs_size );

			// due to 2pi inside Fourier transofrm
			multiply_vector( freqs, 2 * M_PI ); 
			// transforming reverse atomic time units to cm^-1
			multiply_vector( freqs, constants::HZTOCM / constants::ATU );

			cout << ">> Processing " << trajectory_number << " trajectory. npoints = " << npoints << endl;

			vector<double> dipx_out = fft( dipx );	
			vector<double> dipy_out = fft( dipy );	
			vector<double> dipz_out = fft( dipz );	

			// auxiliary variables to store interim variables
			double power;
			vector<double> intensities; 
			
			//show( "dipx_out", dipx_out );
			//show( "dipy_out", dipy_out );
			//show( "dipz_out", dipz_out );

			//cout << "b_weight: " << b_weight << endl;
			//cout << "v_weight: " << v_weight << endl;

			for ( int i = 0; i < freqs_size; i++ )
			{
				power = dipx_out[i] + dipy_out[i] + dipz_out[i];	
				intensities.push_back( power * b_weight * v_weight );	
			}

			//show( "freqs", "ints", freqs, intensities );

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
	//Initialize the MPI environment
	MPI_Init( &argc, &argv );

	//getting id of the current process
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	//getting number of running processes
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
